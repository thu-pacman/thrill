/*******************************************************************************
 * thrill/data/block_pool.cpp
 *
 * Part of Project Thrill - http://project-thrill.org
 *
 * Copyright (C) 2015-2016 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <thrill/common/logger.hpp>
#include <thrill/common/math.hpp>
#include <thrill/common/lowlevel_cache.h>
#include <thrill/data/block.hpp>
#include <thrill/data/block_pool.hpp>
#include <thrill/mem/aligned_allocator.hpp>
#include <thrill/mem/pool.hpp>

#include <foxxll/io/file.hpp>
#include <foxxll/io/iostats.hpp>
#include <tlx/container/lru_cache.hpp>
#include <tlx/die.hpp>
#include <tlx/math/is_power_of_two.hpp>
#include <tlx/string/join_generic.hpp>

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace thrill {
namespace data {

//! debug block life cycle output: create, destroy
static constexpr bool debug_blc = false;

//! debug block memory alloc and dealloc
static constexpr bool debug_alloc = false;

//! debug block pinning:
static constexpr bool debug_pin = false;

//! debug memory requests
static constexpr bool debug_mem = false;

//! debug block eviction: evict, write complete, read complete
static constexpr bool debug_em = false;

/******************************************************************************/
// std::new_handler() which gets called when malloc() returns nullptr

static std::recursive_mutex s_new_mutex;
static std::vector<BlockPool*> s_blockpools;

static std::atomic<bool> in_new_handler {
    false
};

static void OurNewHandler() {
    std::unique_lock<std::recursive_mutex> lock(s_new_mutex);
    if (in_new_handler) {
        printf("new handler called recursively! fixup using mem::Pool!\n");
        abort();
    }

    static bool s_notify_new_handler = false;
    if (!s_notify_new_handler) {
        fprintf(stderr, "Thrill: new_handler called! Program is out of C++ heap memory, trying to\n");
        fprintf(stderr, "Thrill: swap out Blocks to external memory. Check your program's memory usage.\n");
        in_new_handler = true;
    }

    static size_t s_iter = 0;
    foxxll::request_ptr req;

    // first try to find a handle to a currently being written block.
    for (size_t i = 0; i < s_blockpools.size(); ++i) {
        req = s_blockpools[s_iter]->GetAnyWriting();
        ++s_iter %= s_blockpools.size();
        if (req) break;
    }

    if (!req) {
        // if no writing active, evict a block
        for (size_t i = 0; i < s_blockpools.size(); ++i) {
            req = s_blockpools[s_iter]->EvictBlockLRU();
            ++s_iter %= s_blockpools.size();
            if (req) break;
        }
    }

    if (req) {
        req->wait();
        in_new_handler = false;
    }
    else {
        printf("new handler found no ByteBlock to evict.\n");
        for (size_t i = 0; i < s_blockpools.size(); ++i) {
            LOG1 << "BlockPool[" << i << "]"
                 << " total_blocks=" << s_blockpools[i]->total_blocks()
                 << " total_bytes=" << s_blockpools[i]->total_bytes()
                 << " pinned_blocks=" << s_blockpools[i]->pinned_blocks()
                 << " writing_blocks=" << s_blockpools[i]->writing_blocks()
                 << " swapped_blocks=" << s_blockpools[i]->swapped_blocks()
                 << " reading_blocks=" << s_blockpools[i]->reading_blocks();
        }
        mem::malloc_tracker_print_status();
        in_new_handler = false;

        lock.unlock();
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

/******************************************************************************/
// BlockPool::PinCount

struct BlockPool::Counter {
    //! current counter value
    size_t value = 0;
    //! maximum counter value since last stats read
    size_t hmax = 0;

    operator size_t () const { return value; }

    Counter& operator += (size_t v) {
        value += v;
        hmax = std::max(hmax, value);
        return *this;
    }
    Counter& operator -= (size_t v) {
        value -= v;
        return *this;
    }

    //! get last held max value and update to current
    size_t hmax_update() {
        size_t m = hmax;
        hmax = value;
        return m;
    }
};

struct BlockPool::PinCount {
    //! current total number of pins, where each thread pin counts
    //! individually.
    size_t              total_pins_ = 0;

    //! total number of bytes pinned.
    Counter             total_pinned_bytes_;

    //! maximum number of total pins
    size_t              max_pins = 0;

    //! maximum number of pinned bytes
    size_t              max_pinned_bytes = 0;

    //! number of pinned blocks per local worker id - this is used to count
    //! the amount of memory locked per thread.
    std::vector<size_t> pin_count_;

    //! number of bytes pinned per local worker id.
    std::vector<size_t> pinned_bytes_;

    //! ctor: initializes vectors to correct size.
    explicit PinCount(size_t workers_per_host);

    //! increment pin counter for thread_id by given size in bytes
    void Increment(size_t local_worker_id, size_t size);

    //! decrement pin counter for thread_id by given size in bytes
    void Decrement(size_t local_worker_id, size_t size);

    //! assert that it is zero.
    void AssertZero() const;
};

BlockPool::PinCount::PinCount(size_t workers_per_host)
    : pin_count_(workers_per_host),
      pinned_bytes_(workers_per_host) { }

void BlockPool::PinCount::Increment(size_t local_worker_id, size_t size) {
    ++pin_count_[local_worker_id];
    pinned_bytes_[local_worker_id] += size;
    ++total_pins_;
    total_pinned_bytes_ += size;
    max_pins = std::max(max_pins, total_pins_);
    max_pinned_bytes = std::max(max_pinned_bytes, total_pinned_bytes_.value);
}

void BlockPool::PinCount::Decrement(size_t local_worker_id, size_t size) {
    die_unless(pin_count_[local_worker_id] > 0);
    die_unless(pinned_bytes_[local_worker_id] >= size);
    die_unless(total_pins_ > 0);
    die_unless(total_pinned_bytes_ >= size);

    --pin_count_[local_worker_id];
    pinned_bytes_[local_worker_id] -= size;
    --total_pins_;
    total_pinned_bytes_ -= size;
}

void BlockPool::PinCount::AssertZero() const {
    die_unless(total_pins_ == 0);
    die_unless(total_pinned_bytes_ == 0);
    for (const size_t& pc : pin_count_)
        die_unless(pc == 0);
    for (const size_t& pb : pinned_bytes_)
        die_unless(pb == 0);
}

std::ostream& operator << (std::ostream& os, const BlockPool::PinCount& p) {
    os << " total_pins_=" << p.total_pins_
       << " total_pinned_bytes_=" << p.total_pinned_bytes_
       << " pin_count_=[" << tlx::join(',', p.pin_count_) << "]"
       << " pinned_bytes_=[" << tlx::join(',', p.pinned_bytes_) << "]"
       << " max_pin=" << p.max_pins
       << " max_pinned_bytes=" << p.max_pinned_bytes;
    return os;
}

/******************************************************************************/
// BlockPool::Data

//! type of set of ByteBlocks currently begin written to EM.
using WritingMap = std::unordered_map<
    ByteBlock*, foxxll::request_ptr,
    std::hash<ByteBlock*>, std::equal_to<>,
    mem::GPoolAllocator<std::pair<ByteBlock* const, foxxll::request_ptr> > >;

//! type of set of ByteBlocks currently begin read from EM.
using ReadingMap = std::unordered_map<
    ByteBlock*, PinRequestPtr,
    std::hash<ByteBlock*>, std::equal_to<>,
    mem::GPoolAllocator<
        std::pair<ByteBlock* const, PinRequestPtr> > >;

class BlockPool::Data
{
public:
    //! For waiting on hard memory limit
    std::condition_variable cv_memory_change_;

    //! Soft limit for the block pool, blocks will be written to disk if this
    //! limit is reached. 0 for no limit.
    size_t soft_ram_limit_;

    //! Hard limit for the block pool, memory requests will block if this limit
    //! is reached. 0 for no limit.
    size_t hard_ram_limit_;

    //! print a message on the first block evicted to external memory
    bool notify_em_used_ = false;

    //! list of all blocks that are _in_memory_ but are _not_ pinned.
    tlx::LruCacheSet<
        ByteBlock*, mem::GPoolAllocator<ByteBlock*> > unpinned_blocks_;

    //! set of ByteBlocks currently begin written to EM.
    WritingMap writing_;

    //! set of ByteBlocks currently begin read from EM.
    ReadingMap reading_;

    //! set of ByteBlock currently in EM.
    std::unordered_set<
        ByteBlock*, std::hash<ByteBlock*>, std::equal_to<>,
        mem::GPoolAllocator<ByteBlock*> > swapped_;

    //! I/O layer stats when BlockPool was created.
    foxxll::stats_data io_stats_first_;

    //! I/O layer stats of previous profile tick
    foxxll::stats_data io_stats_prev_;

    //! reference to io block manager
    foxxll::block_manager* bm_;

    //! Allocator for ByteBlocks such that they are aligned for faster
    //! I/O. Allocations are counted via mem_manager_.
    mem::AlignedAllocator<Byte, mem::Allocator<char> > aligned_alloc_;

    //! next unique File id
    std::atomic<size_t> next_file_id_ { 0 };

    //! number of unpinned bytes
    Counter unpinned_bytes_;

    //! pin counter class
    PinCount pin_count_;

    //! number of bytes currently begin requested from RAM.
    size_t requested_bytes_ = 0;

    //! number of bytes currently being written to EM.
    Counter writing_bytes_;

    //! total number of bytes in swapped blocks
    Counter swapped_bytes_;

    //! number of bytes currently being read from to EM.
    Counter reading_bytes_;

    //! total number of ByteBlocks allocated
    size_t total_byte_blocks_ = 0;

    //! condition variable to wait on for ByteBlock deallocation
    std::condition_variable cv_total_byte_blocks_;

    //! total number of bytes in all ByteBlocks (in memory or swapped)
    Counter total_bytes_;

    //! maximum number of bytes in all ByteBlocks (in memory or swapped)
    size_t max_total_bytes_ = 0;

    //! total number of bytes used in RAM by pinned and unpinned blocks, and
    //! also additionally reserved memory via BlockPoolMemoryHolder.
    Counter total_ram_bytes_;

    //! last time statistics where outputted
    std::chrono::steady_clock::time_point tp_last_
        = std::chrono::steady_clock::now();

public:
    Data(BlockPool& block_pool,
         size_t soft_ram_limit, size_t hard_ram_limit,
         size_t workers_per_host)
        : soft_ram_limit_(soft_ram_limit),
          hard_ram_limit_(hard_ram_limit),
          bm_(foxxll::block_manager::get_instance()),
          aligned_alloc_(mem::Allocator<char>(block_pool.mem_manager_)),
          pin_count_(workers_per_host) { }

    //! Updates the memory manager for internal memory. If the hard limit is
    //! reached, the call is blocked intil memory is free'd
    void IntRequestInternalMemory(std::unique_lock<std::mutex>& lock, size_t size);

    //! Updates the memory manager for the internal memory, wakes up waiting
    //! BlockPool::RequestInternalMemory calls
    void IntReleaseInternalMemory(size_t size);

    //! Unpins a block. If all pins are removed, the block might be swapped.
    //! Returns immediately. Actual unpinning is async.
    void IntUnpinBlock(
        BlockPool& bp, ByteBlock* block_ptr, size_t local_worker_id);

    //! Evict a block from the lru list into external memory
    foxxll::request_ptr IntEvictBlockLRU();

    //! Evict a block into external memory. The block must be unpinned and not
    //! swapped.
    foxxll::request_ptr IntEvictBlock(ByteBlock* block_ptr);

    //! \name Block Statistics
    //! \{

    //! Total number of allocated blocks of this block pool
    size_t int_total_blocks()  noexcept;

    //! Total number of bytes allocated in blocks of this block pool
    size_t int_total_bytes()  noexcept;

    //! \}
};

/******************************************************************************/
// BlockPool

BlockPool::BlockPool(size_t workers_per_host)
    : BlockPool(0, 0, nullptr, nullptr, workers_per_host) { }

BlockPool::BlockPool(size_t soft_ram_limit, size_t hard_ram_limit,
                     common::JsonLogger* logger, mem::Manager* mem_manager,
                     size_t workers_per_host)
    : logger_(logger),
      mem_manager_(mem_manager, "BlockPool"),
      workers_per_host_(workers_per_host),
      d_(std::make_unique<Data>(
             *this, soft_ram_limit, hard_ram_limit, workers_per_host)) {

    die_unless(hard_ram_limit >= soft_ram_limit);
    {
        std::unique_lock<std::recursive_mutex> lock(s_new_mutex);
        // register BlockPool as method of OurNewHandler to free memory.
        s_blockpools.reserve(32);
        s_blockpools.push_back(this);

        std::set_new_handler(OurNewHandler);
    }

    d_->io_stats_first_ = d_->io_stats_prev_ =
        foxxll::stats_data(*foxxll::stats::get_instance());

    logger_ << "class" << "BlockPool"
            << "event" << "create"
            << "soft_ram_limit" << soft_ram_limit
            << "hard_ram_limit" << hard_ram_limit;
}

BlockPool::~BlockPool() {
    std::unique_lock<std::mutex> lock(mutex_);

    // check that not writing any block.
    while (d_->writing_.begin() != d_->writing_.end()) {

        ByteBlock* block_ptr = d_->writing_.begin()->first;
        foxxll::request_ptr req = d_->writing_.begin()->second;

        LOGC(debug_em)
            << "BlockPool::~BlockPool() block=" << block_ptr
            << " is currently begin written to external memory, canceling.";

        lock.unlock();
        // cancel I/O request
        if (!req->cancel()) {

            LOGC(debug_em)
                << "BlockPool::~BlockPool() block=" << block_ptr
                << " is currently begin written to external memory,"
                << " cancel failed, waiting.";

            // must still wait for cancellation to complete and the I/O handler.
            req->wait();
        }
        lock.lock();

        LOGC(debug_em)
            << "BlockPool::PinBlock block=" << block_ptr
            << " is currently begin written to external memory,"
            << " cancel/wait done.";
    }

    die_unless(d_->writing_bytes_ == 0);

    // check that not reading any block.
    while (d_->reading_.begin() != d_->reading_.end()) {

        ByteBlock* block_ptr = d_->reading_.begin()->first;
        PinRequestPtr read = d_->reading_.begin()->second;

        LOGC(debug_em)
            << "BlockPool::~BlockPool() block=" << block_ptr
            << " is currently begin read from external memory, waiting.";

        lock.unlock();
        // wait for I/O request for completion and the I/O handler.
        read->req_->wait();
        lock.lock();
    }

    die_unless(d_->reading_bytes_ == 0);

    // wait for deletion of last ByteBlocks. this may actually be needed, when
    // the I/O handlers have been finished, and the corresponding references are
    // freed, but DestroyBlock() could not be called yet.
    while (d_->total_byte_blocks_ != 0)
        d_->cv_total_byte_blocks_.wait(lock);

    d_->pin_count_.AssertZero();
    die_unequal(d_->total_ram_bytes_, 0u);
    die_unequal(d_->total_bytes_, 0u);
    die_unequal(d_->unpinned_blocks_.size(), 0u);

    LOGC(debug_pin)
        << "~BlockPool()"
        << " max_pin=" << d_->pin_count_.max_pins
        << " max_pinned_bytes=" << d_->pin_count_.max_pinned_bytes;

    logger_ << "class" << "BlockPool"
            << "event" << "destroy"
            << "max_pins" << d_->pin_count_.max_pins
            << "max_pinned_bytes" << d_->pin_count_.max_pinned_bytes;

    std::unique_lock<std::recursive_mutex> s_new_lock(s_new_mutex);
    s_blockpools.erase(
        std::find(s_blockpools.begin(), s_blockpools.end(), this));
}

PinnedByteBlockPtr
BlockPool::AllocateByteBlock(size_t size, size_t local_worker_id) {
    assert(local_worker_id < workers_per_host_);
    auto page_id = cache_alloc_page();
    // create tlx::CountingPtr, no need for special make_shared()-equivalent
    PinnedByteBlockPtr block_ptr(
        mem::GPool().make<ByteBlock>(this, page_id, size), local_worker_id);
    return block_ptr;
}

ByteBlockPtr BlockPool::MapExternalBlock(
    const foxxll::file_ptr& file, uint64_t offset, size_t size) {
  PANIC_NOT_IMPLEMENTED;
}

//! Pins a block by swapping it in if required.
PinRequestPtr BlockPool::PinBlock(const Block& block, size_t local_worker_id) {
    assert(local_worker_id < workers_per_host_);
    ByteBlock* block_ptr = block.byte_block().get();
    return PinRequestPtr(mem::GPool().make<PinRequest>(
        this, PinnedBlock(block, local_worker_id)));
}

std::pair<size_t, size_t> BlockPool::MaxMergeDegreePrefetch(size_t num_files) {
    size_t avail_bytes = hard_ram_limit() / workers_per_host_ / 2;
    size_t avail_blocks = avail_bytes / default_block_size;

    if (num_files >= avail_blocks) {
        // more files than blocks available -> partial merge of avail_bytes
        // Files with prefetch = 0, which is at most one block per File.
        return std::make_pair(avail_blocks, 0u);
    }
    else {
        // less files than available Blocks -> split prefetch size equally
        // among Files.
        return std::make_pair(num_files, avail_bytes / num_files);
    }
}

void PinRequest::OnComplete(foxxll::request* req, bool success) {
    return block_pool_->OnReadComplete(this, req, success);
}

void BlockPool::OnReadComplete(
    PinRequest* read, foxxll::request* req, bool success) {
}

void BlockPool::IncBlockPinCount(ByteBlock* block_ptr, size_t local_worker_id) {
    block_ptr->SetData((Byte*) cache_pin(block_ptr->page_id_));
}

void BlockPool::IntIncBlockPinCount(ByteBlock* block_ptr, size_t local_worker_id) {
    block_ptr->SetData((Byte*) cache_pin(block_ptr->page_id_));
}

void BlockPool::DecBlockPinCount(ByteBlock* block_ptr, size_t local_worker_id) {
    cache_unpin(block_ptr->page_id_, true);
}

void BlockPool::Data::IntUnpinBlock(
    BlockPool& bp, ByteBlock* block_ptr, size_t local_worker_id) {
}

size_t BlockPool::total_blocks() noexcept {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->int_total_blocks();
}

size_t BlockPool::Data::int_total_blocks() noexcept {

    LOG << "BlockPool::total_blocks()"
        << " pinned_blocks_=" << pin_count_.total_pins_
        << " unpinned_blocks_=" << unpinned_blocks_.size()
        << " writing_.size()=" << writing_.size()
        << " swapped_.size()=" << swapped_.size()
        << " reading_.size()=" << reading_.size();

    return pin_count_.total_pins_
           + unpinned_blocks_.size() + writing_.size()
           + swapped_.size() + reading_.size();
}

size_t BlockPool::hard_ram_limit() noexcept {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->hard_ram_limit_;
}

size_t BlockPool::total_bytes() noexcept {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->int_total_bytes();
}

size_t BlockPool::max_total_bytes() noexcept {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->max_total_bytes_;
}

size_t BlockPool::Data::int_total_bytes() noexcept {
    LOG << "BlockPool::total_bytes()"
        << " pinned_bytes_=" << pin_count_.total_pinned_bytes_
        << " unpinned_bytes_=" << unpinned_bytes_
        << " writing_bytes_=" << writing_bytes_
        << " swapped_bytes_=" << swapped_bytes_
        << " reading_bytes_=" << reading_bytes_;

    return pin_count_.total_pinned_bytes_
           + unpinned_bytes_ + writing_bytes_
           + swapped_bytes_ + reading_bytes_;
}

size_t BlockPool::pinned_blocks() noexcept {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->pin_count_.total_pins_;
}

size_t BlockPool::unpinned_blocks() noexcept {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->unpinned_blocks_.size();
}

size_t BlockPool::writing_blocks() noexcept {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->writing_.size();
}

size_t BlockPool::swapped_blocks() noexcept {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->swapped_.size();
}

size_t BlockPool::reading_blocks() noexcept {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->reading_.size();
}

void BlockPool::DestroyBlock(ByteBlock* block_ptr) {
}

void BlockPool::RequestInternalMemory(size_t size) {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->IntRequestInternalMemory(lock, size);
}

void BlockPool::Data::IntRequestInternalMemory(
    std::unique_lock<std::mutex>& lock, size_t size) {

    requested_bytes_ += size;

    LOGC(debug_mem)
        << "BlockPool::RequestInternalMemory()"
        << " size=" << size
        << " total_ram_bytes_=" << total_ram_bytes_
        << " writing_bytes_=" << writing_bytes_
        << " requested_bytes_=" << requested_bytes_
        << " soft_ram_limit_=" << soft_ram_limit_
        << " hard_ram_limit_=" << hard_ram_limit_
        << pin_count_
        << " unpinned_blocks_.size()=" << unpinned_blocks_.size()
        << " swapped_.size()=" << swapped_.size();

    while (soft_ram_limit_ != 0 &&
           unpinned_blocks_.size() &&
           total_ram_bytes_ + requested_bytes_ > soft_ram_limit_ + writing_bytes_)
    {
        // evict blocks: schedule async writing which increases writing_bytes_.
        IntEvictBlockLRU();
    }

    // wait up to 60 seconds for other threads to free up memory or pins
    static constexpr size_t max_retry = 60;
    size_t retry = max_retry;
    size_t last_writing_bytes = 0;

    // wait for memory change due to blocks begin written and deallocated.
    while (hard_ram_limit_ != 0 && total_ram_bytes_ + size > hard_ram_limit_)
    {
        while (hard_ram_limit_ != 0 &&
               unpinned_blocks_.size() &&
               total_ram_bytes_ + requested_bytes_ > hard_ram_limit_ + writing_bytes_)
        {
            // evict blocks: schedule async writing which increases writing_bytes_.
            IntEvictBlockLRU();
        }

        cv_memory_change_.wait_for(lock, std::chrono::seconds(1));

        LOGC(debug_mem)
            << "BlockPool::RequestInternalMemory() waiting for memory"
            << " total_ram_bytes_=" << total_ram_bytes_
            << " writing_bytes_=" << writing_bytes_
            << " requested_bytes_=" << requested_bytes_
            << " soft_ram_limit_=" << soft_ram_limit_
            << " hard_ram_limit_=" << hard_ram_limit_
            << pin_count_
            << " unpinned_blocks_.size()=" << unpinned_blocks_.size()
            << " swapped_.size()=" << swapped_.size();

        if (writing_bytes_ == 0 &&
            total_ram_bytes_ + requested_bytes_ > hard_ram_limit_) {

            LOG1 << "abort() due to out-of-pinned-memory ???"
                 << " total_ram_bytes_=" << total_ram_bytes_
                 << " writing_bytes_=" << writing_bytes_
                 << " requested_bytes_=" << requested_bytes_
                 << " soft_ram_limit_=" << soft_ram_limit_
                 << " hard_ram_limit_=" << hard_ram_limit_
                 << pin_count_
                 << " unpinned_blocks_.size()=" << unpinned_blocks_.size()
                 << " swapped_.size()=" << swapped_.size();

            if (writing_bytes_ == last_writing_bytes) {
                if (--retry == 0)
                    abort();
            }
            else {
                last_writing_bytes = writing_bytes_;
                retry = max_retry;
            }
        }
    }

    requested_bytes_ -= size;
    total_ram_bytes_ += size;
}

void BlockPool::AdviseFree(size_t size) {
    std::unique_lock<std::mutex> lock(mutex_);

    LOGC(debug_mem)
        << "BlockPool::AdviseFree() advice to free memory"
        << " size=" << size
        << " total_ram_bytes_=" << d_->total_ram_bytes_
        << " writing_bytes_=" << d_->writing_bytes_
        << " requested_bytes_=" << d_->requested_bytes_
        << " soft_ram_limit_=" << d_->soft_ram_limit_
        << " hard_ram_limit_=" << d_->hard_ram_limit_
        << d_->pin_count_
        << " unpinned_blocks_.size()=" << d_->unpinned_blocks_.size()
        << " swapped_.size()=" << d_->swapped_.size();

    while (d_->soft_ram_limit_ != 0 && d_->unpinned_blocks_.size() &&
           d_->total_ram_bytes_ + d_->requested_bytes_ + size > d_->hard_ram_limit_ + d_->writing_bytes_)
    {
        // evict blocks: schedule async writing which increases writing_bytes_.
        d_->IntEvictBlockLRU();
    }
}
void BlockPool::ReleaseInternalMemory(size_t size) {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->IntReleaseInternalMemory(size);
}

void BlockPool::Data::IntReleaseInternalMemory(size_t size) {

    LOGC(debug_mem)
        << "BlockPool::IntReleaseInternalMemory()"
        << " size=" << size
        << " total_ram_bytes_=" << total_ram_bytes_;

    die_unless(total_ram_bytes_ >= size);
    total_ram_bytes_ -= size;

    cv_memory_change_.notify_all();
}

void BlockPool::EvictBlock(ByteBlock* block_ptr) {
    std::unique_lock<std::mutex> lock(mutex_);

    die_unless(block_ptr->in_memory());

    die_unless(d_->unpinned_blocks_.exists(block_ptr));
    d_->unpinned_blocks_.erase(block_ptr);
    d_->unpinned_bytes_ -= block_ptr->size();

    d_->IntEvictBlock(block_ptr);
}

foxxll::request_ptr BlockPool::GetAnyWriting() {
    std::unique_lock<std::mutex> lock(mutex_);
    if (!d_->writing_.size()) return foxxll::request_ptr();
    return d_->writing_.begin()->second;
}

foxxll::request_ptr BlockPool::EvictBlockLRU() {
    std::unique_lock<std::mutex> lock(mutex_);
    return d_->IntEvictBlockLRU();
}

foxxll::request_ptr BlockPool::Data::IntEvictBlockLRU() {

    if (!unpinned_blocks_.size()) return foxxll::request_ptr();

    ByteBlock* block_ptr = unpinned_blocks_.pop();
    die_unless(block_ptr);
    unpinned_bytes_ -= block_ptr->size();

    return IntEvictBlock(block_ptr);
}

foxxll::request_ptr BlockPool::Data::IntEvictBlock(ByteBlock* block_ptr) {
  return foxxll::request_ptr();
}

void BlockPool::OnWriteComplete(
    ByteBlock* block_ptr, foxxll::request* req, bool success) {
}

void BlockPool::RunTask(const std::chrono::steady_clock::time_point& tp) {
    std::unique_lock<std::mutex> lock(mutex_);

    foxxll::stats_data stnow(*foxxll::stats::get_instance());
    foxxll::stats_data stf = stnow - d_->io_stats_first_;
    foxxll::stats_data stp = stnow - d_->io_stats_prev_;
    d_->io_stats_prev_ = stnow;

    double elapsed = static_cast<double>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            tp - d_->tp_last_).count()) / 1e6;
    d_->tp_last_ = tp;

    // LOG0 << stp;
    // LOG0 << stf;

    size_t unpinned_bytes = d_->unpinned_bytes_.hmax_update();
    size_t writing_bytes = d_->writing_bytes_.hmax_update();
    size_t reading_bytes = d_->reading_bytes_.hmax_update();
    size_t pinned_bytes = d_->pin_count_.total_pinned_bytes_.hmax_update();

    logger_ << "class" << "BlockPool"
            << "event" << "profile"
            << "total_blocks" << d_->int_total_blocks()
            << "total_bytes" << d_->total_bytes_.hmax_update()
            << "max_total_bytes" << d_->max_total_bytes_
            << "total_ram_bytes" << d_->total_ram_bytes_.hmax_update()
            << "ram_bytes"
            << (unpinned_bytes + pinned_bytes + writing_bytes + reading_bytes)
            << "pinned_blocks" << d_->pin_count_.total_pins_
            << "pinned_bytes" << pinned_bytes
            << "unpinned_blocks" << d_->unpinned_blocks_.size()
            << "unpinned_bytes" << unpinned_bytes
            << "swapped_blocks" << d_->swapped_.size()
            << "swapped_bytes" << d_->swapped_bytes_.hmax_update()
            << "max_pinned_blocks" << d_->pin_count_.max_pins
            << "max_pinned_bytes" << d_->pin_count_.max_pinned_bytes
            << "writing_blocks" << d_->writing_.size()
            << "writing_bytes" << writing_bytes
            << "reading_blocks" << d_->reading_.size()
            << "reading_bytes" << reading_bytes
            << "rd_ops_total" << stf.get_read_count()
            << "rd_bytes_total" << stf.get_read_bytes()
            << "wr_ops_total" << stf.get_write_count()
            << "wr_bytes_total" << stf.get_write_bytes()
            << "rd_ops" << stp.get_read_count()
            << "rd_bytes" << stp.get_read_bytes()
            << "rd_speed" << static_cast<double>(stp.get_read_bytes()) / elapsed
            << "wr_ops" << stp.get_write_count()
            << "wr_bytes" << stp.get_write_bytes()
            << "wr_speed" << static_cast<double>(stp.get_write_bytes()) / elapsed
            << "disk_allocation" << d_->bm_->current_allocation();
}

size_t BlockPool::next_file_id() {
    return ++d_->next_file_id_;
}

} // namespace data
} // namespace thrill

/******************************************************************************/
