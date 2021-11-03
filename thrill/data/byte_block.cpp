/*******************************************************************************
 * thrill/data/byte_block.cpp
 *
 * Part of Project Thrill - http://project-thrill.org
 *
 * Copyright (C) 2015-2018 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <thrill/data/block_pool.hpp>
#include <thrill/data/byte_block.hpp>
#include <thrill/data/file.hpp>
#include <thrill/mem/pool.hpp>

#include <tlx/string/join_generic.hpp>

#include <sstream>
#include <string>

namespace thrill {
namespace data {

size_t start_block_size = 4 * 1024;
// size_t default_block_size = 2 * 1024 * 1024;
size_t default_block_size = 4 * 1024;

size_t File::default_prefetch_size_ = 2 * default_block_size;

ByteBlock::ByteBlock(BlockPool* block_pool, Byte* data, size_t size)
    : data_(data), size_(size),
      block_pool_(block_pool),
      pin_count_(block_pool_->workers_per_host())
{ }

ByteBlock::ByteBlock(
    BlockPool* block_pool, size_t page_id, size_t size)
    : data_(nullptr), size_(size),
      block_pool_(block_pool),
      pin_count_(block_pool_->workers_per_host()),
      page_id_(page_id)
{ }

void ByteBlock::Deleter::operator () (ByteBlock* bb) const {
    sLOG << "ByteBlock[" << bb << "]::deleter()"
         << "pin_count_" << bb->pin_count_str();
    assert(bb->total_pins() == 0);
    assert(bb->reference_count() == 0);

    // call BlockPool's DestroyBlock() to de-register ByteBlock and free data
    assert(bb->block_pool_);
    bb->block_pool_->DestroyBlock(bb);

    sLOG << "ByteBlock[ " << bb << "]::destroy()";
    mem::GPool().destroy(bb);
}

void ByteBlock::Deleter::operator () (const ByteBlock* bb) const {
    return operator () (const_cast<ByteBlock*>(bb));
}

std::string ByteBlock::pin_count_str() const {
    return "[" + tlx::join(',', pin_count_) + "]";
}

void ByteBlock::IncPinCount(size_t local_worker_id) {
    return block_pool_->IncBlockPinCount(this, local_worker_id);
}

void ByteBlock::DecPinCount(size_t local_worker_id) {
    return block_pool_->DecBlockPinCount(this, local_worker_id);
}

void ByteBlock::OnWriteComplete(foxxll::request* req, bool success) {
    return block_pool_->OnWriteComplete(this, req, success);
}

std::ostream& operator << (std::ostream& os, const ByteBlock& b) {
    os << "[ByteBlock" << " " << &b
       << " data_=" << static_cast<const void*>(b.data_)
       << " size_=" << b.size_
       << " block_pool_=" << b.block_pool_
       << " total_pins_=" << b.total_pins_
       << " page_id_=" << b.page_id_;
    return os << "]";
}

} // namespace data
} // namespace thrill

/******************************************************************************/
