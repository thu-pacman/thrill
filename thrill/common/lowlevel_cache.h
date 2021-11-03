//
// Created by ybw on 2021/11/3.
//

#ifndef THRILL_LOWLEVEL_CACHE_H
#define THRILL_LOWLEVEL_CACHE_H

#pragma once
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif
extern size_t *__client_cpus;
extern size_t __num_client_cpus;
extern const size_t CACHE_PAGE_SIZE;
extern size_t CACHE_NUM_PAGES;

// load env: CACHE_PHY_SIZE, CACHE_VIRT_SIZE, CACHE_CONFIG, CACHE_NUM_CLIENTS
extern __attribute__((constructor)) void init();
extern __attribute__((destructor)) void deinit();
extern void *cache_pin(size_t page_id);
extern void cache_unpin(size_t page_id, bool is_write);
extern size_t cache_alloc_page();
extern void cache_free_page(size_t page_id);
extern size_t cache_pin_count(size_t page_id);

#ifdef __cplusplus
}
#endif

#define PANIC_NOT_IMPLEMENTED doPanic(__FILE__, __LINE__, "Not implemented")

#include <iostream>
#include <memory>

static void doPanic(std::string file, int line, std::string givenMsg) {
  std::cout << "Panic: " << file << ":" << line << " " << givenMsg;
  exit(-1);
}

class MyBlock {
  size_t pageId_;
public:
  MyBlock() : pageId_(cache_alloc_page()) {}
  virtual ~MyBlock() { cache_free_page(pageId_); }
};

class MyPinnedBlock {
  size_t pageId_;
  void* ptr_;
  bool dirty_ = false;

public:
  MyPinnedBlock(size_t pageId) : pageId_(pageId), ptr_(cache_pin(pageId)) {}
  virtual ~MyPinnedBlock() { cache_unpin(pageId_, dirty_); }
};

#endif // THRILL_LOWLEVEL_CACHE_H
