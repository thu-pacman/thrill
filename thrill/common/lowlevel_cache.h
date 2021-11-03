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

#ifdef __cplusplus
}
#endif

#endif // THRILL_LOWLEVEL_CACHE_H
