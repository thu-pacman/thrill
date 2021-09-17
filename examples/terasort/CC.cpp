/*******************************************************************************
 * examples/terasort/terasort.cpp
 *
 * Part of Project Thrill - http://project-thrill.org
 *
 * Copyright (C) 2016 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <thrill/api/generate.hpp>
#include <thrill/api/read_binary.hpp>
#include <thrill/api/group_by_key.hpp>
#include <thrill/api/inner_join.hpp>
#include <thrill/api/reduce_by_key.hpp>
#include <thrill/api/cache.hpp>
#include <thrill/api/union.hpp>
#include <thrill/api/size.hpp>
#include <thrill/api/sort.hpp>
#include <thrill/api/write_binary.hpp>
#include <thrill/common/logger.hpp>
#include <thrill/common/string.hpp>
#include <thrill/api/read_lines.hpp>
#include <tlx/cmdline_parser.hpp>

#include <tlx/string/hexdump.hpp>
#include <tlx/string/parse_si_iec_units.hpp>

#include <algorithm>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <cstdlib>
#include <chrono>

using namespace thrill;               // NOLINT

struct EdgeLink {
    size_t src, tgt;
};

using OutgoingLinks = std::vector<size_t>;
using RankedPage = std::pair<size_t, size_t>;
using OutgoingLinksRank = std::pair<std::vector<size_t>, size_t>;
using LinkedPage = std::pair<size_t, std::vector<size_t>>;

constexpr bool UseLocationDetection = true;

static void RunConnectedComponent(
    api::Context &ctx,
    const std::vector<std::string> &input_filelist) {
    
    ctx.enable_consume();

    auto lines = ReadLines(ctx, input_filelist).Keep();
    auto line_size = lines.Size();

    if (ctx.my_rank() == 0) {
        std::cout << "whole count is " << line_size << std::endl;
    }

    common::StatsTimerStart timer;

    auto edges = lines.FlatMap<EdgeLink>(
        [](const std::string &line, auto emit) -> void {
            auto pos = line.find(' ');
            size_t src = std::atoll(line.c_str());
            size_t dst = std::atoll(line.c_str() + pos + 1);
            emit(EdgeLink{src, dst});
            emit(EdgeLink{dst, src});
        }
    );

    auto links = edges.GroupByKey<LinkedPage>(
        [](const EdgeLink &p) { return p.src; },
        [all = std::vector<size_t> ()](auto &r, const size_t &pid) mutable {
            all.clear();
            size_t id;
            while (r.HasNext()) {
                auto item = r.Next();
                id = item.src;
                all.push_back(item.tgt);
            }
            return std::make_pair(id, all);
        }).Cache().KeepForever();
    
    DIA<RankedPage> ranks = links.Map([](const LinkedPage &lp) {
        return std::make_pair(lp.first, *std::min_element(lp.second.begin(), lp.second.end()));
    }).Keep();

    auto time_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; true; ++i) {
        auto out_rank = InnerJoin(
            LocationDetectionFlag<UseLocationDetection>(),
            links, ranks,
            [](const LinkedPage &lp) { return lp.first; },
            [](const RankedPage &rp) { return rp.first; },
            [](const LinkedPage &lp, const RankedPage &rp) {
                return std::make_tuple(lp.first, lp.second, rp.second);
            }
        );
        auto tdelta = out_rank.FlatMap<std::pair<size_t, size_t>>([](const std::tuple<size_t, std::vector<size_t>, size_t> &item, auto emit) {
            emit(std::make_pair(std::get<0>(item), std::get<2>(item)));
            for (auto v: std::get<1>(item)) {
                if (v > std::get<2>(item)) {
                    emit(std::make_pair(v, std::get<2>(item)));
                }
            }
        });

        auto delta = tdelta.ReduceByKey(
            [](const std::pair<size_t, size_t> &in) -> size_t {
                return in.first;
            },
            [](const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b) {
                size_t res = a.second;
                if (res > b.second) res = b.second;
                return std::make_pair(a.first, res);
            });
        auto new_delta_tuple = InnerJoin(
            LocationDetectionFlag<UseLocationDetection>(),
            delta, ranks,
            [](const RankedPage &lp) { return lp.first; },
            [](const RankedPage &rp) { return rp.first; },
            [](const RankedPage &lp, const RankedPage &rp) {
                return std::make_tuple(lp.first, lp.second, rp.second);
            }
        ).Filter([](const std::tuple<size_t, size_t, size_t> &in) { return std::get<1>(in) < std::get<2>(in); });

        auto new_delta =  new_delta_tuple.Map([](const std::tuple<size_t, size_t, size_t> &in) {
            return std::make_pair(std::get<0>(in), std::get<1>(in));
        }).Keep();
        auto updates = new_delta.Size();

        auto new_ranks = ranks.Union(new_delta).ReduceByKey(
            [](const RankedPage &in) -> size_t {
                return in.first;
            },
            [](const RankedPage &a, const RankedPage &b) {
                if (a.second <= b.second) {
                    return std::make_pair(a.first, a.second);
                } else {
                    return std::make_pair(a.first, b.second);
                }
            });
        ranks = new_ranks.Keep();
        auto time_now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = time_now - time_start;
        if (ctx.my_rank() == 0) {
            std::cout << "step " << i << ", time: " << diff.count() << ", update: " << updates << std::endl;
        }
        if (updates == 0) break;
    }

    //auto num_components = ranks.Size();
    ranks.Execute();

    ctx.net.Barrier();
    timer.Stop();
    if (ctx.my_rank() == 0) {
        //std::cout << "num components: " << num_components << std::endl;
        std::cout << "time= " << timer << std::endl;
    }
}


int main(int argc, char* argv[]) {

    tlx::CmdlineParser clp;

    std::vector<std::string> input;
    clp.add_param_stringlist("input", input,
                             "input file pattern(s)");

    if (!clp.process(argc, argv)) {
        return -1;
    }

    clp.print_result();

    return api::Run(
        [&](api::Context& ctx) {
            RunConnectedComponent(ctx, input);
        });
}

/******************************************************************************/
