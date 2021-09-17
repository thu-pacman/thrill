/*******************************************************************************
 * examples/page_rank/page_rank_run.cpp
 *
 * Part of Project Thrill - http://project-thrill.org
 *
 * Copyright (C) 2016 Timo Bingmann <tb@panthema.net>
 * Copyright (C) 2016 Alexander Noe <aleexnoe@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <thrill/api/cache.hpp>
#include <thrill/api/group_by_key.hpp>
#include <thrill/api/group_to_index.hpp>
#include <thrill/api/max.hpp>
#include <thrill/api/read_lines.hpp>
#include <thrill/api/sum.hpp>
#include <thrill/api/write_lines.hpp>
#include <thrill/api/zip_with_index.hpp>
#include <thrill/common/logger.hpp>
#include <thrill/common/stats_timer.hpp>
#include <tlx/cmdline_parser.hpp>
#include <thrill/api/collapse.hpp>
#include <thrill/api/inner_join.hpp>
#include <thrill/api/reduce_by_key.hpp>
#include <thrill/api/reduce_to_index.hpp>
#include <thrill/api/size.hpp>
#include <thrill/api/zip.hpp>

#include <tlx/string/join_generic.hpp>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

using namespace thrill;              // NOLINT

using PageId = std::size_t;
using Rank = double;

struct PagePageLink {
    PageId src, tgt;

    friend std::ostream& operator << (std::ostream& os, const PagePageLink& a) {
        return os << '(' << a.src << '>' << a.tgt << ')';
    }
} TLX_ATTRIBUTE_PACKED;

//! A pair (page, rank)
struct PageRankPair {
    PageId page;
    Rank   rank;

    friend std::ostream& operator << (std::ostream& os, const PageRankPair& a) {
        return os << '(' << a.page << '|' << a.rank << ')';
    }
} TLX_ATTRIBUTE_PACKED;

using PageRankStdPair = std::pair<PageId, Rank>;
using OutgoingLinks = std::vector<PageId>;
using OutgoingLinksRank = std::pair<std::vector<PageId>, Rank>;
using LinkedPage = std::pair<PageId, OutgoingLinks>;
using RankedPage = std::pair<PageId, Rank>;

static constexpr double dampening = 0.85;

template <const bool UseLocationDetection = true, typename InStack>
auto PageRankJoinSelf(const DIA<LinkedPage, InStack>& links, size_t iterations) {

    api::Context& ctx = links.context();
    // initialize all ranks to 1.0 / n: (url, rank)

    DIA<RankedPage> ranks = links.Map([](const LinkedPage &lp) {
        return std::make_pair(lp.first, Rank(1.0));
    });

    // do iterations
    for (size_t iter = 0; iter < iterations; ++iter) {

        // for all outgoing link, get their rank contribution from all
        // links by doing:
        //
        // 1) group all outgoing links with rank of its parent page: (Zip)
        // ([linked_url, linked_url, ...], rank_parent)
        //
        // 2) compute rank contribution for each linked_url: (FlatMap)
        // (linked_url, rank / outgoing.size)

        auto outs_rank = InnerJoin(
            LocationDetectionFlag<UseLocationDetection>(),
            links, ranks,
            [](const LinkedPage& lp) { return lp.first; },
            [](const RankedPage& rp) { return rp.first; },
            [](const LinkedPage& lp, const RankedPage& rp) {
                return std::make_pair(lp.second, rp.second);
            });

        //if (debug && iter == 1) {
        //    outs_rank
        //    .Map([](const OutgoingLinksRank& ol) {
        //             return tlx::join(',', ol.first)
        //             + " <- " + std::to_string(ol.second);
        //         })
        //    .Print("outs_rank");
        //}

        auto contribs = outs_rank.template FlatMap<PageRankStdPair>(
            [](const OutgoingLinksRank& p, auto emit) {
                if (p.first.size() > 0) {
                    Rank rank_contrib = p.second / static_cast<double>(p.first.size());
                    for (const PageId& tgt : p.first)
                        emit(std::make_pair(tgt, rank_contrib));
                }
            });

        // reduce all rank contributions by adding all rank contributions and
        // compute the new rank: (url, rank)
        ranks =
            contribs
            .ReducePair(
                [](const Rank& p1, const Rank& p2) {
                    return p1 + p2;
                })
            .Map([](const PageRankStdPair& p) {
                     return std::make_pair(
                         p.first,
                         dampening * p.second + (1 - dampening));
                 }).Collapse();

        ranks.Execute();
    }
    return ranks;
}

struct PageRankLineParser {
    PagePageLink operator () (const std::string& input) const {
        // parse "source\ttarget\n" lines
        char* endptr;
        unsigned long src = std::strtoul(input.c_str(), &endptr, 10);
        die_unless(endptr && *endptr == ' ' &&
                   "Could not parse src tgt line");
        unsigned long tgt = std::strtoul(endptr + 1, &endptr, 10);
        die_unless(endptr && *endptr == 0 &&
                   "Could not parse src tgt line");
        return PagePageLink { src, tgt };
    }
};

static void RunJoinPageRankEdgePerLine(
    api::Context& ctx,
    const std::vector<std::string>& input_path, const std::string& output_path,
    size_t iterations) {
    ctx.enable_consume();

    common::StatsTimerStart timer;

    const bool UseLocationDetection = true;

    // read input file and create links in this format:
    //
    // url linked_url
    // url linked_url
    // url linked_url
    // ...
    auto input =
        ReadLines(ctx, input_path)
        .Map(PageRankLineParser());

    // aggregate all outgoing links of a page in this format: by index
    // ([linked_url, linked_url, ...])

    // group outgoing links from input file

    auto links = input.GroupByKey<LinkedPage>(
        [](const PagePageLink& p) { return p.src; },
        [all = std::vector<PageId> ()](auto& r, const PageId& pid) mutable {
            all.clear();
            PageId id;
            while (r.HasNext()) {
                auto item = r.Next();
                id = item.src;
                all.push_back(item.tgt);
            }
            return std::make_pair(id, all);
        }).Cache().KeepForever();

    // perform actual page rank calculation iterations

    auto ranks = PageRankJoinSelf<UseLocationDetection>(
        links, iterations);

    // construct output as "pageid: rank"

    if (output_path.size()) {
        ranks.Map([](const RankedPage& rp) {
                      return tlx::ssprintf("%zu %g", rp.first, rp.second);
                  }).WriteLines(output_path);
    }
    else {
        ranks.Execute();
    }

    //ctx.net.Barrier();
    //timer.Stop();

    //if (ctx.my_rank() == 0) {
    //    if (UseLocationDetection) {
    //        std::cout << "RESULT benchmark=pagerank_gen detection=ON"
    //             << " iterations=" << iterations
    //             << " time=" << timer
    //             << " traffic= " << ctx.net_manager().Traffic()
    //             << " hosts=" << ctx.num_hosts();
    //    }
    //    else {
    //        std::cout << "RESULT benchmark=pagerank_gen detection=OFF"
    //             << " iterations=" << iterations
    //             << " time=" << timer
    //             << " traffic=" << ctx.net_manager().Traffic()
    //             << " hosts=" << ctx.num_hosts();
    //    }
    //}
}



int main(int argc, char* argv[]) {

    tlx::CmdlineParser clp;

    std::string output_path;
    clp.add_string('o', "output", output_path,
                   "output file pattern");

    size_t iter = 10;
    clp.add_size_t('n', "iterations", iter, "PageRank iterations, default: 10");

    std::vector<std::string> input_path;
    clp.add_param_stringlist("input", input_path,
                             "input file pattern(s)");

    if (!clp.process(argc, argv)) {
        return -1;
    }

    clp.print_result();

    //die_unless(!generate || input_path.size() == 1);

    return api::Run(
        [&](api::Context& ctx) {
            return RunJoinPageRankEdgePerLine(
                ctx, input_path, output_path, iter);
        });
}

/******************************************************************************/
