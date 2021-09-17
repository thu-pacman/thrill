/*******************************************************************************
 * examples/word_count/word_count_run.cpp
 *
 * Part of Project Thrill - http://project-thrill.org
 *
 * Copyright (C) 2015 Alexander Noe <aleexnoe@gmail.com>
 * Copyright (C) 2016 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#include <thrill/api/cache.hpp>
#include <thrill/api/read_lines.hpp>
#include <thrill/api/write_lines.hpp>
#include <thrill/common/stats_timer.hpp>
#include <thrill/common/string.hpp>
#include <tlx/cmdline_parser.hpp>
#include <thrill/api/size.hpp>
#include <thrill/api/reduce_by_key.hpp>
#include <tlx/string/split_view.hpp>

#include <algorithm>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace thrill;               // NOLINT

/******************************************************************************/
// Run methods

using WordCountPair = std::pair<std::string, size_t>;

//! The most basic WordCount user program: reads a DIA containing std::string
//! words, and returns a DIA containing WordCountPairs.
template <typename InputStack>
auto WordCount(const DIA<std::string, InputStack>& input) {
    auto words = input.template FlatMap<tlx::string_view>(
        [](const std::string& line, auto emit) -> void {
            /* map lambda: emit each word */
            tlx::split_view(
                ' ', line, [&](const tlx::string_view& sv) {
                    emit(sv);
                });
        });

    auto word_pairs = words.Filter([](const tlx::string_view &sv) { return sv.size() != 0; })
            .template Map([](const tlx::string_view &sv) { return WordCountPair{sv.to_string(), 1}; });

    return word_pairs.ReduceByKey(
        [](const WordCountPair& in) -> std::string {
            /* reduction key: the word string */
            return in.first;
        },
        [](const WordCountPair& a, const WordCountPair& b) -> WordCountPair {
            /* associative reduction operator: add counters */
            return WordCountPair(a.first, a.second + b.second);
        });
}

static void RunWordCount(
    api::Context& ctx,
    const std::vector<std::string>& input_filelist, const std::string& output) {
    ctx.enable_consume();

    //common::StatsTimerStart timer;

    auto lines = ReadLines(ctx, input_filelist);

    auto word_pairs = WordCount(lines);

    if (output.size()) {
        word_pairs
        .Map([](const WordCountPair& wc) {
                 return wc.first + " " + std::to_string(wc.second);
             })
        .WriteLines(output);
    }
    else {
        word_pairs.Execute();
    }
    //ctx.net.Barrier();
    //timer.Stop();
    //if (ctx.my_rank() == 0) {
    //    std::cout << "RESULT"
    //         << " benchmark=wordcount"
    //         << " files=" << input_filelist.size()
    //         << " time=" << timer
    //         << " traffic=" << ctx.net_manager().Traffic()
    //         << " hosts=" << ctx.num_hosts();
    //}

}

/******************************************************************************/

int main(int argc, char* argv[]) {

    tlx::CmdlineParser clp;

    std::string output;
    clp.add_string('o', "output", output,
                   "output file pattern");

    std::vector<std::string> input;
    clp.add_param_stringlist("input", input,
                             "input file pattern(s)");

    if (!clp.process(argc, argv)) {
        return -1;
    }

    clp.print_result();

    return api::Run(
        [&](api::Context& ctx) {
            RunWordCount(ctx, input, output);
        });
}

/******************************************************************************/
