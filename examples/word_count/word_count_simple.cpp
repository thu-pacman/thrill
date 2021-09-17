/*******************************************************************************
 * examples/word_count/word_count_simple.cpp
 *
 * Part of Project Thrill - http://project-thrill.org
 *
 * Copyright (C) 2016 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <thrill/api/read_lines.hpp>
#include <thrill/api/reduce_by_key.hpp>
#include <thrill/api/write_lines.hpp>
#include <thrill/api/generate.hpp>
#include <tlx/string/split_view.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <utility>
#include <chrono>
#include <vector>

void WordCount(thrill::Context& ctx,
               std::vector<std::string>& input, std::string output) {
    using Pair = std::pair<std::string, size_t>;
    auto lines = Generate(ctx, input.size(), [&input](size_t index) { return input[index]; });
    auto words =
        lines.template FlatMap<std::string>(
            // flatmap lambda: split and emit each word
            [](const std::string& line, auto emit) {
        size_t old_pos = 0;
        size_t pos;
        for (pos = 0; pos < line.size(); pos++) {
			if (line[pos] == ' ') {
				if (pos > old_pos) {
					emit(line.substr(old_pos, pos - old_pos));
				}
				old_pos = pos + 1;
			}
		}
        if (pos > old_pos) emit(line.substr(old_pos, pos - old_pos));
            });
    auto word_pairs = words.Filter([](const std::string &s) { return !s.empty(); })
            .template Map([](const std::string &s) {  return Pair{s, 1}; });
    auto ans = word_pairs.ReduceByKey(
        // key extractor: the word string
        [](const Pair& p) { return p.first; },
        // commutative reduction: add counters
        [](const Pair& a, const Pair& b) {
            return Pair(a.first, a.second + b.second);
        })
    .Map([](const Pair& p) {
             return p.first + ": "
             + std::to_string(p.second);
         })
    .Execute();
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <input> <output>" << std::endl;
        return -1;
    }
	std::vector<std::string> lines;
	std::ifstream file(argv[1]);
	std::string line;
	while (std::getline(file, line)) {
		if (line.size() == 1) {
			continue;
		}
		lines.emplace_back(line);
	}
	auto start = std::chrono::high_resolution_clock::now();
   auto ret = thrill::Run(
        [&](thrill::Context& ctx) { WordCount(ctx, lines, argv[2]); });
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "time: " << duration.count() << std::endl;
   return ret;
}

/******************************************************************************/
