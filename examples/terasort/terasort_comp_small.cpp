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
#include <thrill/api/size.hpp>
#include <thrill/api/sort.hpp>
#include <thrill/api/write_binary.hpp>
#include <thrill/common/logger.hpp>
#include <thrill/common/string.hpp>
#include <tlx/cmdline_parser.hpp>

#include <tlx/string/hexdump.hpp>
#include <tlx/string/parse_si_iec_units.hpp>

#include <algorithm>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace thrill;               // NOLINT

struct RecordSigned {
    char key[10];
    char value[90];

    // this sorted by _signed_ characters, which is the same as what some
    // Java/Scala TeraSorts do.
    bool operator < (const RecordSigned& b) const {
        return std::lexicographical_compare(key, key + 10, b.key, b.key + 10);
    }
    friend std::ostream& operator << (std::ostream& os, const RecordSigned& c) {
        return os << tlx::hexdump(c.key, 10);
    }
} TLX_ATTRIBUTE_PACKED;


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
            ctx.enable_consume();

            common::StatsTimerStart timer;

            auto rr = ReadBinary<RecordSigned>(ctx, input).Keep();

            auto r = rr.Sort();

            if (output.size())
                r.WriteBinary(output);
            else
                r.Size();

            ctx.net.Barrier();
            timer.Stop();
            //if (ctx.my_rank() == 0) {
            //    std::cout << "compute time is " << inner_timer << std::endl;
            //}

            //ctx.net.Barrier();
            //timer.Stop();
            //if (ctx.my_rank() == 0) {
            //    auto traffic = ctx.net_manager().Traffic();
            //    std::cout << "RESULT"
            //         << " benchmark=terasort"
            //         << " time=" << timer
            //         << " traffic=" << traffic.total()
            //         << " hosts=" << ctx.num_hosts();
            //}
        });
}

/******************************************************************************/
