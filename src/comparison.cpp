#include <tlx/cmdline_parser.hpp>

#include "ChdContender.hpp"
#include "HashDisplaceContender.hpp"
#include "PaCHashContender.hpp"

int main(int argc, char** argv) {
    double loadFactor = 0.8;
    size_t N = 5e6;
    size_t k = 8;
    bool chd = false;
    bool hashdisplace = false;
    bool pachash = false;

    tlx::CmdlineParser cmd;
    cmd.add_double('l', "loadFactor", loadFactor, "Load Factor");
    cmd.add_bytes('n', "numKeys", N, "Number of objects");
    cmd.add_bytes('q', "numQueries", Contender::numQueries, "Number of queries to perform");
    cmd.add_bytes('t', "numThreads", Contender::numThreads, "Number of threads to run benchmarks with");
    cmd.add_flag('T', "skipTests", Contender::skipTests, "Skip testing PHF for validity");
    cmd.add_bytes('k', "k", k, "Number of collisions per output value");

    cmd.add_flag("chd", chd, "Execute CHD benchmark");
    cmd.add_flag("hashDisplace", hashdisplace, "Execute hash displace benchmark");
    cmd.add_flag("pachash", pachash, "Execute PaCHash benchmark");

    if (!cmd.process(argc, argv)) {
        return 1;
    }
    if (chd) {
        chdContenderRunner(N, k, 0.95);
    }
    if (hashdisplace) {
        hashDisplaceContenderRunner(N, k);
    }
    if (pachash) {
        paCHashContenderRunner(N, k);
    }
    return 0;
}
