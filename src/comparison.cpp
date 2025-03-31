#include <tlx/cmdline_parser.hpp>

#include "ChdContender.hpp"
#include "HashDisplaceContender.hpp"
#include "PaCHashContender.hpp"
#include "ThresholdBasedBumpingContender.hpp"
#include "ThresholdBasedBumpingConsensusContender.hpp"
#include "KRecSplitContender.hpp"

int main(int argc, char** argv) {
    size_t N = 5e6;
    size_t k = 8;
    bool chd = false;
    bool hashdisplace = false;
    bool pachash = false;
    bool thresholdbasedbumping = false;
    bool thresholdBasedBumpingConsensus = false;
    bool kRecSplit = false;

    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "numKeys", N, "Number of objects");
    cmd.add_bytes('q', "numQueries", Contender::numQueries, "Number of queries to perform");
    cmd.add_bytes('t', "numThreads", Contender::numThreads, "Number of threads to run benchmarks with");
    cmd.add_flag('T', "skipTests", Contender::skipTests, "Skip testing PHF for validity");
    cmd.add_bytes('k', "k", k, "Number of collisions per output value");

    cmd.add_flag("chd", chd, "Execute CHD benchmark");
    cmd.add_flag("hashDisplace", hashdisplace, "Execute hash displace benchmark");
    cmd.add_flag("pachash", pachash, "Execute PaCHash benchmark");
    cmd.add_flag("thresholdBased", thresholdbasedbumping, "Execute threshold based bumping benchmark");
    cmd.add_flag("thresholdBasedConsensus", thresholdBasedBumpingConsensus, "Execute threshold based bumping consensus benchmark");
    cmd.add_flag("kRecSplit", kRecSplit, "Execute KRecSplit benchmark");

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
    if (thresholdbasedbumping) {
        thresholdBasedBumpingContenderRunner(N, k);
    }
    if (thresholdBasedBumpingConsensus) {
        thresholdBasedBumpingConsensusContenderRunner(N, k);
    }
    if (kRecSplit) {
        kRecSplitContenderRunner(N, k);
    }
    return 0;
}
