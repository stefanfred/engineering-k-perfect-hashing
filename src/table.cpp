#include <tlx/cmdline_parser.hpp>

#include "ChdContender.hpp"
#include "HashDisplaceContender.hpp"
#include "PaCHashContender.hpp"
#include "ThresholdBasedBumpingContender.hpp"
#include "ThresholdBasedBumpingConsensusContender.hpp"
#include "ThresholdBasedBumpingOldContender.hpp"
#include "KRecSplitContender.hpp"

#include <OptimalBucketFunction.hpp>
#include <CompactEncoding.hpp>
#include <RiceEncoding.hpp>

int main(int argc, char** argv) {
    size_t N = 5e6;
    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "numKeys", N, "Number of objects");
    cmd.add_bytes('q', "numQueries", Contender::numQueries, "Number of queries to perform");
    cmd.add_bytes('t', "numThreads", Contender::numThreads, "Number of threads to run benchmarks with");
    cmd.add_flag('T', "skipTests", Contender::skipTests, "Skip testing PHF for validity");
    cmd.add_bytes('s', "seed", Contender::seed, "Seed for random inputs");
    if (!cmd.process(argc, argv)) {
        return 1;
    }
    { ChdContender(N, 10, 0.97, 8, false).run(); }
    { HashDisplaceContender<10, kphf::HashDisplace::OptimalBucketFunction<10>, kphf::HashDisplace::RiceEncoding>(N, 12).run(); }
    { HashDisplaceContender<10, kphf::HashDisplace::OptimalBucketFunction<10>, kphf::HashDisplace::CompactEncoding>(N, 12).run(); }
    { PaCHashContender(N, 10, 10).run(); }
    { ThresholdBasedBumpingOldContender<10>(N).run(); }
    { ThresholdBasedBumpingContender<10, 5, false>(N, 2.0).run(); }
    { ThresholdBasedBumpingContender<10, 4, true>(N, 2.0).run(); }
    { ThresholdBasedBumpingConsensusContender<10, 4>(N, 2.0).run(); }
    { KRecSplitContender<10, 2>(N, 2000).run(); }

    { HashDisplaceContender<1000, kphf::HashDisplace::OptimalBucketFunction<1000>, kphf::HashDisplace::RiceEncoding>(N, 250).run(); }
    { HashDisplaceContender<1000, kphf::HashDisplace::OptimalBucketFunction<1000>, kphf::HashDisplace::CompactEncoding>(N, 250).run(); }
    { PaCHashContender(N, 1000, 1000).run(); }
    { ThresholdBasedBumpingOldContender<1000>(N).run(); }
    { ThresholdBasedBumpingContender<1000, 9, false>(N, 1.2).run(); }
    { ThresholdBasedBumpingContender<1000, 7, true>(N, 1.2).run(); }
    { ThresholdBasedBumpingConsensusContender<1000, 7>(N, 1.2).run(); }
    { KRecSplitContender<1000, 2>(N, 6000).run(); }
    return 0;
}
