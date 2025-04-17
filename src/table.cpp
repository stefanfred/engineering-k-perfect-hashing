#include <tlx/cmdline_parser.hpp>

#include "DispatchK.h"
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

template <size_t k>
struct TableRunner {
    void operator() (size_t N) const {
        // { ChdContender(N, k, 0.97, keysPerBucket, false).run(); }
        // { HashDisplaceContender<k, kphf::HashDisplace::OptimalBucketFunction<k>, kphf::HashDisplace::RiceEncoding>(N, bucketSize).run(); }
        // { HashDisplaceContender<k, kphf::HashDisplace::OptimalBucketFunction<k>, kphf::HashDisplace::CompactEncoding>(N, bucketSize).run(); }
        // { PaCHashContender(N, k, a).run(); }
        // { ThresholdBasedBumpingOldContender<k>(N).run(); }
        // { ThresholdBasedBumpingContender<k, thresholdSize, false>(N, overload).run(); }
        // { ThresholdBasedBumpingContender<k, thresholdSize, true>(N, overload).run(); }
        // { ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, overload).run(); }
        // { KRecSplitContender<k, leafSize>(N, bucketSize).run(); }
    }
};

int main(int argc, char** argv) {
    size_t N = 5e6;
    size_t k = 8;

    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "numKeys", N, "Number of objects");
    cmd.add_bytes('q', "numQueries", Contender::numQueries, "Number of queries to perform");
    cmd.add_bytes('t', "numThreads", Contender::numThreads, "Number of threads to run benchmarks with");
    cmd.add_flag('T', "skipTests", Contender::skipTests, "Skip testing PHF for validity");
    cmd.add_bytes('s', "seed", Contender::seed, "Seed for random inputs");
    cmd.add_bytes('k', "k", k, "Number of collisions per output value");
    if (!cmd.process(argc, argv)) {
        return 1;
    }
    dispatchK<TableRunner>(k, N);
    return 0;
}
