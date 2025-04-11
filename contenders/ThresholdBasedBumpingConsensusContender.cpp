#include "ThresholdBasedBumpingConsensusContender.hpp"

#include <tlx/math/integer_log2.hpp>
#include <DispatchK.h>

template <size_t k, size_t threshold_size>
void dispatchOverloadsConsensus(size_t N) {
    {ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, 1.05).run();}
    {ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, 1.10).run();}
    {ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, 1.20).run();}
    {ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, 1.30).run();}
    {ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, 1.40).run();}
    {ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, 1.50).run();}
    {ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, 1.60).run();}
    {ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, 1.80).run();}
    {ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, 2.00).run();}
    {ThresholdBasedBumpingConsensusContender<k, threshold_size>(N, 2.20).run();}
}

template <size_t k>
struct ThresholdBasedBumpingConsensusContenderRunner {
    void operator() (size_t N) const {
        constexpr size_t firstSize = tlx::integer_log2_ceil(k) >= 6 ? tlx::integer_log2_ceil(k) - 4 : 2;
        dispatchOverloadsConsensus<k, firstSize>(N);
        dispatchOverloadsConsensus<k, firstSize + 1>(N);
        dispatchOverloadsConsensus<k, firstSize + 2>(N);
        dispatchOverloadsConsensus<k, firstSize + 3>(N);
        dispatchOverloadsConsensus<k, firstSize + 4>(N);
        dispatchOverloadsConsensus<k, firstSize + 5>(N);
        dispatchOverloadsConsensus<k, firstSize + 6>(N);
        dispatchOverloadsConsensus<k, firstSize + 7>(N);
    }
};

void thresholdBasedBumpingConsensusContenderRunner(size_t N, size_t k) {
    dispatchK<ThresholdBasedBumpingConsensusContenderRunner>(k, N);
}
