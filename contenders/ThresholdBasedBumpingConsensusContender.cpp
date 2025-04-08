#include "ThresholdBasedBumpingConsensusContender.hpp"
#include "DispatchK.h"

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
        constexpr size_t x = log2(2 * std::numbers::pi * k) / 2;
        dispatchOverloadsConsensus<k, x>(N);
        dispatchOverloadsConsensus<k, x + 1>(N);
        dispatchOverloadsConsensus<k, x + 2>(N);
    }
};

void thresholdBasedBumpingConsensusContenderRunner(size_t N, size_t k) {
    dispatchK<ThresholdBasedBumpingConsensusContenderRunner>(k, N);
}
