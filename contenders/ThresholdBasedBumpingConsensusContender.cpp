#include "ThresholdBasedBumpingConsensusContender.hpp"
#include "DispatchK.h"

template <size_t k, size_t threshold_size>
void dispatchOverloadsConsensus(size_t N) {
    {ThresholdBasedBumpingConsensusContender<k, 1.05, threshold_size>(N).run();}
    {ThresholdBasedBumpingConsensusContender<k, 1.10, threshold_size>(N).run();}
    {ThresholdBasedBumpingConsensusContender<k, 1.20, threshold_size>(N).run();}
    {ThresholdBasedBumpingConsensusContender<k, 1.30, threshold_size>(N).run();}
    {ThresholdBasedBumpingConsensusContender<k, 1.40, threshold_size>(N).run();}
    {ThresholdBasedBumpingConsensusContender<k, 1.50, threshold_size>(N).run();}
    {ThresholdBasedBumpingConsensusContender<k, 1.60, threshold_size>(N).run();}
    {ThresholdBasedBumpingConsensusContender<k, 1.80, threshold_size>(N).run();}
    {ThresholdBasedBumpingConsensusContender<k, 2.00, threshold_size>(N).run();}
    {ThresholdBasedBumpingConsensusContender<k, 2.20, threshold_size>(N).run();}
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
