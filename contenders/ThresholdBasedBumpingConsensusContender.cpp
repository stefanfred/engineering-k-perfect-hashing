#include "ThresholdBasedBumpingConsensusContender.hpp"

void dispatchThresholdsAndOverloadsConsensus(size_t N, uint64_t k, std::initializer_list<uint64_t> threshold_sizes,
                                    std::initializer_list<double> overloads) {
    for (uint64_t threshold_size: threshold_sizes) {
        for (double overload: overloads) {
            // TODO: These should probably be templates for better compile-time evaluation
            {ThresholdBasedBumpingConsensusContender(N, k, overload, threshold_size).run();}
        }
    }
}

void thresholdBasedBumpingConsensusContenderRunner(size_t N, size_t k) {
    size_t x = log2(2 * std::numbers::pi * k) / 2;
    dispatchThresholdsAndOverloadsConsensus(N, k, {x,x+1,x+2}, {1.05,1.10,1.20,1.30,1.40,1.50,1.60,1.80,2.00,2.20});
}
