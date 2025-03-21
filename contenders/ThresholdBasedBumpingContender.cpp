#include "ThresholdBasedBumpingContender.hpp"

#include <gcem.hpp>
#include "DispatchK.h"

template <size_t k>
struct ThresholdBasedBumpingContenderRunner {
    void operator() (size_t N) const {
        constexpr size_t x = gcem::log2(2 * std::numbers::pi * k) / 2;
        {ThresholdBasedBumpingContender<k, 1.30, x>(N).run();}
    }
};

void thresholdBasedBumpingContenderRunner(size_t N, size_t k) {
    dispatchK<ThresholdBasedBumpingContenderRunner>(k, N);
}
