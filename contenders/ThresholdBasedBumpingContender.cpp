#include "ThresholdBasedBumpingContender.hpp"

#include <tlx/math/integer_log2.hpp>
#include "DispatchK.h"

template <size_t k, size_t thresholdSize>
void dispatchFilter(size_t N, double overload) {
    ThresholdBasedBumpingContender<k, thresholdSize, kphf::ThresholdBasedBumping::DummyFilter>(N, overload).run();
    ThresholdBasedBumpingContender<k, thresholdSize, kphf::ThresholdBasedBumping::RibbonFilter>(N, overload).run();
}

template <size_t k, size_t thresholdSize>
void dispatchOverloads(size_t N) {
    dispatchFilter<k, thresholdSize>(N, 1.05);
    dispatchFilter<k, thresholdSize>(N, 1.10);
    dispatchFilter<k, thresholdSize>(N, 1.20);
    dispatchFilter<k, thresholdSize>(N, 1.30);
    dispatchFilter<k, thresholdSize>(N, 1.40);
    dispatchFilter<k, thresholdSize>(N, 1.50);
    dispatchFilter<k, thresholdSize>(N, 1.60);
    dispatchFilter<k, thresholdSize>(N, 1.80);
    dispatchFilter<k, thresholdSize>(N, 2.00);
    dispatchFilter<k, thresholdSize>(N, 2.20);
}

template <size_t k>
struct ThresholdBasedBumpingContenderRunner {
    void operator() (size_t N) const {
        constexpr size_t firstSize = tlx::integer_log2_ceil(k) >= 5 ? tlx::integer_log2_ceil(k) - 4 : 1;
        dispatchOverloads<k, firstSize>(N);
        dispatchOverloads<k, firstSize + 1>(N);
        dispatchOverloads<k, firstSize + 2>(N);
        dispatchOverloads<k, firstSize + 3>(N);
        dispatchOverloads<k, firstSize + 4>(N);
        dispatchOverloads<k, firstSize + 5>(N);
        dispatchOverloads<k, firstSize + 6>(N);
        dispatchOverloads<k, firstSize + 7>(N);
    }
};

void thresholdBasedBumpingContenderRunner(size_t N, size_t k) {
    dispatchK<ThresholdBasedBumpingContenderRunner>(k, N);
}
