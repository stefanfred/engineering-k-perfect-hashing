#include "ThresholdBasedBumpingContender.hpp"

#include <gcem.hpp>
#include "DispatchK.h"

template <size_t k, size_t xFrom, size_t xTo>
void dispatchX(size_t N) {
    if constexpr (xFrom <= xTo) {
		constexpr double targetOverload = 6/(gcem::sqrt(k) + gcem::log2(k));
        for (size_t i = 0; i < 4; i++) {
        	// TODO: Is "1+" really correct here? Consensus variant doesn't have it.
        	double overload = 1 + targetOverload * 0.95 + i * targetOverload * 0.025;
			ThresholdBasedBumpingContender<k, xFrom, kphf::ThresholdBasedBumping::DummyFilter>(N, overload).run();
			ThresholdBasedBumpingContender<k, xFrom, kphf::ThresholdBasedBumping::RibbonFilter>(N, overload).run();
        }
    }
    if constexpr (xFrom <= xTo) {
        dispatchX<k, xFrom + 1, xTo>(N);
    }
}

template <size_t k>
struct ThresholdBasedBumpingContenderRunner {
    void operator() (size_t N) const {
        // TODO: Where does this come from?
        constexpr size_t x = std::max(1ul, static_cast<size_t>(gcem::log2(2 * std::numbers::pi * k) / 2));
        dispatchX<k, x, x + 4>(N);
    }
};

void thresholdBasedBumpingContenderRunner(size_t N, size_t k) {
    dispatchK<ThresholdBasedBumpingContenderRunner>(k, N);
}
