#include "ThresholdBasedBumpingContender.hpp"

#include <gcem.hpp>
#include "DispatchK.h"

template <size_t k, size_t x, double overloadFrom, double overloadStep, int overloadStepCount>
void dispatchOverload(size_t N) {
	if constexpr (overloadFrom > 1.0) {
		ThresholdBasedBumpingContender<k, overloadFrom, x>(N).run();
	}
	if constexpr (overloadStepCount > 0) {
		dispatchOverload<k, x, overloadFrom + overloadStep, overloadStep, overloadStepCount - 1>(N);
	}
}


template <size_t k, size_t xFrom, size_t xTo>
void dispatchX(size_t N) {
    if constexpr (xFrom <= xTo) {
		constexpr double overload = 6/(gcem::sqrt(k) + gcem::log2(k));
		dispatchOverload<k, xFrom, 1 + overload * 0.95, overload * 0.025, 4>(N);
    }
    if constexpr (xFrom <= xTo) {
        dispatchX<k, xFrom + 1, xTo>(N);
    }
}

template <size_t k>
struct ThresholdBasedBumpingContenderRunner {
    void operator() (size_t N) const {
        // TODO: Where does this come from?
        constexpr size_t x = std::max(2ul, static_cast<size_t>(2 * gcem::log2(2 * std::numbers::pi * k) / 2));
        dispatchX<k, x, x + 8>(N);
    }
};

void thresholdBasedBumpingContenderRunner(size_t N, size_t k) {
    dispatchK<ThresholdBasedBumpingContenderRunner>(k, N);
}
