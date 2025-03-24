#include "ThresholdBasedBumpingContender.hpp"

#include <gcem.hpp>
#include "DispatchK.h"


template <size_t k, size_t xFrom, size_t xTo>
void dispatchX(size_t N) {
    if constexpr (xFrom <= xTo) {
        {ThresholdBasedBumpingContender<k, 1.30, xFrom>(N).run();}
    }
    if constexpr (xFrom <= xTo) {
        dispatchX<k, xFrom + 1, xTo>(N);
    }
}

template <size_t k>
struct ThresholdBasedBumpingContenderRunner {
    void operator() (size_t N) const {
        // TODO: Where does this come from?
        constexpr size_t x = std::max(2ul, static_cast<size_t>(gcem::log2(2 * std::numbers::pi * k) / 2));
        dispatchX<k, x, x + 8>(N);
    }
};

void thresholdBasedBumpingContenderRunner(size_t N, size_t k) {
    dispatchK<ThresholdBasedBumpingContenderRunner>(k, N);
}
