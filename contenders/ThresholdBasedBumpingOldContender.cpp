#include "ThresholdBasedBumpingOldContender.hpp"

#include "DispatchK.h"

template <size_t k>
struct ThresholdBasedBumpingOldContenderRunner {
    void operator() (size_t N) const {
        ThresholdBasedBumpingOldContender<k>(N).run();
    }
};

void thresholdBasedBumpingOldContenderRunner(size_t N, size_t k) {
    dispatchK<ThresholdBasedBumpingOldContenderRunner>(k, N);
}
