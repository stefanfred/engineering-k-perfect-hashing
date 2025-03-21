#include "HashDisplaceContender.hpp"

#include <OptimalBucketFunction.hpp>
#include <CompactEncoding.hpp>
#include <RiceEncoding.hpp>
#include "DispatchK.h"

template <size_t k>
struct HashDisplaceContenderRunner {
    void operator() (size_t N) const {
        using namespace kphf::HashDisplace;
        // TODO There is no good reason for this formula
        const double b = sqrt(k) * log2(2 * std::numbers::pi * k) * 1.25;
        for (size_t i = 1; i <= 10; i++) {
            size_t bucketSize = (b * i) / 10;
            HashDisplaceContender<k, OptimalBucketFunction, CompactEncoding>(N, bucketSize).run();
            HashDisplaceContender<k, OptimalBucketFunction, RiceEncoding>(N, bucketSize).run();
        }
    }
};

void hashDisplaceContenderRunner(const size_t N, const size_t k) {
    if (k == 0) {
        throw std::invalid_argument("k must be greater than 0");
    }
    dispatchK<HashDisplaceContenderRunner>(k, N);
}
