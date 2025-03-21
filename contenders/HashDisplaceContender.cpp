#include "HashDisplaceContender.hpp"

#include <OptimalBucketFunction.hpp>
#include <CompactEncoding.hpp>
#include <RiceEncoding.hpp>

template<typename Encoding>
void dispatchHashDisplaceBucketSize(size_t N, size_t k, size_t bucket_size) {
    { HashDisplaceContender<kphf::HashDisplace::OptimalBucketFunction, Encoding>(N, k, bucket_size).run(); }
}

template<typename ...Encoding>
void goDispatchHashDisplaceEncoder(size_t N, size_t k, size_t bucket_size) {
    (dispatchHashDisplaceBucketSize<Encoding>(N, k, bucket_size), ...);
}

void dispatchHashDisplaceEncoder(size_t N, size_t k, size_t bucket_size) {
    goDispatchHashDisplaceEncoder<kphf::HashDisplace::CompactEncoding, kphf::HashDisplace::RiceEncoding>(N, k, bucket_size);
}

void hashDisplaceContenderRunner(size_t N, size_t k) {
    if (k == 0) {
        throw std::invalid_argument("k must be greater than 0");
    }
    /* There is no good reason for this formula. */
    double b = sqrt(k) * log2(2 * std::numbers::pi * k) * 1.25;
    for (int i = 1; i <= 10; i++) {
        dispatchHashDisplaceEncoder(N, k, size_t(b * i / 10));
    }
}
