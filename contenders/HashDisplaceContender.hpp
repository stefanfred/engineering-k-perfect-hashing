#pragma once

#include <HashDisplace.hpp>
#include "Contender.h"

template<size_t k, typename BucketFunction, typename Encoding>
class HashDisplaceContender : public Contender {
    public:
        kphf::HashDisplace::HashDisplace<k, BucketFunction, Encoding> kphf;
        size_t bucketSize;

        HashDisplaceContender(const size_t N, const size_t bucketSize)
                : Contender(N, k, 1.0), bucketSize(bucketSize) {
        }

        std::string name() override {
            return std::string("HashDisplace")
                    + " bucketSize=" + std::to_string(bucketSize)
                    + " bucketFunction=" + BucketFunction::name()
                    + " encoding=" + Encoding::name();
        }

        void construct(const std::vector<std::string> &keys) override {
            kphf = kphf::HashDisplace::HashDisplace<k, BucketFunction, Encoding>(keys, bucketSize);
        }

        size_t sizeBits() override {
            return kphf.count_bits();
        }

        void performQueries(const std::span<std::string> keys) override {
            auto x = [&] (std::string &key) {
                return kphf(key);
            };
            doPerformQueries(keys, x);
        }

        void performTest(const std::span<std::string> keys) override {
            auto x = [&] (std::string &key) {
                return kphf(key);
            };
            doPerformTest(keys, x);
        }
};

void hashDisplaceContenderRunner(size_t N, size_t k);
