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
                    + " encoding=" + typeid(Encoding).name(); // TODO: Support calling Encoding::name()
        }

        void construct(const std::vector<std::string> &keys) override {
            // TODO: Directly take strings in constructor as well
            std::vector<Hash128> keysHashed;
            keysHashed.reserve(keys.size());
            for (auto &key : keys) {
                keysHashed.emplace_back(Hash128(key));
            }
            kphf = kphf::HashDisplace::HashDisplace<k, BucketFunction, Encoding>(keysHashed, bucketSize);
        }

        size_t sizeBits() override {
            return kphf.count_bits();
        }

        void performQueries(const std::span<std::string> keys) override {
            auto x = [&] (std::string &key) {
                return kphf(Hash128(key));
            };
            doPerformQueries(keys, x);
        }

        void performTest(const std::span<std::string> keys) override {
            auto x = [&] (std::string &key) {
                return kphf(Hash128(key));
            };
            doPerformTest(keys, x);
        }
};

void hashDisplaceContenderRunner(size_t N, size_t k);
