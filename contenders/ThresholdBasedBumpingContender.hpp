#pragma once

#include <ThresholdBasedBumping.hpp>
#include "Contender.h"

template <size_t k, double overload, int thresholdSizeHalfbits = bytehamster::util::ceillog2(k),
        typename Filter = kphf::ThresholdBasedBumping::RibbonFilter>
class ThresholdBasedBumpingContender : public Contender {
    public:
        using ThresholdBasedBumping = kphf::ThresholdBasedBumping::ThresholdBasedBumping<k, overload, thresholdSizeHalfbits, Filter>;
        ThresholdBasedBumping kphf;

        explicit ThresholdBasedBumpingContender(size_t N) : Contender(N, k, 1.0) {
        }

        std::string name() override {
            return std::string("ThresholdBasedBumping")
                    + " overload=" + std::to_string(overload)
                    + " thresholdSize=" + std::to_string(0.5 * thresholdSizeHalfbits)
                    + " filter=" + Filter::name();
        }

        void construct(const std::vector<std::string> &keys) override {
            // TODO: Directly take strings in constructor as well
            std::vector<Hash128> keysHashed;
            keysHashed.reserve(keys.size());
            for (auto &key : keys) {
                keysHashed.emplace_back(Hash128(key));
            }
            kphf = ThresholdBasedBumping(keysHashed);
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

void thresholdBasedBumpingContenderRunner(size_t N, size_t k);
