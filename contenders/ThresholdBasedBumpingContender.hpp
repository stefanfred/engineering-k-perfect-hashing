#pragma once

#include <ThresholdBasedBumping.hpp>
#include "Contender.h"

template <size_t k, int thresholdSize = bytehamster::util::ceillog2(k),
        typename Filter = kphf::ThresholdBasedBumping::RibbonFilter>
class ThresholdBasedBumpingContender : public Contender {
    public:
        using ThresholdBasedBumping = kphf::ThresholdBasedBumping::ThresholdBasedBumping<k, thresholdSize, Filter>;
        ThresholdBasedBumping kphf;
        double overload;

        explicit ThresholdBasedBumpingContender(size_t N, double overload)
                : Contender(N, k, 1.0), overload(overload) {
        }

        std::string name() override {
            return std::string("ThresholdBasedBumping")
                    + " overload=" + std::to_string(overload)
                    + " thresholdSize=" + std::to_string(thresholdSize)
                    + " filter=" + Filter::name();
        }

        void construct(const std::vector<std::string> &keys) override {
            kphf = ThresholdBasedBumping(keys, overload);
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

void thresholdBasedBumpingContenderRunner(size_t N, size_t k);
