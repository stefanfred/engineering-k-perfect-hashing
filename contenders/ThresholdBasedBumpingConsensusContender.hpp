#pragma once

#include <ThresholdBasedBumpingConsensus.hpp>
#include "Contender.h"

template <uint64_t k, int threshold_size>
class ThresholdBasedBumpingConsensusContender : public Contender {
    public:
        using kphf_t = kphf::ThresholdBasedBumpingConsensus::ThresholdBasedBumpingConsensus<k, threshold_size>;
        kphf_t kphf;
        double overload;

        explicit ThresholdBasedBumpingConsensusContender(size_t N, double overload)
                : Contender(N, k, 1.0), overload(overload) {
        }

        std::string name() override {
            return std::string("ThresholdBasedBumpingConsensus")
                    + " overload=" + std::to_string(overload)
                    + " thresholdSize=" + std::to_string(threshold_size);
        }

        void construct(const std::vector<std::string> &keys) override {
            kphf = kphf_t(keys, overload);
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

void thresholdBasedBumpingConsensusContenderRunner(size_t N, size_t k);
