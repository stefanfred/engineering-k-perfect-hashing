#pragma once

#include <ThresholdBasedBumpingConsensus.hpp>
#include "Contender.h"

class ThresholdBasedBumpingConsensusContender : public Contender {
    public:
        kphf::ThresholdBasedBumpingConsensus::ThresholdBasedBumpingConsensus kphf;
        double overload;
        size_t thresholdSize;

        ThresholdBasedBumpingConsensusContender(size_t N, size_t k, double overload, size_t thresholdSize)
                : Contender(N, k, 1.0), overload(overload), thresholdSize(thresholdSize) {
        }

        std::string name() override {
            return std::string("ThresholdBasedBumpingConsensus")
                    + " overload=" + std::to_string(overload)
                    + " thresholdSize=" + std::to_string(thresholdSize);
        }

        void construct(const std::vector<std::string> &keys) override {
            // TODO: Directly take strings in constructor as well
            std::vector<Hash128> keysHashed;
            keysHashed.reserve(keys.size());
            for (auto &key : keys) {
                keysHashed.emplace_back(Hash128(key));
            }
            kphf = kphf::ThresholdBasedBumpingConsensus::ThresholdBasedBumpingConsensus(k_contender, keysHashed, overload, thresholdSize);
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

void thresholdBasedBumpingConsensusContenderRunner(size_t N, size_t k);
