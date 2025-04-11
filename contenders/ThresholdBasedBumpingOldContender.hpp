#pragma once

#include <consensus/BumpedKPerfectHashFunction.h>
#include "Contender.h"

template <size_t k>
class ThresholdBasedBumpingOldContender : public Contender {
    public:
        consensus::BumpedKPerfectHashFunction<k> *kphf = nullptr;

        explicit ThresholdBasedBumpingOldContender(size_t N)
                : Contender(N - N % k, k, 1.0) {
            // The technique calculates a 1-perfect hash function if N is not a multiple of k, so round down
        }

        ~ThresholdBasedBumpingOldContender() override {
            free(kphf);
        }

        std::string name() override {
            return std::string("ThresholdBasedBumpingOld");
        }

        void construct(const std::vector<std::string> &keys) override {
            // Other competitors hash keys internally. This competitor was
            // never meant to be used as standalone, so we have to do it from the outside.
            std::vector<uint64_t> hashedKeys;
            for (const auto &key : keys) {
                hashedKeys.emplace_back(bytehamster::util::MurmurHash64(key));
            }
            kphf = new consensus::BumpedKPerfectHashFunction<k>(std::span(hashedKeys));
        }

        size_t sizeBits() override {
            return kphf->getBits();
        }

        void performQueries(const std::span<std::string> keys) override {
            auto x = [&] (std::string &key) {
                return kphf->operator()(bytehamster::util::MurmurHash64(key));
            };
            doPerformQueries(keys, x);
        }

        void performTest(const std::span<std::string> keys) override {
            auto x = [&] (std::string &key) {
                return kphf->operator()(bytehamster::util::MurmurHash64(key));
            };
            doPerformTest(keys, x);
        }
};

void thresholdBasedBumpingOldContenderRunner(size_t N, size_t k);
