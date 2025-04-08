#pragma once

#include <kRecSplit.hpp>
#include "Contender.h"

template <size_t k, size_t leafSize>
class KRecSplitContender : public Contender {
    public:
    using kPHF = kphf::RecSplit::RecSplit<k, leafSize>;
    kPHF *kphf = nullptr;
    size_t bucketSize;

        KRecSplitContender(size_t N, size_t bucketSize)
                : Contender(N, k, 1.0), bucketSize(bucketSize) {
        }

        ~KRecSplitContender() override {
            delete kphf;
        }

        std::string name() override {
            return std::string("kRecSplit")
                    + " leafSize=" + std::to_string(leafSize)
                    + " bucketSize=" + std::to_string(bucketSize);
        }

        void construct(const std::vector<std::string> &keys) override {
            kphf = new kPHF(keys, bucketSize);
        }

        size_t sizeBits() override {
            return kphf->bitCount();
        }

        void performQueries(const std::span<std::string> keys) override {
            auto x = [&] (std::string &key) {
                return kphf->operator()(key);
            };
            doPerformQueries(keys, x);
        }

        void performTest(const std::span<std::string> keys) override {
            auto x = [&] (std::string &key) {
                return kphf->operator()(key);
            };
            doPerformTest(keys, x);
        }
};

void kRecSplitContenderRunner(size_t N, size_t k);
