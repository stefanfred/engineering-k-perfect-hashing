#pragma once

#include <PaCHash.hpp>
#include "Contender.h"

class PaCHashContender : public Contender {
    public:
        kphf::PaCHash::PaCHash kphf;
        double a;

        PaCHashContender(size_t N, size_t k, double a)
                : Contender(N, k, 1.0), a(a) {
        }

        std::string name() override {
            return std::string("PaCHash")
                    + " a=" + std::to_string(a);
        }

        void construct(const std::vector<std::string> &keys) override {
            // TODO: Directly take strings in constructor as well
            std::vector<Hash128> keysHashed;
            keysHashed.reserve(keys.size());
            for (auto &key : keys) {
                keysHashed.emplace_back(Hash128(key));
            }
            kphf = kphf::PaCHash::PaCHash(k_contender, a, keysHashed);
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

void paCHashContenderRunner(size_t N, size_t k);
