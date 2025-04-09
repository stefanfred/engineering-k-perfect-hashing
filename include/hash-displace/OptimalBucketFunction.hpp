#pragma once

#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <gcem.hpp>

namespace kphf::HashDisplace {
namespace optimal_bucket_function {
double poission_pmf(int i, double lambda) {
    return exp(i * log(lambda) - lambda - lgamma(i + 1.0));
}

double expected_truncated(int maxEx, double lambda, bool weighted) {
    double exp = 0.0;
    for (int i = 0; i < maxEx; ++i) {
        exp += poission_pmf(i, lambda) * (weighted ? i : 1.0);
    }
    return exp;
}

std::pair<double, double> alphaAndRelSize(int k, double lambda) {
    double alphaNotFull = expected_truncated(k, lambda, true) / k;
    double prob = expected_truncated(k, lambda, false);
    return {1.0 - prob + alphaNotFull, -log(prob)};
}

void buildRec(std::vector<std::pair<double, double> > &samples, double deltaX, double llamb, double rlamb, int k) {
    std::pair<double, double> left = alphaAndRelSize(k, llamb);
    std::pair<double, double> right = alphaAndRelSize(k, rlamb);

    if (right.first - left.first < deltaX)
        return;
    double midPoint = (rlamb + llamb) / 2.0;
    std::pair<double, double> midSample = alphaAndRelSize(k, midPoint);
    buildRec(samples, deltaX, llamb, midPoint, k);
    if (std::isnormal(midSample.second)) {
        samples.push_back(midSample);
        buildRec(samples, deltaX, midPoint, rlamb, k);
    }
}

template<size_t fulcs>
std::array<uint64_t, fulcs> getBucketFunctionFulcrums(int k) {
    std::vector<std::pair<double, double> > samples;
    double leftLimit = k * 0.0001;
    double rightLimit = k * 10;

    buildRec(samples, 0.1 / fulcs, leftLimit, rightLimit, k);

    std::vector<std::pair<double, double> > integral;
    std::pair<double, double> lastPair = {0, 0};
    integral.push_back(lastPair);
    for (auto pair: samples) {
        pair.second = std::max(lastPair.second, lastPair.second + pair.second * (pair.first - lastPair.first));
        integral.push_back(pair);
        lastPair = pair;
    }
    for (auto &pair: integral) {
        pair.second /= lastPair.second;
    }

    std::array<uint64_t, fulcs> fulcrums;
    int index = 0;
    fulcrums[0]=0;
    for (int i = 1; i < fulcs - 1; ++i) {
        double x = i / (double) (fulcs - 1);
        while (integral[index].first < x) {
            index++;
        }
        fulcrums[i]=integral[index].second * uint64_t(-1);
    }
    fulcrums[fulcs-1]=uint64_t(-1);

    return fulcrums;
}

template <size_t fulcs_n>
inline uint64_t queryFulcrums(const std::array<uint64_t, fulcs_n> &fulcs, uint64_t hash, uint64_t bucketCount) {
    __uint128_t z = hash * __uint128_t(fulcs.size() - 1);
    uint64_t index = z >> 64;
    uint64_t part = z;
    uint64_t v1 = (uint64_t((fulcs[index + 0] * __uint128_t(bucketCount)) >> 32) * __uint128_t(uint64_t(-1) - part)) >> 64;
    uint64_t v2 = (uint64_t((fulcs[index + 1] * __uint128_t(bucketCount)) >> 32) * __uint128_t(part)) >> 64;
    return (v1 + v2) >> 32;
}
}

template <size_t k, size_t fulcs_n = 1024>
class OptimalBucketFunction {
    private:
        // TODO: Make static constexpr when compiling got faster
        std::array<uint64_t, fulcs_n> fulcrums = optimal_bucket_function::getBucketFunctionFulcrums<fulcs_n>(k);
        uint64_t lf64;
        uint64_t nbuckets;
    public:
        OptimalBucketFunction(uint64_t nbuckets_, double loadFactor)
                : nbuckets(nbuckets_) {
            if (loadFactor == 1.0) {
                lf64 = 0;
            } else {
                lf64 = uint64_t(ldexp(loadFactor, 64));
                assert(lf64 > 0);

                auto is_ok = [&](uint64_t nb) -> bool {
                    return optimal_bucket_function::queryFulcrums<fulcs_n>(fulcrums, lf64 - 1, nb) < nbuckets;
                };
                uint64_t l = nbuckets, r;
                while (is_ok(r = 2 * l)) l = r;
                while (l + 1 != r) {
                    uint64_t mid = std::midpoint(l, r);
                    if (is_ok(mid)) l = mid;
                    else r = mid;
                }
                nbuckets = l;
            }
        }

        uint64_t operator()(uint64_t hash) const {
            if (lf64 != 0) {
                hash = ((unsigned __int128) hash * lf64) >> 64;
            }
            return optimal_bucket_function::queryFulcrums(fulcrums, hash, nbuckets);
        }

        static std::string name() {
            return "optimal";
        }
};
}
