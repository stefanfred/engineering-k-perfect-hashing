#pragma once

#include <cmath>
#include <utility>
#include <vector>
#include <gcem.hpp>

namespace kphf::HashDisplace {
namespace optimal_bucket_function {

template <size_t MAX>
constexpr std::array<double, MAX> build_lgamma_lookup() {
    std::array<double, MAX> lgamma_lookup;
    for (size_t i = 0; i < MAX; ++i) {
        lgamma_lookup[i] = gcem::lgamma(i);
    }
    return lgamma_lookup;
}

template <size_t MAX>
std::array<double, MAX> lgamma_lookup = build_lgamma_lookup<MAX>();

template <size_t MAX_I>
double poission_pmf(int i, double lambda, double loglambda) {
    return std::exp(i * loglambda - lambda - lgamma_lookup<MAX_I + 1>[i + 1]);
}

template <size_t maxEx>
std::pair<double, double> expected_truncated(double lambda) {
    const double loglambda = log(lambda);
    double expWeighted = 0.0;
    double expUnweighted = 0.0;
    for (int i = 0; i < maxEx; ++i) {
        double pmf = poission_pmf<maxEx>(i, lambda, loglambda);
        expWeighted += pmf * i;
        expUnweighted += pmf;
    }
    return std::make_pair(expWeighted, expUnweighted);
}

template <size_t k>
std::pair<double, double> alphaAndRelSize(double lambda) {
    auto [expWeighted, expUnweighted] = expected_truncated<k>(lambda);
    double alphaNotFull = expWeighted / k;
    double prob = expUnweighted;
    return {1.0 - prob + alphaNotFull, -log(prob)};
}

template <size_t k>
void buildRec(std::vector<std::pair<double, double> > &samples, double deltaX, double llamb, double rlamb) {
    std::pair<double, double> left = alphaAndRelSize<k>(llamb);
    std::pair<double, double> right = alphaAndRelSize<k>(rlamb);

    if (right.first - left.first < deltaX)
        return;
    double midPoint = (rlamb + llamb) / 2.0;
    std::pair<double, double> midSample = alphaAndRelSize<k>(midPoint);
    buildRec<k>(samples, deltaX, llamb, midPoint);
    if (std::isnormal(midSample.second)) {
        samples.push_back(midSample);
        buildRec<k>(samples, deltaX, midPoint, rlamb);
    }
}

template <size_t k>
std::vector<uint64_t> getBucketFunctionFulcrums(int fulcs) {
    std::vector<std::pair<double, double> > samples;
    double leftLimit = k * 0.0001;
    double rightLimit = k * 10;

    buildRec<k>(samples, 0.1 / fulcs, leftLimit, rightLimit);

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

    std::vector<uint64_t> fulcrums;
    int index = 0;
    fulcrums.push_back(0);
    for (int i = 1; i < fulcs - 1; ++i) {
        double x = i / (double) (fulcs - 1);
        while (integral[index].first < x) {
            index++;
        }
        fulcrums.push_back(integral[index].second * uint64_t(-1));
    }
    fulcrums.push_back(uint64_t(-1));

    return fulcrums;
}

inline uint64_t queryFulcrums(const std::vector<uint64_t> &fulcs, uint64_t hash, uint64_t bucketCount) {
    __uint128_t z = hash * __uint128_t(fulcs.size() - 1);
    uint64_t index = z >> 64;
    uint64_t part = z;
    uint64_t v1 = (uint64_t((fulcs[index + 0] * __uint128_t(bucketCount)) >> 32) * __uint128_t(uint64_t(-1) - part)) >> 64;
    uint64_t v2 = (uint64_t((fulcs[index + 1] * __uint128_t(bucketCount)) >> 32) * __uint128_t(part)) >> 64;
    return (v1 + v2) >> 32;
}
}

template <size_t k>
class OptimalBucketFunction {
    private:
        // TODO: Make static constexpr when compiling got faster
        std::vector<uint64_t> fulcrums = optimal_bucket_function::getBucketFunctionFulcrums<k>(1000);
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
                    return optimal_bucket_function::queryFulcrums(fulcrums, lf64 - 1, nb) < nbuckets;
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
