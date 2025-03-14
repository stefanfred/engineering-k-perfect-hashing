#pragma once

#include <iostream>
#include <cmath>
#include <utility>
#include <vector>

namespace kphf {
namespace HashDisplace {
namespace optimal_bucket_function {

constexpr double poission_pmf(int i, double lambda) {
    return pow(M_E, i * log(lambda) - lambda - lgamma(i + 1.0));
}

constexpr double expected_truncated(int maxEx, double lambda, bool weighted) {
    double exp = 0.0;
    for (int i = 0; i < maxEx; ++i) {
        exp += poission_pmf(i, lambda) * (weighted ? i : 1.0);
    }
    return exp;
}


constexpr std::pair<double, double> alphaAndRelSize(int k, double lambda) {
    double alphaNotFull = expected_truncated(k, lambda, true) / k;
    double prob = expected_truncated(k, lambda, false);
    return {1.0 - prob + alphaNotFull, -log(prob)};
}

constexpr void buildRec(std::vector<std::pair<double, double>> &samples, double deltaX, double llamb, double rlamb, int k) {
    std::pair<double, double> left = alphaAndRelSize(k, llamb);
    std::pair<double, double> right = alphaAndRelSize(k, rlamb);

    //double area = (left.second + right.second) / 2.0 * (right.first - left.first);
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

constexpr std::vector<uint64_t> getBucketFunctionFulcrums(int k, int fulcs) {
    std::vector<std::pair<double, double>> samples;
    double leftLimit = k * 0.0001;
    double rightLimit = k * 10;
    buildRec(samples, 0.1 / fulcs, leftLimit, rightLimit, k);

    std::vector<std::pair<double, double>> integral;
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

uint64_t queryFulcrums(const std::vector<uint64_t> &fulcs, uint64_t hash, uint64_t bucketCount) {
    __uint128_t z = hash * __uint128_t(fulcs.size() - 1);
    uint64_t index = z >> 64;
    uint64_t part = z;
    uint64_t v1 = (uint64_t((fulcs[index + 0] * __uint128_t(bucketCount)) >> 32) * __uint128_t(uint64_t(-1) - part)) >> 64;
    uint64_t v2 = (uint64_t((fulcs[index + 1] * __uint128_t(bucketCount)) >> 32) * __uint128_t(part)) >> 64;
    return (v1 + v2) >> 32;
}

}

class OptimalBucketFunction {
private:
	std::shared_ptr<std::vector<uint64_t>> fulcrums;

public:
	class Instance {
	private:
		std::shared_ptr<std::vector<uint64_t>> fulcrums;
		uint64_t bucket_count, load_factor;
		Instance(std::shared_ptr<std::vector<uint64_t>> fulcrums,
		  uint64_t bucket_count, uint64_t load_factor):
			fulcrums(move(fulcrums)), bucket_count(bucket_count),
		    load_factor(load_factor) {}

	public:
		uint64_t operator()(uint64_t hash) const {
			if (load_factor != 0) {
				hash = ((unsigned __int128) hash * load_factor) >> 64;
			}
			return optimal_bucket_function::queryFulcrums(*fulcrums, hash, bucket_count);
		}

		size_t count_bits() const {
			return 8 * sizeof(*this);
		}

		friend class OptimalBucketFunction;
	};

	OptimalBucketFunction(uint64_t k):
	  fulcrums(std::make_shared<std::vector<uint64_t>>(optimal_bucket_function::getBucketFunctionFulcrums(k, 1000))) {
	}

	Instance operator()(uint64_t n, uint64_t nbuckets, double load_factor)
	  const {
		(void) n;
		uint64_t lf64;
		if (load_factor == 1.0) {
			lf64 = 0;
		} else {
			lf64 = uint64_t(ldexp(load_factor, 64));
			assert(lf64 > 0);

			auto is_ok = [&](uint64_t nb) -> bool {
				return optimal_bucket_function::queryFulcrums(*fulcrums, lf64-1, nb) < nbuckets;
			};
			uint64_t l = nbuckets, r;
			while (is_ok(r = 2*l)) l = r;
			while (l+1 != r) {
				uint64_t mid = std::midpoint(l, r);
				if (is_ok(mid)) l = mid;
				else r = mid;
			}
			nbuckets = l;
		}
		return Instance(fulcrums, nbuckets, lf64);
	}

	static std::string name() { return "optimal"; }
};

}
}
