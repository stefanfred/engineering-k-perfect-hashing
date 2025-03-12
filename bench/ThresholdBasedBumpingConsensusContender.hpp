#pragma once

#include "bench.hpp"
#include <ThresholdBasedBumpingConsensus.hpp>
#include "HashDisplaceContender.hpp"

namespace ThresholdBasedBumpingConsensusContender {

struct ThresholdBasedBumpingConsensusContender {
private:
	uint64_t k;
	double overload;
	int threshold_size;
	using FallbackPhf = HashDisplaceContender::HashDisplaceContender<kphf::HashDisplace::OptimalBucketFunction, kphf::HashDisplace::RiceEncoding>;
	FallbackPhf fallback_phf;
	kphf::ThresholdBasedBumpingConsensus::ThresholdBasedBumpingConsensusBuilder<FallbackPhf> builder;
public:
	using Hash = kphf::ThresholdBasedBumpingConsensus::ThresholdBasedBumpingConsensus<FallbackPhf::Hash>;

	ThresholdBasedBumpingConsensusContender(uint64_t k, double overload, int threshold_size):
	  k(k), overload(overload), threshold_size(threshold_size),
	  fallback_phf(1, 6, 0.95),
	  builder(k, overload, threshold_size, std::move(fallback_phf)) {}

	Hash operator()(const std::vector<Hash128> &items) const {
		return builder(items);
	}

	uint64_t get_k() const { return k; }

	static std::string name() {
		return "ThresholdBasedBumpingConsensus";
	}

	std::vector<std::pair<std::string, std::string>> meta() const {
		std::vector<std::pair<std::string, std::string>> res = {
			{ "overload", std::to_string(overload) },
			{ "threshold_size", std::to_string(threshold_size) },
			{ "phf", fallback_phf.name() },
		};
		for (auto &[name, value]: fallback_phf.meta()) {
			res.emplace_back("phf_" + name, value);
		}
		return res;
	}
};

void dispatchThresholdsAndOverloads(size_t numKeys, uint64_t k,
		std::initializer_list<uint64_t> threshold_sizes, std::initializer_list<double> overloads) {
	for (uint64_t threshold_size: threshold_sizes) {
		for (double overload: overloads) {
			TestAndBenchmark(numKeys, ThresholdBasedBumpingConsensusContender(k, overload, threshold_size)).run();
		}
	}
}

void benchmark(size_t numKeys) {
	//f<10, 3, 4, 5>(bench);
	//f<100, 6, 7, 8>(bench);
	//f<1000, 8, 9, 10, 11>(bench);
	//f<100, 7,8,9>(bench);
	//f(bench, 100, {6});
	dispatchThresholdsAndOverloads(numKeys, 10, {3,4,5}, {1.85,1.90,1.95,2.00,2.05});
	dispatchThresholdsAndOverloads(numKeys, 100, {5,6,7}, {1.30,1.35,1.40});
	//f(bench, 1000, {8,9,10,11});
}
};
