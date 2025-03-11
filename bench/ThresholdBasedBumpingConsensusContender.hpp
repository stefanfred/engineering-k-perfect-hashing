#pragma once

#include "bench.hpp"
#include <ThresholdBasedBumpingConsensus.hpp>
#include "HashDisplaceContender.hpp"

namespace ThresholdBasedBumpingConsensusContender {

template<typename Builder>
struct ThresholdBasedBumpingConsensusContender {
private:
	uint64_t k;
	double overload;
	int threshold_size;
	Builder phf_builder;
	kphf::ThresholdBasedBumpingConsensus::ThresholdBasedBumpingConsensusBuilder<Builder> builder;

public:
	using Hash = kphf::ThresholdBasedBumpingConsensus::ThresholdBasedBumpingConsensus<typename Builder::Hash>;

	ThresholdBasedBumpingConsensusContender(uint64_t k, double overload, int threshold_size, Builder &&phf_builder):
	  k(k), overload(overload), threshold_size(threshold_size),
	  phf_builder(phf_builder),
	  builder(k, overload, threshold_size, std::move(phf_builder)) {}

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
			{ "phf", phf_builder.name() },
		};
		for (auto &[name, value]: phf_builder.meta()) {
			res.emplace_back("phf_" + name, value);
		}
		return res;
	}
};

using H = HashDisplaceContender::HashDisplaceContender<kphf::HashDisplace::OptimalBucketFunction, kphf::HashDisplace::RiceEncoding>;

H hash_displace(1, 6, 0.95);

void f(Benchmarks &bench, uint64_t k,
  std::initializer_list<uint64_t> threshold_sizes,
  std::initializer_list<double> overloads) {
	for (uint64_t threshold_size: threshold_sizes) {
		for (double overload: overloads) {
			bench.add(TestAndBenchmark(ThresholdBasedBumpingConsensusContender(k, overload, threshold_size, H(hash_displace))));
			//std::cerr << '.';
		}
		//std::cerr << '\n';
	}
}

void benchmark(Benchmarks &bench) {
	//f<10, 3, 4, 5>(bench);
	//f<100, 6, 7, 8>(bench);
	//f<1000, 8, 9, 10, 11>(bench);
	//f<100, 7,8,9>(bench);
	//f(bench, 100, {6});
	f(bench, 10, {3,4,5}, {1.85,1.90,1.95,2.00,2.05});
	f(bench, 100, {5,6,7}, {1.30,1.35,1.40});
	//f(bench, 1000, {8,9,10,11});
}
};
