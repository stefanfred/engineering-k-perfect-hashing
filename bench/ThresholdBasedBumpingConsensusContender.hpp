#pragma once

#include "bench.hpp"
#include <kphf/kphf.hpp>
#include <kphf/ThresholdBasedBumpingConsensus.hpp>
#include "HashDisplaceContender.hpp"

namespace ThresholdBasedBumpingConsensusContender {

template<uint64_t K, double OVERLOAD,
  uint64_t THRESHOLD_SIZE, typename Builder>
struct ThresholdBasedBumpingConsensusContender {
private:
	Builder phf_builder;

public:
	using Hash = kphf::ThresholdBasedBumpingConsensus::ThresholdBasedBumpingConsensus<K, OVERLOAD, THRESHOLD_SIZE, typename Builder::Hash>;

	ThresholdBasedBumpingConsensusContender(Builder &&phf_builder): phf_builder(std::move(phf_builder)) {}

	Hash operator()(const std::vector<Hash128> &items) const {
		return Hash(items, std::as_const(phf_builder));
	}

	uint64_t get_k() const { return K; }

	static std::string name() {
		return "ThresholdBasedBumpingConsensus";
	}

	std::vector<std::pair<std::string, std::string>> meta() const {
		std::vector<std::pair<std::string, std::string>> res = {
			{ "overload", std::to_string(OVERLOAD) },
			{ "threshold_size", std::to_string(THRESHOLD_SIZE) },
			{ "phf", phf_builder.name() },
		};
		for (auto &[name, value]: phf_builder.meta()) {
			res.emplace_back("phf_" + name, value);
		}
		return res;
	}
};

constexpr uint64_t N = 10'000'000;

using H = HashDisplaceContender::HashDisplaceContender<kphf::HashDisplace::OptimalBucketFunction, kphf::HashDisplace::RiceEncoding>;

template<uint64_t K, uint64_t THRESHOLD_SIZE, double OVERLOAD>
void i(Benchmarks &bench) {
	for (uint64_t bucket_size: {6}) {
		for (double load_factor: {0.95}) {
			bench.add(TestAndBenchmark(ThresholdBasedBumpingConsensusContender<K, OVERLOAD, THRESHOLD_SIZE, H>(H(1, bucket_size, load_factor)), N));
		}
	}
}

template<uint64_t K, uint64_t THRESHOLD_SIZE, double ...OVERLOAD>
void h(Benchmarks &bench) {
	(i<K, THRESHOLD_SIZE, OVERLOAD>(bench), ...);
}

template<uint64_t K, uint64_t ...THRESHOLD_SIZE>
void f(Benchmarks &bench) {
	//(h<K, THRESHOLD_SIZE, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00>(bench), ...);
	//(h<K, THRESHOLD_SIZE, 1.30, 1.35, 1.40>(bench), ...);
	(h<K, THRESHOLD_SIZE, 1.30>(bench), ...);
}

void benchmark(Benchmarks &bench) {
	//f<10, 3, 4, 5>(bench);
	//f<100, 6, 7, 8>(bench);
	//f<1000, 8, 9, 10, 11>(bench);
	//f<100, 7,8,9>(bench);
	f<100, 6>(bench);
}
};
