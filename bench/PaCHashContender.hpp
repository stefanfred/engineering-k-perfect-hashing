#pragma once

#include <PaCHash.hpp>
#include "bench.hpp"

namespace PaCHashContender {

struct PaCHashContender {
private:
	int k;
	double bucket_size;

public:
	using Hash = kphf::PaCHash::PaCHash;

	PaCHashContender(int k, double bucket_size): k(k), bucket_size(bucket_size) {}

	Hash operator()(const std::vector<Hash128> &items) const {
		return Hash(k, bucket_size, items);
	}

	uint64_t get_k() const { return k; }

	static std::string name() {
		return "PaCHash";
	}

	std::vector<std::pair<std::string, std::string>> meta() const {
		return std::vector<std::pair<std::string, std::string>> {
			{"bucket_size", std::to_string(bucket_size)},
		};
	}
};

void benchmark() {
	for (int i = 1; i < 200; i++) {
		TestAndBenchmark(PaCHashContender(10, i / 20.0)).run();
		TestAndBenchmark(PaCHashContender(100, i / 20.0)).run();
		TestAndBenchmark(PaCHashContender(1000, i / 20.0)).run();
	}
}

}
