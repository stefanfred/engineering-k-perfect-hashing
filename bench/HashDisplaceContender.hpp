#pragma once

#include <HashDisplace.hpp>
#include <OptimalBucketFunction.hpp>
#include <CompactEncoding.hpp>
#include <RiceEncoding.hpp>
#include "bench.hpp"

namespace HashDisplaceContender {

template<typename BucketFunction, typename Encoding>
struct HashDisplaceContender {
private:
	uint64_t k, bucket_size;
	double load_factor;
	kphf::HashDisplace::HashDisplaceBuilder<BucketFunction, Encoding> builder;

public:
	using Hash = kphf::HashDisplace::HashDisplace<BucketFunction, Encoding>;

	HashDisplaceContender(uint64_t k, uint64_t bucket_size, double load_factor = 1.0):
	  k(k), bucket_size(bucket_size), load_factor(load_factor), builder(k) {}

	Hash operator()(const std::vector<Hash128> &items) const {
		return builder(items, bucket_size, load_factor);
	}

	uint64_t get_k() const { return k; }

	static std::string name() {
		return "HashDisplace";
	}

	std::vector<std::pair<std::string, std::string>> meta() const {
		return std::vector<std::pair<std::string, std::string>> {
			{ "bucket_function", BucketFunction::name() },
			{ "encoding", Encoding::name() },
			{ "bucket_size", std::to_string(bucket_size) },
			{ "load_factor", std::to_string(load_factor) },
		};
	}
};

template<typename Encoding>
void h(size_t k, size_t bucket_size) {
	TestAndBenchmark(HashDisplaceContender<kphf::HashDisplace::OptimalBucketFunction, Encoding>(k, bucket_size)).run();
}

template<typename ...Encoding>
void g(size_t k, size_t bucket_size) {
	(h<Encoding>(k, bucket_size), ...);
}

void f(size_t k, size_t bucket_size) {
	g<kphf::HashDisplace::CompactEncoding, kphf::HashDisplace::RiceEncoding>(k, bucket_size);
}

void benchmark(size_t k) {
	double b = sqrt(k) * log2(2 * std::numbers::pi * k) * 1.25;
	for (int i = 1; i <= 10; i++) f(k, size_t(b * i / 10));
/*
	f(1000, {30,40,50,60,70,80,90,100,150,200,250,300,350,400});
	f(100, {10,20,30,40,50,60,70,80,90,100,110,120});
	f(10, {5,10,15,20,25,30});
*/
}

}
