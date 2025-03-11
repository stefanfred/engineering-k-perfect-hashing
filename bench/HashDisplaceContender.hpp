#pragma once

#include <HashDisplace.hpp>
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

constexpr size_t N = 10'000'000;

template<typename Encoding>
void h(Benchmarks &bench, size_t k, std::initializer_list<size_t> bucket_sizes) {
	for (size_t bucket_size: bucket_sizes) {
		bench.add(TestAndBenchmark(HashDisplaceContender<kphf::HashDisplace::OptimalBucketFunction, Encoding>(k, bucket_size), N));
	}
}

template<typename ...Encoding>
void g(Benchmarks &bench, size_t k, std::initializer_list<size_t> bucket_sizes) {
	(h<Encoding>(bench, k, bucket_sizes), ...);
}

void f(Benchmarks &bench, size_t k, std::initializer_list<size_t> bucket_sizes) {
	g<kphf::HashDisplace::CompactEncoding, kphf::HashDisplace::RiceEncoding>(bench, k, bucket_sizes);
}

void benchmark(Benchmarks &bench) {
	f(bench, 1000, {30,40,50,60,70,80,90,100,150,200,250,300,350,400});
	f(bench, 100, {10,20,30,40,50,60,70,80,90,100,110,120});
	f(bench, 10, {5,10,15,20,25,30});
}

}
