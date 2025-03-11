#pragma once

#include <cstdint>
#include <chrono>
#include <random>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <functional>
#include <vector>
#include <utility>
#include <csignal>
#include <memory>
#include <format>

#include <hash128.hpp>
#include "stat.hpp"

using Clock = std::chrono::steady_clock;

template<typename T>
static inline T &&black_box(T &&value) {
	asm volatile("" : "+rm" (value) : : "memory");
	return std::forward<T>(value);
}

template<typename F, typename ...Args>
std::invoke_result_t<F, Args...>
  measure(stat::Statistics::Builder &builder,
          F &&f, Args &&...args) {
	Clock::time_point start = Clock::now();
	auto result =
	  black_box(std::invoke(std::forward<F>(f), std::forward<Args>(black_box(std::forward<Args>(args)))...));
	builder.add(std::chrono::duration<double>(Clock::now() - start).count());
	return result;
}

template<typename F, typename R>
void measure_multi(stat::Statistics::Builder &builder, F &&f, R &&range) {
	uint64_t count = 0;
	Clock::time_point start = Clock::now();
	for (auto &&item: range) {
		black_box(std::invoke(std::forward<F>(f), black_box(item)));
		count++;
	}
	builder.add(std::chrono::duration<double>(Clock::now() - start).count()
	             / count);
}

class bad_hash: public std::runtime_error {
	using std::runtime_error::runtime_error;
};

class Benchmark {
public:
	static size_t numKeys;
	virtual void iteration() = 0;
	virtual void dump() = 0;
	virtual ~Benchmark() = default;
};
size_t Benchmark::numKeys = 1'000'000;

template<typename Builder>
class TestAndBenchmark: public Benchmark {
private:
	std::mt19937_64 rng;
	stat::Statistics::Builder cons_stats, query_stats, size_stats;
	Builder builder;

public:
	TestAndBenchmark(Builder &&builder): builder(builder) {}

	void iteration() override {
		std::vector<Hash128> items(numKeys);
		for (auto &[hi,lo]: items) hi = rng(), lo = rng();

		const typename Builder::Hash hash =
		  measure(cons_stats, builder, items);
		size_stats.add(hash.count_bits());

		uint64_t k = builder.get_k();
		std::vector<uint64_t> counts((numKeys + k - 1) / k);
#ifndef NDEBUG
		std::vector<std::vector<Hash128>> inv((numKeys + k - 1) / k);
		for (auto &v: inv) v.reserve(k);
#endif
		// TODO: Testing (with counting) and measuring the actual throughput should be two independent steps
		measure_multi(query_stats, [&](Hash128 item) {
			uint64_t r = hash(item);
			if (r >= counts.size()) throw bad_hash("hash value out of bounds");
			counts[r]++;
#ifndef NDEBUG
			inv[r].push_back(item);
#endif
			return r;
		}, items);

		for (size_t i = 0; i < counts.size(); i++) {
			if (counts[i] > k) {
#ifndef NDEBUG
				std::raise(SIGTRAP);
				for (Hash128 item: inv[i]) {
					black_box(hash(item));
				}
#endif
				throw bad_hash("hash not k-perfect");
			}
		}
	}

	void run() {
		for (unsigned int i = 0; i < 10; i++) {
			iteration();
		}
		dump();
	}

	void dump() override {
		auto construction = cons_stats.build();
		auto query = query_stats.build();
		auto size = size_stats.build();
		std::string name = builder.name();
		uint64_t k = builder.get_k();
		std::vector<std::pair<std::string, std::string>> meta = builder.meta();
		std::cout << name << "(k=" << k;
		for (auto [k,v]: meta) std::cout << ", " << k << "=" << v;
		std::cout << "):\n";
		std::cout << "  Construction time [s]: " << construction << "\n";
		std::cout << "  Query time [s]:        " << query << "\n";
		std::cout << "  Size [bit]:            " << size << "\n";
		std::cout << "RESULT"
			<< "\talgo=" << name
			<< "\tn=" << numKeys
			<< "\tk=" << k
			<< "\tconstruction=" << construction.mean
			<< "\tquery=" << query.mean
			<< "\tsize=" << size.mean;
		for (auto [k,v]: meta) std::cout << "\t" << k << "=" << v;
		std::cout << "\n";
	}
};
