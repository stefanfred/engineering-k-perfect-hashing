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
	virtual void iteration() = 0;
	virtual void dump() = 0;
	virtual ~Benchmark() = default;
};

template<typename Builder>
class BenchmarkSample: public Benchmark {
private:
	std::mt19937_64 rng;
	stat::Statistics::Builder cons_stats, query_stats, size_stats;
	Builder builder;
	uint64_t n, sample;

public:
	BenchmarkSample(Builder &&builder, uint64_t n, uint64_t sample):
	  builder(builder), n(n), sample(sample) {}

	void iteration() override {
		std::vector<Hash128> items(n);
		for (auto &[hi,lo]: items) hi = rng(), lo = rng();

		const typename Builder::Hash hash =
		  measure(cons_stats, builder, items);
		size_stats.add(hash.count_bits());

		std::vector<Hash128> samples(sample);
		for (Hash128 &h: samples) {
			h = items[std::uniform_int_distribution<uint64_t>(0, n-1)(rng)];
		}
		measure_multi(query_stats, hash, samples);
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
			<< "\tn=" << n
			<< "\tk=" << k
			<< "\tconstruction=" << construction.mean
			<< "\tquery=" << query.mean
			<< "\tsize=" << size.mean;
		for (auto [k,v]: meta) std::cout << "\t" << k << "=" << v;
		std::cout << "\n";
	}
};

template<typename Builder>
class TestAndBenchmark: public Benchmark {
private:
	std::mt19937_64 rng;
	stat::Statistics::Builder cons_stats, query_stats, size_stats;
	Builder builder;
	uint64_t n;

public:
	TestAndBenchmark(Builder &&builder, uint64_t n): builder(builder), n(n) {}

	void iteration() override {
		std::vector<Hash128> items(n);
		for (auto &[hi,lo]: items) hi = rng(), lo = rng();

		const typename Builder::Hash hash =
		  measure(cons_stats, builder, items);
		size_stats.add(hash.count_bits());

		uint64_t k = builder.get_k();
		std::vector<uint64_t> counts((n + k - 1) / k);
#ifndef NDEBUG
		std::vector<std::vector<Hash128>> inv((n + k - 1) / k);
		for (auto &v: inv) v.reserve(k);
#endif
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
			<< "\tn=" << n
			<< "\tk=" << k
			<< "\tconstruction=" << construction.mean
			<< "\tquery=" << query.mean
			<< "\tsize=" << size.mean;
		for (auto [k,v]: meta) std::cout << "\t" << k << "=" << v;
		std::cout << "\n";
	}
};

class Benchmarks {
private:
	std::vector<std::unique_ptr<Benchmark>> benchmarks;

public:
	Benchmarks() {}

	template<typename B>
	void add(B &&bench) {
		benchmarks.push_back(std::make_unique<B>(bench));
	}

	void run(unsigned int runs, bool quiet = false) {
		{
			const std::time_t now = std::time(nullptr);
			std::cerr <<
			  std::put_time(std::localtime(&now), "Started at %T\n");
		}
		size_t total = runs * size(benchmarks);
		size_t done = 0;
		for (unsigned int i = 0; i < runs; i++) {
			for (auto &b: benchmarks) {
				std::cerr << std::format("  {:2}%\r", 100u * done++ / total);

				b->iteration();
			}
		}
		std::cerr << "      \r";
		if (!quiet) for (auto &b: benchmarks) b->dump();
	}
};
