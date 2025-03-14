#pragma once

#include <cstdint>
#include <chrono>
#include <random>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

#include <hash128.hpp>
#include "stat.hpp"

#define DO_NOT_OPTIMIZE(value) asm volatile ("" : : "r,m"(value) : "memory")
using Clock = std::chrono::steady_clock;

class bad_hash: public std::runtime_error {
	using std::runtime_error::runtime_error;
};

size_t numKeys = 1'000'000, numQueries = 1'000'000;

template<typename Builder>
class TestAndBenchmark {
	stat::Statistics::Builder cons_stats, query_stats, size_stats;
	Builder builder;
public:
	explicit TestAndBenchmark(Builder &&builder)
		: builder(builder) {
	}

	void iteration() {
		std::mt19937_64 rng(std::random_device{}());
		std::vector<Hash128> items(numKeys);
		for (auto &[hi,lo] : items) {
			hi = rng();
			lo = rng();
		}

		sleep(1); // Cooldown
		Clock::time_point start = Clock::now();
		auto hash = builder(items);
		cons_stats.add(std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start).count());
		size_stats.add(1.0 * hash.count_bits() / items.size());

		{
			std::vector<uint64_t> counts((numKeys + builder.get_k() - 1) / builder.get_k());
			for (auto &&item : items) {
				uint64_t r = hash(item);
				if (r >= counts.size()) throw bad_hash("hash value out of bounds");
				if (counts[r] >= builder.get_k()) throw bad_hash("hash not k-perfect");
				counts[r]++;
			}
		}

		std::vector<Hash128> queries(numQueries);
		queries.erase(std::ranges::sample(items, queries.begin(), numQueries, rng), queries.end());
		std::ranges::shuffle(queries, rng);

		sleep(1); // Cooldown
		start = Clock::now();
		for (auto &&item : queries) {
			size_t result = hash(item);
			DO_NOT_OPTIMIZE(result);
		}
		query_stats.add(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count() / queries.size());
	}

	void run() {
		std::cout << "Running " << builder.name();
		for (auto [k,v]: builder.meta()) std::cout << "\t" << k << "=" << v;
		std::cout << std::endl;
		for (unsigned int i = 0; i < 5; i++) {
			iteration();
		}
		dump();
	}

	void dump() {
		auto construction = cons_stats.build();
		auto query = query_stats.build();
		auto size = size_stats.build();

		std::cout << builder.name() << "(k=" << builder.get_k();
		for (auto [k,v]: builder.meta()) std::cout << ", " << k << "=" << v;
		std::cout << "):\n";
		std::cout << "  Construction time [ms]: " << construction << "\n";
		std::cout << "  Query time [ns]:        " << query << "\n";
		std::cout << "  Size [bit/key]:         " << size << "\n";
		std::cout << "RESULT"
			<< "\talgo=" << builder.name()
			<< "\tn=" << numKeys
			<< "\tk=" << builder.get_k()
			<< "\tconstructionMs=" << construction.mean
			<< "\tqueryNs=" << query.mean
			<< "\tsizePerKey=" << size.mean;
		for (auto [k,v]: builder.meta()) std::cout << "\t" << k << "=" << v;
		std::cout << std::endl;
	}
};
