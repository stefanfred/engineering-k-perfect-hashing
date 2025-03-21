#pragma once

#include <vector>
#include <bit>
#include <ranges>
#ifdef STATS
#include <map>
#endif
#include <SimpleRibbon.h>
#include <ips2ra.hpp>
#include <bytehamster/util/Function.h>

#include <hash128.hpp>
#include "EliasFano.hpp"

namespace kphf::PaCHash {

#ifdef STATS
struct Config {
	int k;
	double bucket_size;

	bool operator<(const Config &rhs) const {
		return k < rhs.k || (k == rhs.k && bucket_size < rhs.bucket_size);
	}
};

struct Stat {
	uint64_t ribbon_total = 0;
	uint64_t keys_total = 0;
};

std::map<Config, Stat> stats;

void dump_stats() {
	for (auto [conf, stat]: stats) {
		std::cout << "RESULT"
			<< "\tk=" << conf.k
			<< "\tbucket_size=" << conf.bucket_size
			<< "\tribbon_per_key="
				<< (double)stat.ribbon_total / stat.keys_total
			<< "\n";
	}
}
#endif

class PaCHash {
private:
	uint64_t n_buckets;
	EliasFano ef;
	mutable SimpleRibbon<1> ribbon;

	PaCHash(uint64_t n_buckets, EliasFano &&ef, SimpleRibbon<1> &&ribbon):
	  n_buckets(n_buckets), ef(std::move(ef)), ribbon(std::move(ribbon)) {}

public:
	PaCHash() : n_buckets(0) {
	}

	PaCHash(int k, double bucket_size, const std::vector<Hash128> &keys):
	  PaCHash(std::move(build(k, bucket_size, keys))) {}

	uint64_t operator()(Hash128 key) const {
		uint64_t bucket = bytehamster::util::fastrange64(key.hi, n_buckets);
		auto [start, end] = ef.search(bucket+1);
		int bits = std::bit_width(end - start);
		uint64_t x = 0;
		for (int j = 0; j < bits; j++) {
			x |= uint64_t(ribbon.retrieve(key.lo + uint64_t(j))) << j;
		}
		return start + x - 1;
	}

	size_t count_bits() const {
		return sizeof(*this) * 8
			+ ef.count_bits() - sizeof(ef) * 8
			+ ribbon.sizeBytes() * 8 - sizeof(ribbon) * 8
		;
	}

private:
	static PaCHash build(int k, double bucket_size,
	  std::vector<Hash128> keys) {
		uint64_t n_bins = (keys.size() + k - 1) / k;
		uint64_t n_buckets = std::ceil(keys.size() / bucket_size);

		for (Hash128 &key: keys) key.hi = bytehamster::util::fastrange64(key.hi, n_buckets);
		ips2ra::sort(keys.begin(), keys.end(),
			[](const Hash128 &key) -> uint64_t { return key.hi; });
		std::vector<uint64_t> bucket_sizes(n_buckets);
		for (uint64_t i = 0; i < keys.size(); i++) {
			bucket_sizes[keys[i].hi]++;
		}

		std::vector<uint64_t> threshold(n_bins+1);
		threshold.front() = 0;
		for (size_t i = 1; i < n_bins; i++) {
			uint64_t b = keys[k*i].hi;
			if (keys[k*i-1].hi != b && bucket_sizes[b-1] < bucket_sizes[b]) {
				b--;
			}
			threshold[i] = b + 1;
		}
		threshold.back() = n_buckets + 1;

		std::vector<std::pair<uint64_t, uint8_t>> ribbon_data;
		uint64_t offset = 0, next = 1;
		for (uint64_t i = 0; i < n_buckets; i++) {
			uint64_t start = next;
			while (threshold[next] == i+1) next++;
			uint64_t end = next;

			int bits = std::bit_width(end - start);
			uint64_t x = 0;
			for (Hash128 key; (key = keys[offset]).hi == i; offset++) {
				if (offset == (start+x) * k) x++;
				for (int j = 0; j < bits; j++) {
					ribbon_data.emplace_back(key.lo + uint64_t(j),
					                         (x >> j) & 1);
				}
			}
		}

#ifdef STATS
		{
			Stat &stat = stats[Config { k, bucket_size }];
			stat.ribbon_total += ribbon_data.size();
			stat.keys_total += keys.size();
		}
#endif

		return PaCHash(n_buckets,
		  std::move(EliasFano(threshold)),
		  std::move(SimpleRibbon<1>(ribbon_data)));
	}
};

}

