#pragma once

#include <bit>
#include <vector>
#include <ranges>
#ifdef STATS
#include <map>
#endif

#include <ips2ra.hpp>
#include <bytehamster/util/Function.h>
#include <sux/bits/SimpleSelectHalf.hpp>

#include <hash128.hpp>

namespace kphf::HashDisplace {

#ifdef STATS
struct Config {
	uint64_t n, k, bucket_size;

	bool operator<(const Config &rhs) const {
		return n != rhs.n ? n < rhs.n :
		       k != rhs.k ? k < rhs.k : bucket_size < rhs.bucket_size;
	}
};

template<typename BucketFunction>
std::map<Config, std::pair<uint64_t, std::vector<uint64_t>>> stats;

template<typename BucketFunction>
void dump_stats() {
	for (const auto &[config, data]: stats<BucketFunction>) {
		const auto &[count, totals] = data;
		for (size_t i = 0; i < totals.size(); i++) {
			std::cout << "RESULT"
			  << "\tbucket_function=" << BucketFunction::name()
			  << "\tn=" << config.n
			  << "\tk=" << config.k
			  << "\tbucket_size=" << config.bucket_size
			  << "\tbucket_idx=" << i
			  << "\tavg_seed=" << (double) totals[i]/count
			  << "\n";
		}
	}
}
#endif

template<typename BucketFunction, typename Encoding>
class HashDisplaceBuilder;

template<typename BucketFunction, typename Encoding>
class HashDisplace {
private:
	uint64_t nbins;
	BucketFunction::Instance bucket_function;
	Encoding seeds;

	HashDisplace(uint64_t nbins, BucketFunction::Instance &&bucket_function,
	  Encoding &&seeds): nbins(nbins), bucket_function(std::move(bucket_function)),
	  seeds(std::move(seeds)) {}

public:
	uint64_t operator()(Hash128 item) const {
		uint64_t bucket = bucket_function(item.hi);
		uint64_t seed = seeds[bucket];
		return bytehamster::util::fastrange64(bytehamster::util::remix(item.lo + seed), nbins);
	}

	size_t count_bits() const {
		return 8 * sizeof(*this)
		  + (bucket_function.count_bits() - 8 * sizeof(bucket_function))
		  + (seeds.count_bits() - 8 * sizeof(seeds));
	}

	friend class HashDisplaceBuilder<BucketFunction, Encoding>;
};

template<typename BucketFunction, typename Encoding>
class HashDisplaceBuilder {
private:
	uint64_t k;
	BucketFunction bucket_function;

public:
	HashDisplaceBuilder(uint64_t k): k(k), bucket_function(k) {}

	HashDisplace<BucketFunction, Encoding> operator()(
		const std::vector<Hash128> &items,
		uint64_t bucket_size,
		double load_factor = 1.0
	) const {
		assert(bucket_size > 0);
		assert(0.0 <= load_factor && load_factor <= 1.0);

		if (items.empty()) {
			return HashDisplace<BucketFunction, Encoding>(0,
			  bucket_function(0, 0, 1.0),
			  Encoding(std::vector<uint64_t> {}));
		}

		uint64_t nbuckets = (items.size() + bucket_size - 1) / bucket_size;
		uint64_t nbins = ceil(items.size() / load_factor / k);
		typename BucketFunction::Instance bf_instance =
		  bucket_function(items.size(), nbuckets, load_factor);

		auto r = std::views::transform(items,
		  [&bf_instance](Hash128 item) -> std::pair<uint64_t, uint64_t> {
			return {bf_instance(item.hi), item.lo};
		});
		std::vector<std::pair<uint64_t, uint64_t>>
		  sorted_items(r.begin(), r.end());
		ips2ra::sort(sorted_items.begin(), sorted_items.end(),
		  [](std::pair<uint64_t, uint64_t> k) -> uint64_t { return k.first; });

		using Iter = decltype(sorted_items)::iterator;
		using Range = std::ranges::subrange<Iter>;
		struct Bucket {
			uint64_t idx;
			Range range;
		};
		std::vector<Bucket> buckets(nbuckets);
		{
			Iter it = sorted_items.begin();
			for (uint64_t i = 0; i < nbuckets; i++) {
				Iter start = it;
				while (it != sorted_items.end() && it->first == i) ++it;
				buckets[i] = {i, Range(start, it)};
			}
		}
		ips2ra::sort(buckets.rbegin(), buckets.rend(),
		  [](const Bucket &b) -> size_t {
			return b.range.size();
		});
		std::vector<uint32_t> counts(nbins, 0);
		std::vector<uint64_t> seeds_vec(nbuckets);
#ifdef STATS
		auto &[stats_count, stats_vec] =
		  stats<BucketFunction>.try_emplace(Config { items.size(), k, bucket_size }, std::piecewise_construct, std::make_tuple(0), std::make_tuple(nbuckets)).first->second;
		stats_count++;
		auto stats_it = stats_vec.begin();
#endif
		for (auto [i, bucket]: buckets) {
			uint64_t seed;
			for (seed = 0;; seed++) {
				uint64_t j;
				for (j = 0; j < bucket.size(); j++) {
					uint64_t hash = bytehamster::util::fastrange64(bytehamster::util::remix(bucket[j].second + seed), nbins);
					if (counts[hash] == k) break;
					counts[hash]++;
				}
				if (j == bucket.size()) break;
				while (j > 0) {
					uint64_t hash = bytehamster::util::fastrange64(bytehamster::util::remix(bucket[--j].second + seed), nbins);
					counts[hash]--;
				}
			}
			seeds_vec[i] = seed;
#ifdef STATS
			*stats_it++ += seed;
#endif
		}
		Encoding seeds(seeds_vec);

		return HashDisplace<BucketFunction, Encoding>(nbins, std::move(bf_instance), std::move(seeds));
	}
};

}
