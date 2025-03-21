#pragma once

#include <vector>
#include <limits>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <ranges>

#include <sux/bits/EliasFano.hpp>
#include <ips2ra.hpp>
#include <bytehamster/util/Function.h>
#include <HashDisplace.hpp>
#include <OptimalBucketFunction.hpp>
#include <CompactEncoding.hpp>

#include <hash128.hpp>
#include "consensus.hpp"

#if 0
/* highly illegal hack */
namespace tlx {
	template<>
	inline unsigned clz<unsigned __int128>(unsigned __int128 i) {
		return std::countl_zero(i);
	}

	template<>
	inline unsigned ctz<unsigned __int128>(unsigned __int128 i) {
		return std::countr_zero(i);
	}
};
#endif

namespace kphf {
// TODO: Add non-consensus variant as well
namespace ThresholdBasedBumpingConsensus {

template<uint64_t n_thresholds>
constexpr std::array<uint64_t, n_thresholds>
  compute_thresholds(uint64_t _k, double bucket_size) {
	double lo = 0, hi = _k;

	for (int i = 0; i < 100; i++) {
		double mid = (lo+hi)/2;
		double prev = 0, cur = mid;
		for (uint64_t j = 0; j < n_thresholds-1; j++) {
			double x = cur / _k - prev / _k * std::pow(prev / cur, _k-1);
			prev = cur;
			if (x >= 1) {
				cur = std::numeric_limits<double>::infinity();
				break;
			} else {
				cur -= std::log1p(-x);
			}
		}
		if (cur > bucket_size
		  || cur / _k - prev / _k * std::pow(prev / cur, _k-1) > 1) {
			hi = mid;
		} else {
			lo = mid;
		}
	}

	auto convert = [](double x) -> uint64_t {
		x = std::round(std::ldexp(x, 64));
		if (x >= std::ldexp(1, 64)) {
			return std::numeric_limits<uint64_t>::max();
		} else {
			return uint64_t(x);
		}
	};

	std::array<uint64_t, n_thresholds> res;
	double prev = 0, cur = lo;
	for (uint64_t idx = 0; idx < n_thresholds; idx++) {
		res[idx] = convert(cur / bucket_size);
		double x = cur / _k - prev / _k * std::pow(prev / cur, _k-1);
		prev = cur;
		if (x == 1.0) cur = std::numeric_limits<double>::infinity();
		else cur -= std::log1p(-x);
	}

	return res;
}

constexpr double TUNE = 1.05;
auto compute_thresholds(uint64_t k, double lamda, int threshold_size) ->
  std::pair<std::vector<uint64_t>, std::vector<std::pair<uint64_t, uint64_t>>> {
	auto poisson = [](double lamda, uint64_t k) -> double {
		return exp(k * log(lamda) - lamda - lgamma(k+1));
	};

	auto success_prob =
	  [](double threshold, uint64_t n, uint64_t k) -> double {
		if (k > n) return 0;
		double lbinom = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
		double l1 = log(threshold);
		double l2 = log(1 - threshold);
		if (std::isinf(l1)) return k == 0;
		if (std::isinf(l2)) return k == n;
		return exp(lbinom + k * l1 + (n-k) * l2);
	};

	uint64_t n_thresholds = uint64_t(1) << threshold_size;
	std::vector<double> thresholds(n_thresholds);
	for (uint64_t i = 0; i < n_thresholds; i++) {
		thresholds[i] = (double) i / (n_thresholds-1);
	}

	uint64_t max_n = (uint64_t)(2 * k * lamda);
	std::vector<std::vector<double>> success_exp(max_n+1, std::vector<double>(k+1));

	auto recomp = [&]() {
		for (uint64_t n = 0; n <= max_n; n++) {
			auto &v = success_exp[n];
			for (uint64_t i = 0; i < v.size(); i++) {
				double e = 0;
				for (double t: thresholds) e += success_prob(t, n, i);
				v[i] = e;
			}
		}
	};

	auto set_threshold = [&](uint64_t i, double value) {
		double old = std::exchange(thresholds[i], value);
		for (uint64_t n = 0; n <= max_n; n++) {
			auto &v = success_exp[n];
			for (uint64_t i = 0; i < v.size(); i++) {
				v[i] += success_prob(value, n, i) - success_prob(old, n, i);
			}
		}
	};

	auto best_error = [&](uint64_t n) -> std::pair<double, double> {
		uint64_t goal = std::min(n, k);
		double need = TUNE;
		need -= success_exp[n][goal];
		if (need <= 0) return {0,0};
		double avg = 0;
		for (uint64_t error = 1; error <= goal; error++) {
			double s = success_exp[n][goal-error];
			if (s >= need) {
				return {error - 1 + need/s, avg + error * need};
			}
			need -= s;
			avg += error * s;
		}
		return {goal, std::numeric_limits<double>::infinity()};
	};

	auto eval = [&]() -> double {
		if (!std::ranges::contains(thresholds, 0.0)) {
			return std::numeric_limits<double>::infinity();
		}
		double exp = k * lamda;
		double badness = 0;
		for (uint64_t n = 0; n <= max_n; n++) {
			badness += poisson(exp, n) * best_error(n).second;
		}
		return badness;
	};

	double temp = 0.8;
	double cur_badness;
	while (temp > 0.01 / k) {
		recomp();
		cur_badness = eval();
		for (uint64_t j = 0; j < n_thresholds; j++) {
			double b = thresholds[j];
			double a = std::max(b - temp, 0.0);
			double c = std::min(b + temp, 1.0);
			set_threshold(j, a);
			double a_badness = eval();
			set_threshold(j, c);
			double c_badness = eval();
			double best = std::min({cur_badness, a_badness, c_badness});
			if (best == c_badness) ;
			else if (best == a_badness) set_threshold(j, a);
			else if (best == cur_badness) set_threshold(j, b);
			cur_badness = best;
		}
		temp *= 0.8;
	}
	std::ranges::sort(thresholds);
	recomp();

	auto double_to_u64 = [](double x) -> uint64_t {
		return x == 1.0 ? uint64_t(-1) : uint64_t(ldexp(x, 64));
	};

	std::vector<uint64_t> thresholds_final(n_thresholds);
	for (uint64_t i = 0; i < n_thresholds; i++) {
		thresholds_final[i] = double_to_u64(thresholds[i]);
	}

	std::vector<std::pair<uint64_t, uint64_t>> errors(max_n+1);
	for (uint64_t n = 0; n <= max_n; n++) {
		double e = best_error(n).first;
		uint64_t q = e;
		errors[n] = {q+1,double_to_u64(e-q)};
	}

	return {thresholds_final, errors};
}

#ifdef STATS
uint64_t perfect_thresholds = 0, total_thresholds = 0, extra_bumped = 0;
uint64_t bumped_keys = 0, total_keys = 0;
uint64_t filled_buckets = 0, overfull_buckets = 0;

void reset_stats() {
	perfect_thresholds = total_thresholds = extra_bumped = 0;
	bumped_keys = total_keys = 0;
	filled_buckets = overfull_buckets = 0;
}

void dump_stats() {
	std::cout << "Perfect thresholds: "
		<< 100.0 * perfect_thresholds / total_thresholds << "%\n";
	std::cout << "Avg. extra bumped: "
		<< (double) extra_bumped / total_thresholds << "\n";
	std::cout << "Bumped keys: " << 100.0 * bumped_keys / total_keys << "%\n";
	std::cout << "Filled buckets: " << 100.0 * filled_buckets / total_thresholds << "%\n";
	std::cout << "Overfull buckets: " << 100.0 * overfull_buckets / total_thresholds << "%\n";
}
#endif

struct Key {
	uint64_t bucket, fingerprint;
	Hash128 hash;
};

namespace {
	template<typename R>
	void sort_buckets(R &&r) {
		ips2ra::sort(r.begin(), r.end(),
		  [](const Key &key) -> uint64_t { return key.bucket; });
	}

	template<typename R>
	void sort_fingerprints(R &&r) {
		ips2ra::sort(r.begin(), r.end(),
		  [](const Key &key) -> uint64_t { return key.fingerprint; });
	}
}

struct SharedData {
	uint64_t k;
	double overload;
	int threshold_size;
	std::vector<uint64_t> thresholds;
	std::vector<std::pair<uint64_t, uint64_t>> errors;

	SharedData(uint64_t k, double overload, int threshold_size):
	  k(k), overload(overload), threshold_size(threshold_size) {
		assert(overload > 1.0);
		assert(threshold_size >= 1);
		assert(k > 0);
		std::tie(thresholds, errors) =
		  compute_thresholds(k, overload, threshold_size);
	}
};

class ThresholdBasedBumpingConsensus {
private:
	std::shared_ptr<SharedData> shared;
	uint64_t n;
	std::vector<std::pair<uint64_t, Consensus>> layers;
    using PHF = kphf::HashDisplace::HashDisplace<1, kphf::HashDisplace::OptimalBucketFunction, kphf::HashDisplace::CompactEncoding>;
	PHF phf;
	mutable sux::bits::EliasFano<> gaps;

public:

	uint64_t operator()(Hash128 key) const {
		uint64_t offset = 0;
		for (uint64_t i = 0; i < layers.size(); i++) {
			auto &[cur_buckets, consensus] = layers[i];
			uint64_t h = bytehamster::util::remix(key.hi + i);
			uint64_t b = bytehamster::util::fastrange64(h, cur_buckets);
			auto [seed,tidx] =
			  consensus.get(b * shared->threshold_size, shared->threshold_size);
			tidx = decrypt(seed, tidx, shared->threshold_size);
			uint64_t f = bytehamster::util::remix(key.lo + seed + i);

			if (f < shared->thresholds[tidx]) return offset + b;

			offset += cur_buckets;
		}

		return gaps.select(phf(key));
	}

	size_t count_bits() const {
		size_t s = 0;
		for (auto &[nb,c]: layers) s += c.count_bits() - sizeof(c) * 8;
		return sizeof(*this) * 8
			+ layers.capacity() * sizeof(layers[0])
		    + s
			+ phf.count_bits() - sizeof(phf) * 8
			+ gaps.bitCount() - sizeof(gaps) * 8
		;
	}

private:
	static inline uint64_t encrypt(uint64_t key, uint64_t tidx, uint64_t threshold_size) {
		return (tidx ^ bytehamster::util::remix(key)) & ((uint64_t(1) << threshold_size) - 1);
	}

	static inline uint64_t decrypt(uint64_t key, uint64_t tidx, uint64_t threshold_size) {
		return (tidx ^ bytehamster::util::remix(key)) & ((uint64_t(1) << threshold_size) - 1);
	}

	struct Builder {
		std::vector<Key> keys;
		std::ranges::subrange<std::vector<Key>::iterator> current;
		uint64_t cur_bucket, total_buckets, layer, offset;
		std::vector<uint64_t> *spots;
		std::vector<Hash128> *bumped;

		uint64_t k, threshold_size;
		std::span<uint64_t> thresholds;
		std::span<std::pair<uint64_t, uint64_t>> errors;

		Builder(SharedData *shared, std::vector<Hash128> &hashes,
		  uint64_t layer, uint64_t buckets, uint64_t offset,
		  std::vector<uint64_t> *spots):
		  keys(std::move(prepare(hashes, layer, buckets))),
		  current(keys.begin(), keys.begin()),
		  cur_bucket(0), total_buckets(buckets), layer(layer), offset(offset),
		  spots(spots), bumped(&hashes),
		  k(shared->k), threshold_size(shared->threshold_size),
		  thresholds(shared->thresholds), errors(shared->errors) {
			bumped->clear();
		}

		static std::vector<Key> prepare(const std::vector<Hash128> &hashes,
		  uint64_t seed, uint64_t n_buckets) {
			auto r = std::views::transform(hashes,
			  [seed, n_buckets](Hash128 key) {
				uint64_t b = bytehamster::util::fastrange64(bytehamster::util::remix(key.hi + seed), n_buckets);
				return Key { b, 0, key };
			});
			std::vector<Key> keys(r.begin(), r.end());
			sort_buckets(keys);
			return keys;
		}

		std::optional<uint64_t> advance(uint64_t seed) {
			if (cur_bucket == total_buckets) return {};

			auto it = current.end();
			while (it != keys.end() && it->bucket == cur_bucket) {
				it->fingerprint = bytehamster::util::remix(it->hash.lo + seed + layer);
				++it;
			}

			current = std::ranges::subrange(current.end(), it);
			sort_fingerprints(current);

			cur_bucket++;

			return threshold_size;
		}

		uint64_t backtrack() {
			cur_bucket--;

			if (cur_bucket == 0) {
				current = std::ranges::subrange(keys.begin(), keys.begin());
				return 0;
			}

			auto it = current.begin();
			while (it != keys.begin() && prev(it)->bucket == cur_bucket-1) {
				--it;
			}
			current = std::ranges::subrange(it, current.begin());

			uint64_t cnt = 0;
			while (!spots->empty() &&
			  spots->back() == offset + cur_bucket - 1) {
				spots->pop_back();
				cnt++;
			}
			bumped->resize(bumped->size() - (current.size() + cnt - k));

			return threshold_size;
		}

		std::optional<uint64_t> find_first(uint64_t seed) {
			if (current.size() > k) {
				uint64_t fp = current[k].fingerprint;
				auto it = std::ranges::upper_bound(thresholds, fp);
				return find(seed, it - begin(thresholds));
			} else {
				return find(seed, thresholds.size());
			}
		}

		std::optional<uint64_t> find_next(uint64_t seed, uint64_t prev) {
			return find(seed, decrypt(seed, prev, threshold_size));
		}

		std::optional<uint64_t> find(uint64_t seed, uint64_t prev) {
			uint64_t goal = std::min(k, current.size());
			uint64_t error_bound, error_prob;
			if (current.size() < errors.size()) {
				std::tie(error_bound, error_prob) =
				  errors[current.size()];
			} else {
				error_bound = k+1, error_prob = 0;
			}
			for (uint64_t tidx = prev; tidx > 0;) {
				tidx--;
				auto it = std::ranges::lower_bound(current, thresholds[tidx], {}, [](Key k) -> std::uint64_t { return k.fingerprint; });
				uint64_t idx = it - current.begin();
				uint64_t error = goal - idx;
				if (error > error_bound) break;
				if (error < error_bound || bytehamster::util::remix(seed + tidx) < error_prob) {
					for (Key &k: current | std::views::drop(idx)) {
						bumped->push_back(k.hash);
					}
					for (uint64_t x = idx; x < k; x++) {
						spots->push_back(offset + cur_bucket - 1);
					}
					return encrypt(seed, tidx, threshold_size);
				}
			}
			return {};
		}
	};

	public:
	ThresholdBasedBumpingConsensus()
		: n(0), shared(new SharedData(8, 2.0, 3)), gaps({}, 0) {
	}

    ThresholdBasedBumpingConsensus(size_t k, std::vector<Hash128> &keys, double overload, size_t thresholdSize)
            : n(keys.size()), shared(new SharedData(k, overload, thresholdSize)), gaps({}, 0) {
		double overload_bucket_size = k * shared->overload;
		uint64_t total_buckets = (n + k - 1) / k;

#ifdef STATS
		total_keys += n;
#endif

		uint64_t offset = 0;
		std::vector<uint64_t> spots;
		for (uint64_t i = 0; offset != total_buckets; i++) {
			uint64_t remaining = total_buckets - offset;
			uint64_t cur_buckets =
			  std::ceil(keys.size() / overload_bucket_size);
			if (cur_buckets >= remaining) cur_buckets = remaining;
			else if ((double)(keys.size() - k * cur_buckets) / (remaining - cur_buckets) > overload_bucket_size) cur_buckets = remaining;

			Builder builder(shared.get(), keys, i, cur_buckets, offset, &spots);
			Consensus consensus(builder);

			layers.emplace_back(cur_buckets, std::move(consensus));

			offset += cur_buckets;
		}
		layers.shrink_to_fit();

#ifdef STATS
		bumped_keys += keys.size();
#endif

		phf = PHF(keys, 5);
		std::vector<uint64_t> actual_spots;
		for (Hash128 k: keys) {
			uint64_t h = phf(k);
			if (h >= actual_spots.size()) actual_spots.resize(h+1);
			actual_spots[h] = 1;
		}
		auto it = spots.begin();
		for (uint64_t &x: actual_spots) {
			uint64_t v = *it;
			if (x) ++it;
			x = v;
		}
		gaps = sux::bits::EliasFano<>(actual_spots, total_buckets);
	}
};

}
}
