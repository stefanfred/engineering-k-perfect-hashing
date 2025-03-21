/*
 * Sux: Succinct data structures
 *
 * Copyright (C) 2019-2020 Emmanuel Esposito and Sebastiano Vigna
 *
 *  This library is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * Under Section 7 of GPL version 3, you are granted additional permissions
 * described in the GCC Runtime Library Exception, version 3.1, as published by
 * the Free Software Foundation.
 *
 * You should have received a copy of the GNU General Public License and a copy of
 * the GCC Runtime Library Exception along with this program; see the files
 * COPYING3 and COPYING.RUNTIME respectively.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <sux/support/SpookyV2.hpp>
#include <sux/util/Vector.hpp>
#include <sux/bits/Rank9.hpp>
#include <sux/bits/SimpleSelectHalf.hpp>
#include <sux/function/DoubleEF.hpp>
#include <sux/function/RiceBitVector.hpp>

namespace sux::function::krecsplit {

using namespace std;
using namespace std::chrono;
using sux::function::RiceBitVector;
using sux::function::DoubleEF;
using sux::bits::Rank9;
using sux::bits::SimpleSelectHalf;
using sux::util::Vector;

template<size_t K, size_t LEAF_SIZE>
static constexpr bool is_small = (K <= 20 && LEAF_SIZE <= 2) || K < 10;

// Assumed *maximum* size of a bucket.
template<size_t K, size_t LEAF_SIZE>
static constexpr size_t MAX_BUCKET_SIZE = is_small<K, LEAF_SIZE> ? 3000 : 30000;

static const int MAX_FANOUT = 32;

#if defined(MORESTATS) && !defined(STATS)
#define STATS
#endif

#ifdef MORESTATS

#define MAX_LEVEL_TIME (20)

static constexpr double log2e = 1.44269504089;
static uint64_t num_split_trials;
static uint64_t num_split_evals;
static uint64_t split_count;
static uint64_t expected_split_trials, expected_split_evals;
static uint64_t split_unary, split_fixed, split_unary_golomb, split_fixed_golomb;
static uint64_t max_split_code, min_split_code, sum_split_codes;
static uint64_t sum_depths;
static uint64_t time_split[MAX_LEVEL_TIME];
#endif

static const uint64_t start_seed[] = {0x106393c187cae21a, 0x6453cec3f7376937, 0x643e521ddbd2be98, 0x3740c6412f6572cb, 0x717d47562f1ce470, 0x4cd6eb4c63befb7c, 0x9bfd8c5e18c8da73,
									  0x082f20e10092a9a3, 0x2ada2ce68d21defc, 0xe33cb4f3e7c6466b, 0x3980be458c509c59, 0xc466fd9584828e8c, 0x45f0aabe1a61ede6, 0xf6e7b8b33ad9b98d,
									  0x4ef95e25f4b4983d, 0x81175195173b92d3, 0x4e50927d8dd15978, 0x1ea2099d1fafae7f, 0x425c8a06fbaaa815, 0xcd4216006c74052a};

/** David Stafford's (http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html)
 * 13th variant of the 64-bit finalizer function in Austin Appleby's
 * MurmurHash3 (https://github.com/aappleby/smhasher).
 *
 * @param z a 64-bit integer.
 * @return a 64-bit integer obtained by mixing the bits of `z`.
 */
// TODO: Remove functions that are duplicates of things that can be found in RecSplit/Sux,
//  no need to define them again
uint64_t inline remix(uint64_t z) {
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}

/** 128-bit hashes.
 *
 * In the construction of RecSplit, keys are replaced with instances
 * of this class using SpookyHash, first thing.
 * Moreover, it is possible to build and query RecSplit instances using 128-bit
 * random hashes only (mainly for benchmarking purposes).
 */

typedef struct __hash128_t {
	uint64_t first, second;
	bool operator<(const __hash128_t &o) const { return first < o.first || second < o.second; }
	__hash128_t(const uint64_t first, const uint64_t second) {
		this->first = first;
		this->second = second;
	}
} hash128_t;

/** Convenience function hashing a key a returning a __hash128_t
 *
 * @param data a pointer to the key.
 * @param length the length in bytes of the key.
 * @param seed an additional seed.
 */

hash128_t inline spooky(const void *data, const size_t length, const uint64_t seed) {
	uint64_t h0 = seed, h1 = seed;
	SpookyHash::Hash128(data, length, &h0, &h1);
	return {h1, h0};
}

// Quick replacements for min/max on not-so-large integers.

static constexpr inline uint64_t min(int64_t x, int64_t y) { return y + ((x - y) & ((x - y) >> 63)); }
static constexpr inline uint64_t max(int64_t x, int64_t y) { return x - ((x - y) & ((x - y) >> 63)); }

/** A class emboding the splitting strategy of k-RecSplit.
 *
 *  Note that this class is used _for statistics only_. The splitting strategy is embedded
 *  into the generation code, which uses only the public fields SplittingStrategy::lower_aggr and SplittingStrategy::upper_aggr.
 */

template <size_t K, size_t LEAF_SIZE> class SplittingStrategy {
	static constexpr size_t _k = K;
	static constexpr size_t _leaf_size = LEAF_SIZE;
	static constexpr size_t _leaf_keys = _k * _leaf_size;
	static_assert(_k >= 1);
	static_assert(_leaf_size >= 1);
	size_t m, curr_unit, curr_index, last_unit;
	size_t _fanout;
	size_t unit;

	inline size_t part_size() const { return (curr_index < _fanout - 1) ? unit : last_unit; }

	static inline constexpr size_t split_difficulty(size_t part_count, size_t part_size) {
		double d = 1;
		for (size_t j = 0; j < part_size; j++) {
			for (size_t i = 0; i < part_count; i++) {
				d *= (double) part_count * (j+1) / (j*part_count+i+1);
			}
		}
		return d;
	}

	static constexpr size_t difficulty = split_difficulty(_leaf_size, _k);

	static inline constexpr size_t next_aggr(size_t cur) {
		double prev = split_difficulty(2, cur);
		for (size_t i = 3;; i++) {
			double d = split_difficulty(i, cur);
			if (d > difficulty) {
				if (difficulty < (prev + d)/2) return cur * (i-1);
				else return cur * i;
			}
			prev = d;
		}
	}

  public:
	/** The lower bound for primary (lower) key aggregation. */
	static constexpr size_t lower_aggr = next_aggr(_leaf_keys);
	/** The lower bound for secondary (upper) key aggregation. */
	static constexpr size_t upper_aggr = next_aggr(lower_aggr);

	static inline constexpr void split_params(const size_t m, size_t &fanout, size_t &unit) {
		if (m > upper_aggr) { // High-level aggregation (fanout 2)
			unit = upper_aggr * (uint16_t(m / 2 + upper_aggr - 1) / upper_aggr);
			fanout = 2;
		} else if (m > lower_aggr) { // Second-level aggregation
			unit = lower_aggr;
			fanout = uint16_t(m + lower_aggr - 1) / lower_aggr;
		} else if (m > _leaf_keys) { // First-level aggregation
			unit = _leaf_keys;
			fanout = uint16_t(m + _leaf_keys - 1) / _leaf_keys;
		} else { // Leaf
			unit = _k;
			fanout = uint16_t(m + _k - 1) / _k;
		}
	}

	// Note that you can call this iterator only *once*.
	class split_iterator {
		SplittingStrategy *strat;

	  public:
		using value_type = size_t;
		using difference_type = ptrdiff_t;
		using pointer = size_t *;
		using reference = size_t &;
		using iterator_category = input_iterator_tag;

		split_iterator(SplittingStrategy *strat) : strat(strat) {}
		size_t operator*() const { return strat->curr_unit; }
		size_t *operator->() const { return &strat->curr_unit; }
		split_iterator &operator++() {
			++strat->curr_index;
			strat->curr_unit = strat->part_size();
			strat->last_unit -= strat->curr_unit;
			return *this;
		}
		bool operator==(const split_iterator &other) const { return strat == other.strat; }
		bool operator!=(const split_iterator &other) const { return !(*this == other); }
	};

	explicit SplittingStrategy(size_t m) : m(m), last_unit(m), curr_index(0), curr_unit(0) {
		split_params(m, _fanout, unit);
		this->curr_unit = part_size();
		this->last_unit -= this->curr_unit;
	}

	split_iterator begin() { return split_iterator(this); }
	split_iterator end() { return split_iterator(nullptr); }

	inline size_t fanout() { return this->_fanout; }
};

// Generates the precomputed table of 32-bit values holding the Golomb-Rice code
// of a splitting (upper 5 bits), the number of nodes in the associated subtree
// (following 11 bits) and the sum of the Golomb-Rice codelengths in the same
// subtree (lower 16 bits).

template <size_t K, size_t LEAF_SIZE> static constexpr void _fill_golomb_rice(const int m, array<uint32_t, MAX_BUCKET_SIZE<K, LEAF_SIZE>> *memo) {
	array<int, MAX_FANOUT> k{0};

	size_t fanout = 0, unit = 0;
	SplittingStrategy<K, LEAF_SIZE>::split_params(m, fanout, unit);

	k[fanout - 1] = m;
	for (size_t i = 0; i < fanout - 1; ++i) {
		k[i] = unit;
		k[fanout - 1] -= k[i];
	}

	double sqrt_prod = 1;
	for (size_t i = 0; i < fanout; ++i) {
		if (k[i] > 1) sqrt_prod *= sqrt(k[i]);
		else sqrt_prod *= numbers::e / sqrt(2 * numbers::pi);
	}

	const double p = sqrt(m) / (pow(2 * M_PI, (fanout - 1.) / 2) * sqrt_prod);
	auto golomb_rice_length = (uint32_t)ceil(log2(-log((sqrt(5) + 1) / 2) / log1p(-p))); // log2 Golomb modulus

	assert(golomb_rice_length <= 0x1F); // Golomb-Rice code, stored in the 5 upper bits
	(*memo)[m] = golomb_rice_length << 27;
	for (size_t i = 0; i < fanout; ++i) golomb_rice_length += (*memo)[k[i]] & 0xFFFF;
	assert(golomb_rice_length <= 0xFFFF); // Sum of Golomb-Rice codeslengths in the subtree, stored in the lower 16 bits
	(*memo)[m] |= golomb_rice_length;

	uint32_t nodes = 1;
	for (size_t i = 0; i < fanout; ++i) nodes += ((*memo)[k[i]] >> 16) & 0x7FF;
	assert(nodes <= 0x7FF); // Number of nodes in the subtree, stored in the middle 11 bits
	(*memo)[m] |= nodes << 16;
}

template <size_t K, size_t LEAF_SIZE> static constexpr array<uint32_t, MAX_BUCKET_SIZE<K, LEAF_SIZE>> fill_golomb_rice() {
	array<uint32_t, MAX_BUCKET_SIZE<K, LEAF_SIZE>> memo{0};
	for (size_t s = K+1; s < MAX_BUCKET_SIZE<K, LEAF_SIZE>; ++s) _fill_golomb_rice<K, LEAF_SIZE>(s, &memo);
	return memo;
}

// Computes the Golomb modulu of a splitting (for statistics purposes only)

template <size_t K, size_t LEAF_SIZE> static constexpr uint64_t split_golomb_b(const int m) {
	array<int, MAX_FANOUT> k{0};

	size_t fanout = 0, unit = 0;
	SplittingStrategy<K, LEAF_SIZE>::split_params(m, fanout, unit);

	k[fanout - 1] = m;
	for (size_t i = 0; i < fanout - 1; ++i) {
		k[i] = unit;
		k[fanout - 1] -= k[i];
	}

	double sqrt_prod = 1;
	for (size_t i = 0; i < fanout; ++i) sqrt_prod *= sqrt(k[i]);

	const double p = sqrt(m) / (pow(2 * M_PI, (fanout - 1.) / 2) * sqrt_prod);
	return ceil(-log(2 - p) / log1p(-p)); // Golomb modulus
}

#define first_hash(k, len) spooky(k, len, 0)
#define golomb_param(m) (memo[m] >> 27)
#define skip_bits(m) (memo[m] & 0xFFFF)
#define skip_nodes(m) ((memo[m] >> 16) & 0x7FF)

enum class BumpStrategy {
	RANK_RECURSE,
	RANK_SPLIT,
	INTERSPERSE_SPLIT,
	INTERSPERSE_SPLIT_NOBITS,
	BINS_MOD_RECURSE,
	SELECT_RECURSE,
};

template <size_t K, size_t LEAF_SIZE, BumpStrategy BS, util::AllocType AT = util::AllocType::MALLOC> class RecSplit;

template<size_t K, size_t LEAF_SIZE, BumpStrategy BS, util::AllocType AT>
struct ExtraFields {
	size_t extraBitCount() { return 0; }
#ifdef STATS
	void stats(auto &show) {
		(void) show;
	}
#endif
};

template<size_t K, size_t LEAF_SIZE, util::AllocType AT>
struct ExtraFields<K, LEAF_SIZE, BumpStrategy::RANK_RECURSE, AT> {
	size_t unbumped_bins;
	unique_ptr<RecSplit<K, LEAF_SIZE, BumpStrategy::RANK_RECURSE, AT>> recurse_bumped = nullptr;
	Rank9<AT> *buckets_bumped_rank = nullptr;
	Vector<uint64_t, AT> buckets_bumped;

	size_t extraBitCount() {
		return
			(recurse_bumped == nullptr ? 0 : recurse_bumped->bitCount())
			+ buckets_bumped_rank->bitCount() - 8 * sizeof(*buckets_bumped_rank)
			+ buckets_bumped.bitCount() - 8 * sizeof(buckets_bumped);
	}

#ifdef STATS
	void stats(auto &show) {
		show("Recursive structure",
		  recurse_bumped == nullptr ? 0 : recurse_bumped->bitCount());
		show("Bumped buckets",
		  buckets_bumped.bitCount() - 8 * sizeof(buckets_bumped));
		show("Rank on bumped buckets",
		  buckets_bumped_rank.bitCount() - 8 * sizeof(buckets_bumped_rank));
	}
#endif
};

template<size_t K, size_t LEAF_SIZE, util::AllocType AT>
struct ExtraFields<K, LEAF_SIZE, BumpStrategy::RANK_SPLIT, AT> {
	size_t unbumped_bins;
	Rank9<AT> buckets_bumped_rank;
	Vector<uint64_t, AT> buckets_bumped;

	size_t extraBitCount() {
		return
			buckets_bumped_rank.bitCount() - 8 * sizeof(buckets_bumped_rank)
			+ buckets_bumped.bitCount() - 8 * sizeof(buckets_bumped);
	}

#ifdef STATS
	void stats(auto &show) {
		show("Bumped buckets",
		  buckets_bumped.bitCount() - 8 * sizeof(buckets_bumped));
		show("Rank on bumped buckets",
		  buckets_bumped_rank.bitCount() - 8 * sizeof(buckets_bumped_rank));
	}
#endif
};

template<size_t K, size_t LEAF_SIZE, util::AllocType AT>
struct ExtraFields<K, LEAF_SIZE, BumpStrategy::INTERSPERSE_SPLIT, AT> {
	Vector<uint64_t, AT> buckets_bumped;

	size_t extraBitCount() {
		return buckets_bumped.bitCount() - 8 * sizeof(buckets_bumped);
	}

#ifdef STATS
	void stats(auto &show) {
		show("Bumped buckets",
		  buckets_bumped.bitCount() - 8 * sizeof(buckets_bumped));
	}
#endif
};

template<size_t K, size_t LEAF_SIZE, util::AllocType AT>
struct ExtraFields<K, LEAF_SIZE, BumpStrategy::BINS_MOD_RECURSE, AT> {
	size_t unbumped_bins;
	unique_ptr<RecSplit<K, LEAF_SIZE, BumpStrategy::BINS_MOD_RECURSE, AT>>
	  recurse_bumped;
	Vector<uint8_t, AT> bucket_mods;

	size_t extraBitCount() {
		return
			(recurse_bumped == nullptr ? 0 : recurse_bumped->bitCount())
			+ bucket_mods.bitCount() - 8 * sizeof(bucket_mods);
	}

#ifdef STATS
	void stats(auto &show) {
		show("Recursive structure",
		  recurse_bumped == nullptr ? 0 : recurse_bumped->bitCount());
		show("Bucket size modulo K",
		  bucket_mods.bitCount() - 8 * sizeof(bucket_mods));
	}
#endif
};

template<size_t K, size_t LEAF_SIZE, util::AllocType AT>
struct ExtraFields<K, LEAF_SIZE, BumpStrategy::SELECT_RECURSE, AT> {
	unique_ptr<RecSplit<K, LEAF_SIZE, BumpStrategy::SELECT_RECURSE, AT>>
	  recurse_bumped;
	SimpleSelectHalf<AT> buckets_bumped_select;
	Vector<uint64_t, AT> buckets_bumped;

	size_t extraBitCount() {
		return
			(recurse_bumped == nullptr ? 0 : recurse_bumped->bitCount())
			+ buckets_bumped_select.bitCount()
			  - 8 * sizeof(buckets_bumped_select)
			+ buckets_bumped.bitCount() - 8 * sizeof(buckets_bumped);
	}

#ifdef STATS
	void stats(auto &show) {
		show("Recursive structure",
		  recurse_bumped == nullptr ? 0 : recurse_bumped->bitCount());
		show("Bumped buckets",
		  buckets_bumped.bitCount() - 8 * sizeof(buckets_bumped));
		show("Select on bumped buckets",
		  buckets_bumped_select.bitCount() - 8 * sizeof(buckets_bumped_select));
	}
#endif
};

/**
 *
 * A class for storing minimal perfect hash functions. The template
 * parameter decides how large a leaf will be. Larger leaves imply
 * slower construction, but less space and faster evaluation.
 *
 * @tparam LEAF_SIZE the size of a leaf; typicals value range from 6 to 8
 * for fast, small maps, or up to 16 for very compact functions.
 * @tparam AT a type of memory allocation out of sux::util::AllocType.
 */

template <size_t K, size_t LEAF_SIZE, BumpStrategy BS, util::AllocType AT> class RecSplit: ExtraFields<K, LEAF_SIZE, BS, AT> {
	using SplitStrat = SplittingStrategy<K, LEAF_SIZE>;

	static constexpr size_t _k = K;
	static constexpr size_t _leaf_size = LEAF_SIZE;
	static constexpr size_t leaf_keys = _k * _leaf_size;
	static constexpr size_t lower_aggr = SplitStrat::lower_aggr;
	static constexpr size_t upper_aggr = SplitStrat::upper_aggr;

	// For each bucket size, the Golomb-Rice parameter (upper 8 bits) and the number of bits to
	// skip in the fixed part of the tree (lower 24 bits).
	static constexpr array<uint32_t, MAX_BUCKET_SIZE<K, LEAF_SIZE>> memo = fill_golomb_rice<K, LEAF_SIZE>();

	size_t bucket_size;
	size_t nbuckets;
	size_t keys_count;
	RiceBitVector<AT> descriptors;
	DoubleEF<AT> ef;

  public:
	RecSplit() {}

	/** Builds a RecSplit instance using a given list of keys and bucket size.
	 *
	 * **Warning**: duplicate keys will cause this method to never return.
	 *
	 * @param keys a vector of strings.
	 * @param bucket_size the desired bucket size; typical sizes go from
	 * 100 to 2000, with smaller buckets giving slightly larger but faster
	 * functions.
	 */
	RecSplit(const vector<string> &keys, const size_t bucket_size) {
		this->bucket_size = bucket_size;
		this->keys_count = keys.size();
		hash128_t *h = (hash128_t *)malloc(this->keys_count * sizeof(hash128_t));
		for (size_t i = 0; i < this->keys_count; ++i) {
			h[i] = first_hash(keys[i].c_str(), keys[i].size());
		}
		hash_gen(h);
		free(h);
	}

	/** Builds a RecSplit instance using a given list of 128-bit hashes and bucket size.
	 *
	 * **Warning**: duplicate keys will cause this method to never return.
	 *
	 * Note that this constructor is mainly useful for benchmarking.
	 * @param keys a vector of 128-bit hashes.
	 * @param bucket_size the desired bucket size; typical sizes go from
	 * 100 to 2000, with smaller buckets giving slightly larger but faster
	 * functions.
	 */
	RecSplit(vector<hash128_t> &keys, const size_t bucket_size) {
		this->bucket_size = bucket_size;
		this->keys_count = keys.size();
		hash_gen(&keys[0]);
	}

	/** Builds a RecSplit instance using a list of keys returned by a stream and bucket size.
	 *
	 * **Warning**: duplicate keys will cause this method to never return.
	 *
	 * @param input an open input stream returning a list of keys, one per line.
	 * @param bucket_size the desired bucket size.
	 */
	RecSplit(ifstream &input, const size_t bucket_size) {
		this->bucket_size = bucket_size;
		vector<hash128_t> h;
		for (string key; getline(input, key);) h.push_back(first_hash(key.c_str(), key.size()));
		this->keys_count = h.size();
		hash_gen(&h[0]);
	}

	/** Returns the value associated with the given 128-bit hash.
	 *
	 * Note that this method is mainly useful for benchmarking.
	 * @param hash a 128-bit hash.
	 * @return the associated value.
	 */
	size_t operator()(const hash128_t &hash) {
		size_t bucket = hash128_to_bucket(hash);
		uint64_t cum_bins, cum_bumped, cum_keys, bit_pos;
		size_t m;
		if constexpr (BS == BumpStrategy::BINS_MOD_RECURSE) {
			uint64_t cum_bins_next;
			ef.get(bucket, cum_bins, cum_bins_next, bit_pos);
			size_t bit = bucket * bit_width(_k - 1);
			uint64_t mod;
			memcpy(&mod, &this->bucket_mods + bit / 8, 8);
			mod >>= (bit&7);
			mod &= (1 << bit_width(_k - 1)) - 1;
			m = (cum_bins_next - cum_bins) * _k + mod;
		} else {
			uint64_t cum_keys_next;
			ef.get(bucket, cum_keys, cum_keys_next, bit_pos);
			if constexpr (BS == BumpStrategy::RANK_RECURSE
			  || BS == BumpStrategy::RANK_SPLIT) {
				cum_bumped = this->buckets_bumped_rank->rank(bucket);
			} else {
				cum_bumped = 0;
			}
			cum_bins = cum_keys / _k - cum_bumped;
			m = cum_keys_next - cum_keys;
		}

		auto reader = descriptors.reader();
		reader.readReset(bit_pos, skip_bits(m));
		int level = 0;

		while (m > upper_aggr) { // fanout = 2
			const auto d = reader.readNext(golomb_param(m));
			const size_t hmod = remap16(remix(hash.second + d + start_seed[level]), m);

			const uint32_t split_bins = ((uint16_t(m / 2 + upper_aggr - 1) / upper_aggr)) * (upper_aggr / _k);
			const uint32_t split = split_bins * _k;
			if (hmod < split) {
				m = split;
			} else {
				reader.skipSubtree(skip_nodes(split), skip_bits(split));
				m -= split;
				cum_bins += split_bins;
			}
			level++;
		}
		if (m > lower_aggr) {
			const auto d = reader.readNext(golomb_param(m));
			const size_t hmod = remap16(remix(hash.second + d + start_seed[level]), m);

			const int part = uint16_t(hmod) / lower_aggr;
			m = min(lower_aggr, m - part * lower_aggr);
			cum_bins += (lower_aggr/_k) * part;
			if (part) reader.skipSubtree(skip_nodes(lower_aggr) * part, skip_bits(lower_aggr) * part);
			level++;
		}

		if (m > leaf_keys) {
			const auto d = reader.readNext(golomb_param(m));
			const size_t hmod = remap16(remix(hash.second + d + start_seed[level]), m);

			const int part = uint16_t(hmod) / leaf_keys;
			m = min(leaf_keys, m - part * leaf_keys);
			cum_bins += _leaf_size * part;
			if (part) reader.skipSubtree(part, skip_bits(leaf_keys) * part);
			level++;
		}

		if (m > _k) {
			const auto d = reader.readNext(golomb_param(m));
			const size_t hmod = remap16(remix(hash.second + d + start_seed[level]), m);

			const int part = uint16_t(hmod) / _k;
			m = min(_k, m - part * _k);
			cum_bins += part;
			level++;
		}

		if (m == _k) return cum_bins;

		if constexpr (BS == BumpStrategy::RANK_RECURSE
		  || BS == BumpStrategy::BINS_MOD_RECURSE
		  || BS == BumpStrategy::SELECT_RECURSE) {
			size_t rec;
			if (this->recurse_bumped == nullptr) rec = 0;
			else rec = (*this->recurse_bumped)(hash);
			if constexpr (BS == BumpStrategy::SELECT_RECURSE) {
				bucket = this->buckets_bumped_select.select(rec);
				if (bucket == nbuckets) return keys_count / _k;
				uint64_t cum_keys_next;
				ef.get(bucket, cum_keys, cum_keys_next, bit_pos);
				return cum_keys_next / _k - 1;
			} else {
				return this->unbumped_bins + rec;
			}
		} else if constexpr (BS == BumpStrategy::RANK_SPLIT
		  || BS == BumpStrategy::INTERSPERSE_SPLIT
		  || BS == BumpStrategy::INTERSPERSE_SPLIT_NOBITS) {
			const size_t split = _k - cum_keys % _k;

			bool split_right;
			if (split >= m) split_right = false;
			else {
				const double p = sqrt(m / (2 * M_PI * split * (m - split)));
				const uint64_t seed = end(start_seed)[-1] + reader.readGolomb(p);
				const size_t hmod = remap16(remix(hash.second + seed), m);
				split_right = hmod >= split;
			}

			if constexpr (BS == BumpStrategy::RANK_SPLIT) {
				return this->unbumped_bins + cum_bumped + split_right;
			} else {
				bucket += split_right;
				uint64_t cum_keys_next;
				if constexpr (BS == BumpStrategy::INTERSPERSE_SPLIT) {
					uint64_t v = this->buckets_bumped[bucket/64] >> (bucket&63);
					if (v == 0) {
						bucket &= ~uint64_t(63);
						do {
							bucket += 64;
							v = this->buckets_bumped[bucket/64];
						} while (v == 0);
					}
					bucket += rho(v);
					if (bucket == nbuckets) return keys_count / _k;
					ef.get(bucket, cum_keys, cum_keys_next, bit_pos);
				} else if constexpr (BS == BumpStrategy::INTERSPERSE_SPLIT_NOBITS) {
					do {
						if (bucket == nbuckets) return keys_count / _k;
						ef.get(bucket++, cum_keys, cum_keys_next, bit_pos);
					} while (cum_keys % _k <= cum_keys_next % _k);
				}
				return cum_keys_next / _k - 1;
			}

		}
	}

	/** Returns the value associated with the given key.
	 *
	 * @param key a key.
	 * @return the associated value.
	 */
	size_t operator()(const string &key) { return operator()(first_hash(key.c_str(), key.size())); }

	/** Returns the number of keys used to build this RecSplit instance. */
	inline size_t size() { return this->keys_count; }

	size_t bitCount() {
		return 8 * sizeof(*this)
			+ descriptors.getBits()
			+ ef.bitCountCumKeys() + ef.bitCountPosition()
			+ this->extraBitCount();
		;
	}

  private:
	// Maps a 128-bit to a bucket using the first 64-bit half.
	inline uint64_t hash128_to_bucket(const hash128_t &hash) const { return remap128(hash.first, nbuckets); }

	// Computes and stores the splittings and bijections of a bucket.
	void recSplit(vector<hash128_t> &bucket, typename RiceBitVector<AT>::Builder &builder, vector<uint32_t> &unary, vector<hash128_t> &bumped) {
		const auto m = bucket.size();
		vector<hash128_t> temp(m, hash128_t(0,0));
		return recSplit(bucket, temp, 0, m, builder, unary, bumped, 0);
	}

	void recSplit(vector<hash128_t> &bucket, vector<hash128_t> &temp, size_t start, size_t end, typename RiceBitVector<AT>::Builder &builder, vector<uint32_t> &unary, vector<hash128_t> &bumped, const int level) {
		const auto m = end - start;
		assert(m > 0);
		uint64_t x = start_seed[level];

		if (m <= _k) {
#ifdef MORESTATS
			sum_depths += m * level;
#endif
			if (m != _k) {
				bumped.insert(bumped.end(),
				  bucket.begin() + start, bucket.begin() + end);
			}
		} else {
#ifdef MORESTATS
			auto start_time = high_resolution_clock::now();
#endif
			if (m > upper_aggr) { // fanout = 2
				const size_t split = ((uint16_t(m / 2 + upper_aggr - 1) / upper_aggr)) * upper_aggr;

				size_t count[2];
				for (;;) {
					count[0] = 0;
					for (size_t i = start; i < end; i++) {
						count[remap16(remix(bucket[i].second + x), m) >= split]++;
#ifdef MORESTATS
						++num_split_evals;
#endif
					}
					if (count[0] == split) break;
					x++;
				}

				count[0] = 0;
				count[1] = split;
				for (size_t i = start; i < end; i++) {
					temp[count[remap16(remix(bucket[i].second + x), m) >= split]++] = bucket[i];
				}
				copy(&temp[0], &temp[m], &bucket[start]);
				x -= start_seed[level];
				const auto log2golomb = golomb_param(m);
				builder.appendFixed(x, log2golomb);
				unary.push_back(x >> log2golomb);

#ifdef MORESTATS
				time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
				recSplit(bucket, temp, start, start + split, builder, unary, bumped, level + 1);
				recSplit(bucket, temp, start + split, end, builder, unary, bumped, level + 1);
			} else if (m > lower_aggr) { // 2nd aggregation level
				const size_t fanout = uint16_t(m + lower_aggr - 1) / lower_aggr;
				size_t count[fanout]; // Note that we never read count[fanout-1]
				for (;;) {
					memset(count, 0, sizeof count - sizeof *count);
					for (size_t i = start; i < end; i++) {
						count[uint16_t(remap16(remix(bucket[i].second + x), m)) / lower_aggr]++;
#ifdef MORESTATS
						++num_split_evals;
#endif
					}
					size_t broken = 0;
					for (size_t i = 0; i < fanout - 1; i++) broken |= count[i] - lower_aggr;
					if (!broken) break;
					x++;
				}

				for (size_t i = 0, c = 0; i < fanout; i++, c += lower_aggr) count[i] = c;
				for (size_t i = start; i < end; i++) {
					temp[count[uint16_t(remap16(remix(bucket[i].second + x), m)) / lower_aggr]++] = bucket[i];
				}
				copy(&temp[0], &temp[m], &bucket[start]);

				x -= start_seed[level];
				const auto log2golomb = golomb_param(m);
				builder.appendFixed(x, log2golomb);
				unary.push_back(x >> log2golomb);

#ifdef MORESTATS
				time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
				size_t i;
				for (i = 0; i < m - lower_aggr; i += lower_aggr) {
					recSplit(bucket, temp, start + i, start + i + lower_aggr, builder, unary, bumped, level + 1);
				}
				recSplit(bucket, temp, start + i, end, builder, unary, bumped, level + 1);
			} else if (m > leaf_keys) { // First aggregation level, m <= lower_aggr
				const size_t fanout = uint16_t(m + leaf_keys - 1) / leaf_keys;
				size_t count[fanout]; // Note that we never read count[fanout-1]
				for (;;) {
					memset(count, 0, sizeof count - sizeof *count);
					for (size_t i = start; i < end; i++) {
						count[uint16_t(remap16(remix(bucket[i].second + x), m)) / leaf_keys]++;
#ifdef MORESTATS
						++num_split_evals;
#endif
					}
					size_t broken = 0;
					for (size_t i = 0; i < fanout - 1; i++) broken |= count[i] - leaf_keys;
					if (!broken) break;
					x++;
				}
				for (size_t i = 0, c = 0; i < fanout; i++, c += leaf_keys) count[i] = c;
				for (size_t i = start; i < end; i++) {
					temp[count[uint16_t(remap16(remix(bucket[i].second + x), m)) / leaf_keys]++] = bucket[i];
				}
				copy(&temp[0], &temp[m], &bucket[start]);

				x -= start_seed[level];
				const auto log2golomb = golomb_param(m);
				builder.appendFixed(x, log2golomb);
				unary.push_back(x >> log2golomb);

#ifdef MORESTATS
				time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
				size_t i;
				for (i = 0; i < m - leaf_keys; i += leaf_keys) {
					recSplit(bucket, temp, start + i, start + i + leaf_keys, builder, unary, bumped, level + 1);
				}
				recSplit(bucket, temp, start + i, end, builder, unary, bumped, level + 1);
			} else { // Leaf
				const size_t fanout = uint16_t(m + _k - 1) / _k;
				size_t count[fanout]; // Note that we never read count[fanout-1]
				for (;;) {
					memset(count, 0, sizeof count - sizeof *count);
					for (size_t i = start; i < end; i++) {
						count[uint16_t(remap16(remix(bucket[i].second + x), m)) / _k]++;
#ifdef MORESTATS
						++num_split_evals;
#endif
					}
					size_t broken = 0;
					for (size_t i = 0; i < fanout - 1; i++) broken |= count[i] - _k;
					if (!broken) break;
					x++;
				}
				for (size_t i = 0, c = 0; i < fanout; i++, c += _k) count[i] = c;
				for (size_t i = start; i < end; i++) {
					temp[count[uint16_t(remap16(remix(bucket[i].second + x), m)) / _k]++] = bucket[i];
				}
				copy(&temp[0], &temp[m], &bucket[start]);

				x -= start_seed[level];
				const auto log2golomb = golomb_param(m);
				builder.appendFixed(x, log2golomb);
				unary.push_back(x >> log2golomb);

#ifdef MORESTATS
				time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
				size_t i;
				for (i = 0; i < m - _k; i += _k) {
					recSplit(bucket, temp, start + i, start + i + _k, builder, unary, bumped, level + 1);
				}
				recSplit(bucket, temp, start + i, end, builder, unary, bumped, level + 1);
			}

#ifdef MORESTATS
			++split_count;
			num_split_trials += x + 1;
			double e_trials = 1;
			size_t aux = m;
			SplitStrat strat{m};
			auto v = strat.begin();
			for (int i = 0; i < strat.fanout(); ++i, ++v) {
				e_trials *= pow((double)m / *v, *v);
				for (size_t j = *v; j > 0; --j, --aux) {
					e_trials *= (double)j / aux;
				}
			}
			expected_split_trials += (size_t)e_trials;
			expected_split_evals += (size_t)e_trials * m;
			const auto log2golomb = golomb_param(m);
			split_unary += 1 + (x >> log2golomb);
			split_fixed += log2golomb;

			min_split_code = min(min_split_code, x);
			max_split_code = max(max_split_code, x);
			sum_split_codes += x;

			auto b = split_golomb_b<K, LEAF_SIZE>(m);
			auto log2b = lambda(b);
			split_unary_golomb += x / b + 1;
			split_fixed_golomb += x % b < ((1ULL << log2b + 1) - b) ? log2b : log2b + 1;
#endif
		}
	}

	void hash_gen(hash128_t *hashes) {
#ifdef MORESTATS
		memset(time_split, 0, sizeof time_split);
		split_unary = split_fixed = 0;
		min_split_code = 1UL << 63;
		max_split_code = sum_split_codes = 0;
		sum_depths = 0;
		split_unary_golomb = split_fixed_golomb = 0;
		size_t minsize = keys_count, maxsize = 0;
		double ub_split_bits = 0;
		double ub_split_evals = 0;
#endif

#ifndef __SIZEOF_INT128__
		if (keys_count > (1ULL << 32)) {
			fprintf(stderr, "For more than 2^32 keys, you need 128-bit integer support.\n");
			abort();
		}
#endif
		nbuckets = max(1, (keys_count + bucket_size - 1) / bucket_size);
		auto bucket_size_acc = vector<int64_t>(nbuckets + 1);
		auto bucket_pos_acc = vector<int64_t>(nbuckets + 1);
		if constexpr (BS == BumpStrategy::RANK_RECURSE
		  || BS == BumpStrategy::RANK_SPLIT
		  || BS == BumpStrategy::INTERSPERSE_SPLIT
		  || BS == BumpStrategy::SELECT_RECURSE) {
			this->buckets_bumped.size((nbuckets + 1 + 63) / 64);
		}

		sort(hashes, hashes + keys_count, [this](const hash128_t &a, const hash128_t &b) { return hash128_to_bucket(a) < hash128_to_bucket(b); });
		typename RiceBitVector<AT>::Builder builder;

		bucket_size_acc[0] = bucket_pos_acc[0] = 0;
		vector<hash128_t> bumped_keys;
		if constexpr (BS == BumpStrategy::RANK_RECURSE
		  || BS == BumpStrategy::RANK_SPLIT
		  || BS == BumpStrategy::BINS_MOD_RECURSE) {
			this->unbumped_bins = 0;
		}
		if constexpr (BS == BumpStrategy::BINS_MOD_RECURSE) {
			this->bucket_mods.size((nbuckets * bit_width(_k - 1) + 7) / 8 + 7);
		}
		for (size_t i = 0, last = 0; i < nbuckets; i++) {
			vector<hash128_t> bucket;
			for (; last < keys_count && hash128_to_bucket(hashes[last]) == i; last++) bucket.push_back(hashes[last]);

			const size_t s = bucket.size();
			assert((s < MAX_BUCKET_SIZE<K, LEAF_SIZE>));
			if constexpr (BS == BumpStrategy::BINS_MOD_RECURSE) {
				bucket_size_acc[i + 1] = bucket_size_acc[i] + s / _k;
				size_t rem = s % _k;
				size_t bit = i * bit_width(_k - 1);
				uint64_t word;
				memcpy(&word, &this->bucket_mods + bit/8, 8);
				word |= rem << (bit&7);
				memcpy(&this->bucket_mods + bit/8, &word, 8);
			} else {
				bucket_size_acc[i + 1] = bucket_size_acc[i] + s;
			}
			if constexpr (BS == BumpStrategy::RANK_RECURSE
			  || BS == BumpStrategy::RANK_SPLIT
			  || BS == BumpStrategy::BINS_MOD_RECURSE) {
				this->unbumped_bins += s / _k;
			}
			vector<uint32_t> unary;
			recSplit(bucket, builder, unary, bumped_keys);
			builder.appendUnaryAll(unary);
			if constexpr (BS == BumpStrategy::RANK_RECURSE
			  || BS == BumpStrategy::RANK_SPLIT
			  || BS == BumpStrategy::INTERSPERSE_SPLIT
			  || BS == BumpStrategy::SELECT_RECURSE) {
				if (bucket_size_acc[i] % _k > bucket_size_acc[i+1] % _k) {
					this->buckets_bumped[i / 64] |= uint64_t(1) << (i&63);
				}
			}
			if constexpr (BS == BumpStrategy::RANK_SPLIT
			  || BS == BumpStrategy::INTERSPERSE_SPLIT
			  || BS == BumpStrategy::INTERSPERSE_SPLIT_NOBITS) {
				const size_t m = bumped_keys.size();
				const size_t split = _k - bucket_size_acc[i] % _k;
				if (split < m) {
					uint64_t x = end(start_seed)[-1];

					size_t count[2];
					for (;;) {
						count[0] = 0;
						for (hash128_t item: bumped_keys) {
							count[remap16(remix(item.second + x), m) >= split]++;
						}
						if (count[0] == split) break;
						x++;
					}

					x -= end(start_seed)[-1];
					const double p = sqrt(m / (2 * M_PI * split * (m - split)));
					builder.appendGolomb(x, p);
				}
				bumped_keys.clear();
			}
			bucket_pos_acc[i + 1] = builder.getBits();
#ifdef MORESTATS
			auto upper_leaves = (s + _k - 1) / _k;
			auto upper_height = ceil(log(upper_leaves) / log(2)); // TODO: check
			auto upper_s = _k * pow(2, upper_height);
			ub_split_bits += (double)upper_s / (_k * 2) * log2(2 * M_PI * _k) - .5 * log2(2 * M_PI * upper_s);
			ub_split_evals += 4 * upper_s * sqrt(pow(2 * M_PI * upper_s, 2 - 1) / pow(2, 2));
			minsize = min(minsize, s);
			maxsize = max(maxsize, s);
#endif
		}
		if constexpr (BS == BumpStrategy::RANK_RECURSE
		  || BS == BumpStrategy::RANK_SPLIT
		  || BS == BumpStrategy::INTERSPERSE_SPLIT
		  || BS == BumpStrategy::SELECT_RECURSE) {
			this->buckets_bumped[nbuckets / 64] |=
			  uint64_t(1) << (nbuckets & 63);
		}
		builder.appendFixed(1, 1); // Sentinel (avoids checking for parts of size 1)
		descriptors = builder.build();
		ef = DoubleEF<AT>(vector<uint64_t>(bucket_size_acc.begin(), bucket_size_acc.end()), vector<uint64_t>(bucket_pos_acc.begin(), bucket_pos_acc.end()));
		if constexpr (BS == BumpStrategy::RANK_RECURSE
		  || BS == BumpStrategy::RANK_SPLIT) {
			this->buckets_bumped_rank =
			  new Rank9<AT>(&this->buckets_bumped, nbuckets+1);
		} else if constexpr (BS == BumpStrategy::SELECT_RECURSE) {
			this->buckets_bumped_select =
			  SimpleSelectHalf<AT>(&this->buckets_bumped, nbuckets+1);
		}
		if constexpr (BS == BumpStrategy::RANK_RECURSE
		  || BS == BumpStrategy::BINS_MOD_RECURSE
		  || BS == BumpStrategy::SELECT_RECURSE) {
			if (bumped_keys.size() > _k) {
				this->recurse_bumped =
				  make_unique<RecSplit<K, LEAF_SIZE, BS, AT>>(bumped_keys,
				                                              bucket_size);
			}
		}

#ifdef STATS
		// Evaluation purposes only
		double total = 0;
		auto show = [&](std::string name, double bits) {
			bits /= keys_count;
			name.push_back(':');
			std::cout << format("{:<25}{:.6f} bits/key\n", name, bits);
			total += bits;
		};
		printf("Elias-Fano cumul sizes:  %f bits/bucket\n", (double)ef.bitCountCumKeys() / nbuckets);
		printf("Elias-Fano cumul bits:   %f bits/bucket\n", (double)ef.bitCountPosition() / nbuckets);
		show("Elias-Fano cumul sizes", ef.bitCountCumKeys());
		show("Elias-Fano cumul bits", ef.bitCountPosition());
		show("Rice-Golomb descriptors", builder.getBits());
		this->stats(show);
		show("Data structure", sizeof(*this) * 8);
		printf("Total bits:              %f bits/key\n", total);
#endif
#ifdef MORESTATS

		printf("\n");
		printf("Min bucket size: %lu\n", minsize);
		printf("Max bucket size: %lu\n", maxsize);

		printf("\n");
		printf("Split count:       %16zu\n", split_count);

		printf("Total split evals: %16lld\n", num_split_evals);

		printf("\n");
		printf("Average depth:        %f\n", (double)sum_depths / keys_count);
		printf("\n");
		printf("Trials per split:     %16.3f\n", (double)num_split_trials / split_count);
		printf("Exp trials per split: %16.3f\n", (double)expected_split_trials / split_count);
		printf("Evals per split:      %16.3f\n", (double)num_split_evals / split_count);
		printf("Exp evals per split:  %16.3f\n", (double)expected_split_evals / split_count);

		printf("\n");
		printf("Unary bits per split: %10.5f\n", (double)split_unary / split_count);
		printf("Fixed bits per split: %10.5f\n", (double)split_fixed / split_count);
		printf("Total bits per split: %10.5f\n", (double)(split_unary + split_fixed) / split_count);
		printf("Total bits per key:   %10.5f\n", (double)(split_unary + split_fixed) / keys_count);

		printf("\n");
		printf("Unary bits per split (Golomb): %10.5f\n", (double)split_unary_golomb / split_count);
		printf("Fixed bits per split (Golomb): %10.5f\n", (double)split_fixed_golomb / split_count);
		printf("Total bits per split (Golomb): %10.5f\n", (double)(split_unary_golomb + split_fixed_golomb) / split_count);
		printf("Total bits per key (Golomb):   %10.5f\n", (double)(split_unary_golomb + split_fixed_golomb) / keys_count);

		printf("\n");

		printf("Total split bits        %16.3f\n", (double)split_fixed + split_unary);
		printf("Upper bound split bits: %16.3f\n", ub_split_bits);
#endif
	}

#if 0
	friend ostream &operator<<(ostream &os, const RecSplit<K, BS, AT> &rs) {
		const size_t k = K;
		os.write((char *)&k, sizeof(k));
		os.write((char *)&rs.bucket_size, sizeof(rs.bucket_size));
		os.write((char *)&rs.keys_count, sizeof(rs.keys_count));
		os << rs.descriptors;
		os << rs.ef;
		return os;
	}

	friend istream &operator>>(istream &is, RecSplit<K, BS, AT> &rs) {
		size_t k;
		is.read((char *)&k, sizeof(k));
		if (k != K) {
			fprintf(stderr, "Serialized k %d, code k %d\n", int(k), int(K));
			abort();
		}
		is.read((char *)&rs.bucket_size, sizeof(bucket_size));
		is.read((char *)&rs.keys_count, sizeof(keys_count));
		rs.nbuckets = max(1, (rs.keys_count + rs.bucket_size - 1) / rs.bucket_size);

		is >> rs.descriptors;
		is >> rs.ef;
		return is;
	}
#endif
};

} // namespace sux::function::krecsplit
