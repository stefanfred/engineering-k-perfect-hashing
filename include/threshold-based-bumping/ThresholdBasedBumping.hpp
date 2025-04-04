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

#include <SimpleRibbon.h>
#include <ips2ra.hpp>
#include <hash128.hpp>
#include <sux/bits/EliasFano.hpp>
#include <bytehamster/util/Function.h>
#include <Fips.h>

#include "optimalThresholds.hpp"

#include "common.hpp"

namespace kphf::ThresholdBasedBumping {

template<uint64_t n_thresholds>
std::array<uint64_t, n_thresholds> get_thresholds(uint64_t _k, double bucket_size) {
    auto res = compute_thresholds(_k,bucket_size,n_thresholds);
    std::array<uint64_t, n_thresholds> raw;
    for (size_t i = 0; i < n_thresholds; ++i) {
        raw[i] = double_to_u64(res[i]);
    }
    return raw;
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

class DummyFilter {
private:
    DummyFilter() {}

public:
    class Builder {
    public:
        Builder() {}

        void add(int level, std::vector<Hash128> &elems, size_t bump) {
            (void) level, (void) elems, (void) bump;
        }

        DummyFilter build() {
            return DummyFilter();
        }

    };

    bool bump(int level, Hash128 hash) const {
        (void) level, (void) hash;
        return true;
    }

    size_t count_bits() const {
        return sizeof(*this) * 8;
    }

    static std::string name() {
        return "none";
    }
};

class RibbonFilter {
private:
    using Ribbon = SimpleRibbon<1>;

    mutable Ribbon ribbon;
    RibbonFilter(Ribbon &&src): ribbon(std::move(src)) {}

    static uint64_t hash(Hash128 h, int level) {
        return h.lo + uint64_t(level);
    }

public:
    class Builder {
    private:
        std::vector<std::pair<uint64_t, uint8_t>> items;

    public:
        Builder() {}

        void add(int level, std::vector<Hash128> &elems, size_t bump) {
            for (size_t i = 0; i < elems.size(); i++) {
                items.emplace_back(hash(elems[i], level), i < bump);
            }
            elems.resize(bump);
        }

        RibbonFilter build() {
            return RibbonFilter(std::move(Ribbon(items)));
        }

    };

    bool bump(int level, Hash128 h) const {
        return ribbon.retrieve(hash(h, level));
    }

    size_t count_bits() const {
        return ribbon.sizeBytes() * 8;
    }

    static std::string name() {
        return "ribbon";
    }
};

template<uint64_t K, double OVERLOAD, int THRESHOLD_SIZE_HALFBITS, typename Filter = RibbonFilter>
class ThresholdBasedBumping {
    static_assert(OVERLOAD > 1.0);
    static_assert(THRESHOLD_SIZE_HALFBITS >= 2);
    static_assert(K > 0);

private:
    static constexpr uint64_t _k = K;
    static constexpr double overload = OVERLOAD;
    static constexpr uint64_t threshold_size_halfbits = THRESHOLD_SIZE_HALFBITS;
    static constexpr uint64_t n_regions =
      uint64_t(sqrt(uint64_t(1) << threshold_size_halfbits));
    static constexpr uint64_t n_thresholds = n_regions - 1;
    static constexpr double overload_bucket_size = _k * overload;

    std::array<uint64_t, n_thresholds> avail_thresholds = get_thresholds<n_thresholds>(_k, overload_bucket_size);

    uint64_t n;
    std::vector<uint64_t> nbuckets;
    std::vector<uint8_t> thresholds;
    [[no_unique_address]] Filter filter;
    fips::FiPS<> phf;
    mutable sux::bits::EliasFano<> gaps;

    ThresholdBasedBumping(uint64_t n, std::vector<uint64_t> &&nbuckets, std::vector<uint8_t> &&thresholds,
        Filter &&filter, fips::FiPS<> &&phf, sux::bits::EliasFano<> &&gaps)
            : n(n), nbuckets(std::move(nbuckets)),
              thresholds(std::move(thresholds)), filter(std::move(filter)),
              phf(std::move(phf)), gaps(std::move(gaps)) {
    }

public:
    ThresholdBasedBumping()
        : n(0), filter(typename Filter::Builder().build()), gaps(std::vector<uint64_t>(), 0) {
    }

    template<typename F>
    ThresholdBasedBumping(const std::vector<Hash128> &keys, F &&build_phf, Filter::Builder filter = {})
        : ThresholdBasedBumping(build(keys, std::forward<F>(build_phf), std::move(filter))) {
    }

    uint64_t operator()(Hash128 key) const {
        uint64_t offset = 0;
        for (uint64_t i = 0; i < nbuckets.size(); i++) {
            uint64_t cur_buckets = nbuckets[i];
            uint64_t h = bytehamster::util::remix(key.hi + i);
            uint64_t b = offset + bytehamster::util::fastrange64(h, cur_buckets);
            uint64_t f = bytehamster::util::remix(key.lo + i);
            uint64_t tidx;
            if constexpr (threshold_size_halfbits & 1) {
                uint64_t off = (b/2) * threshold_size_halfbits;
                memcpy(&tidx, thresholds.data() + off/8, 8);
                tidx >>= off%8;
                tidx &= (uint64_t(1) << threshold_size_halfbits) - 1;
                if (b%2 == 0) tidx %= n_regions;
                else tidx /= n_regions;
            } else {
                uint64_t off = b * (threshold_size_halfbits/2);
                memcpy(&tidx, thresholds.data() + off/8, 8);
                tidx >>= off%8;
                tidx &= (uint64_t(1) << (threshold_size_halfbits/2)) - 1;
            }

            if (tidx != 0 && f < avail_thresholds[tidx-1]) return b;
            if (tidx == n_thresholds || f < avail_thresholds[tidx]) {
                if (!filter.bump(i, key)) return b;
            }

            offset += cur_buckets;
        }

        return gaps.select(phf(key.hi ^ key.lo));
    }

    size_t count_bits() const {
        return sizeof(*this) * 8
            + nbuckets.capacity() * 64
            + thresholds.capacity() * 8
            + filter.count_bits() - sizeof(filter) * 8
            + phf.getBits() - sizeof(phf) * 8
            + gaps.bitCount() - sizeof(gaps) * 8
        ;
    }

    ThresholdBasedBumping(std::vector<Hash128> keys)
            : n(keys.size()), filter(typename Filter::Builder().build()), gaps(std::vector<uint64_t>(), 0) {
        typename Filter::Builder filter;
        uint64_t total_buckets = (n + _k - 1) / _k;
        thresholds.resize(((total_buckets+1)/2 * threshold_size_halfbits + 7) / 8 + 7);
        thresholds.shrink_to_fit();

#ifdef STATS
        total_keys += n;
#endif

        uint64_t offset = 0;
        std::vector<uint64_t> spots;
        for (uint64_t i = 0; offset != total_buckets; i++) {
            uint64_t cur_buckets = std::min(total_buckets - offset, uint64_t(std::ceil(keys.size() / overload_bucket_size)));
            nbuckets.push_back(cur_buckets);

            auto r = std::views::transform(keys,
              [i, cur_buckets](Hash128 key) {
                uint64_t b = bytehamster::util::fastrange64(bytehamster::util::remix(key.hi + i), cur_buckets);
                uint64_t f = bytehamster::util::remix(key.lo + i);
                return Key { b, f, key };
            });
            std::vector<Key> sorted_keys(r.begin(), r.end());
            keys.clear();

            sort_buckets(sorted_keys);
            auto next = sorted_keys.begin();
            for (uint64_t j = 0; j < cur_buckets; j++) {
                auto start = next;
                while (next != sorted_keys.end() && next->bucket == j) ++next;
                auto bucket = std::ranges::subrange(start, next);
                sort_fingerprints(bucket);
                uint64_t tidx;
                typename std::vector<Key>::iterator perfect;
                if (bucket.size() <= _k) {
                    tidx = n_thresholds;
                    perfect = bucket.end();
                } else {
                    tidx = std::ranges::upper_bound(avail_thresholds, bucket[_k].fingerprint) - avail_thresholds.begin();
                    perfect = bucket.begin() + _k;
                }

                auto bound = tidx == 0 ? bucket.begin() :
                    std::ranges::lower_bound(
                        bucket,
                        avail_thresholds[tidx-1],
                        std::less<uint64_t>(),
                        [](auto &x) { return x.fingerprint; }
                    );

#ifdef STATS
                total_thresholds++;
                if (bound == perfect) perfect_thresholds++;
                if (bucket.size() >= _k) overfull_buckets++;
#endif

                uint64_t in_bucket;
                /* bound == bucket.begin() may occur if bucket.empty() */
                if (tidx >= 2 && bound == perfect && bound != bucket.begin()
                        && prev(bound)->fingerprint < avail_thresholds[tidx-2]) {
                    tidx--;
                    in_bucket = bound - bucket.begin();
                    for (Key &k: std::ranges::subrange(bound, bucket.end())) {
                        keys.push_back(k.hash);
                    }
                } else {
                    auto upper = bound;
                    std::vector<Hash128> bump;
                    while (upper != bucket.end()
                            && (tidx == n_thresholds || upper->fingerprint < avail_thresholds[tidx])) {
                        bump.push_back(upper->hash);
                        ++upper;
                    }
                    size_t need_bump = upper - perfect;
                    filter.add(i, bump, need_bump);
                    assert(bump.size() >= need_bump);
#ifdef STATS
                    extra_bumped += (uint64_t) (bump.size() - need_bump);
#endif
                    keys.insert(keys.end(), bump.begin(), bump.end());
                    for (Key &k: std::ranges::subrange(upper, bucket.end())) {
                        keys.push_back(k.hash);
                    }
                    in_bucket = (upper - bucket.begin()) - bump.size();
                }

#ifdef STATS
                if (in_bucket == _k) filled_buckets++;
#endif

                uint64_t b = offset + j;
                if constexpr (threshold_size_halfbits&1) {
                    uint64_t off = (b/2) * threshold_size_halfbits;
                    uint64_t word;
                    memcpy(&word, thresholds.data() + off/8, 8);
                    word += (b%2 == 0 ? tidx : tidx * n_regions) << off%8;
                    memcpy(thresholds.data() + off/8, &word, 8);
                } else {
                    uint64_t off = b * (threshold_size_halfbits/2);
                    uint64_t word;
                    memcpy(&word, thresholds.data() + off/8, 8);
                    word |= tidx << off%8;
                    memcpy(thresholds.data() + off/8, &word, 8);
                }

                for (uint64_t c = 0; c < _k - in_bucket; c++) {
                    spots.push_back(offset + j);
                }
            }
            offset += cur_buckets;
        }
        nbuckets.shrink_to_fit();

#ifdef STATS
        bumped_keys += keys.size();
#endif

        std::vector<uint64_t> fallbackKeys;
        for (Hash128 key : keys) {
            fallbackKeys.push_back(key.hi ^ key.lo);
        }
        phf = fips::FiPS<>(fallbackKeys);
        std::vector<uint64_t> actual_spots;
        for (uint64_t key : fallbackKeys) {
            uint64_t h = phf(key);
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
        this->filter = std::move(filter.build());
    }
};

}
