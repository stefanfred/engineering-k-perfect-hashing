#pragma once

#include <vector>
#include <limits>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <ranges>

#include <sux/bits/EliasFano.hpp>
#include <ips2ra.hpp>
#include <bytehamster/util/Function.h>
#include <Fips.h>

#include <hash128.hpp>
#include "common.hpp"
#include "consensus.hpp"

#include "optimalThresholds.hpp"

namespace kphf::ThresholdBasedBumpingConsensus {

struct ThresholdInfo {
    double threshold;
    double l1, l2;

    ThresholdInfo(double t): threshold(t), l1(log(t)), l2(log(1-t)) {}
};

inline std::vector<ThresholdInfo> compute_threshold_info(const std::vector<double> &thresholds) {
    std::vector<ThresholdInfo> threshold_info;
    for (double t: thresholds) {
        threshold_info.emplace_back(t);
    }
    return threshold_info;
}

inline std::vector<double> compute_lgamma_cache(uint64_t limit) {
    std::vector<double> r(limit + 2);
    for (size_t i = 0; i < r.size(); i++) {
        r[i] = lgamma(i);
    }
    return r;
}

inline double success_exp(uint64_t n, uint64_t k, const std::vector<ThresholdInfo> &thresholds,
        const std::vector<double> &lgamma_cache) {
    if (k > n) return 0;
    double lbinom = lgamma_cache[n + 1] - lgamma_cache[k + 1] - lgamma_cache[n - k + 1];
    double e = 0;
    for (const auto &[threshold, l1, l2] : thresholds) {
        if (std::isinf(l1)) {
            e += (k == 0);
        } else if (std::isinf(l2)) {
            e += (k == n);
        } else {
            e += exp(lbinom + k * l1 + (n - k) * l2);
        }
    }
    return e;
}

inline double best_error(uint64_t n, size_t k, const std::vector<ThresholdInfo> &thresholds, const std::vector<double> &lgamma_cache) {
    uint64_t goal = std::min(n, k);
    constexpr double TUNE = 1.05;
    double need = TUNE;
    need -= success_exp(n, goal, thresholds, lgamma_cache);
    if (need <= 0) {
        return 0;
    }
    double avg = 0;
    for (uint64_t error = 1; error <= goal; error++) {
        double s = success_exp(n, goal-error, thresholds, lgamma_cache);
        if (s >= need) {
            return error - 1 + need/s;
        }
        need -= s;
        avg += error * s;
    }
    return goal;
}

struct AllowedImperfectBumpingEntry {
    uint64_t allowedBumped;
    uint64_t probability;
};

template <uint64_t k, int threshold_size>
struct AllowedImperfectBumping {
    static constexpr uint64_t n_thresholds = 1ul << threshold_size;
    const uint64_t max_n;
    std::vector<double> thresholds;
    std::vector<ThresholdInfo> threshold_info;
    std::vector<double> lgamma_cache;
    std::vector<AllowedImperfectBumpingEntry> errors;

    AllowedImperfectBumping(double lambda)
        : max_n(static_cast<uint64_t>(2ul * k * lambda)) {
        thresholds = ThresholdBasedBumping::compute_thresholds_normalized(k, lambda*k, n_thresholds);
        threshold_info = compute_threshold_info(thresholds);
        lgamma_cache = compute_lgamma_cache(max_n);
        errors.resize(max_n + 1, {~0ul, 0});
    }

    AllowedImperfectBumpingEntry getAllowed(size_t actualBucketSize) {
        if (errors[actualBucketSize].allowedBumped == ~0ul) {
            double e = best_error(actualBucketSize, k, threshold_info, lgamma_cache);
            uint64_t q = static_cast<uint64_t>(e);
            errors[actualBucketSize] = {q+1, ThresholdBasedBumping::double_to_u64(e-q)};
        }
        return errors[actualBucketSize];
    }

    size_t maxBucketSize() {
        return max_n;
    }
};

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

template <uint64_t k, int threshold_size>
class ThresholdBasedBumpingConsensus {
private:
    static_assert(threshold_size > 1);
    static_assert(k > 0);
    static constexpr uint64_t n_thresholds = 1ul << threshold_size;
    uint64_t n;
    std::array<uint64_t, n_thresholds> avail_thresholds;
    std::vector<std::pair<uint64_t, Consensus>> layers;
    fips::FiPS<> phf;
    mutable sux::bits::EliasFano<> gaps;
    using Key = kphf::ThresholdBasedBumping::Key;
public:

    uint64_t operator()(const std::string &key) const {
        return operator()(Hash128(key));
    }

    uint64_t operator()(const Hash128 hash) const {
        uint64_t offset = 0;
        for (uint64_t i = 0; i < layers.size(); i++) {
            auto &[cur_buckets, consensus] = layers[i];
            uint64_t b = calculateBucket(hash, i, cur_buckets);
            auto [seed, tidx] = consensus.get(b * threshold_size, threshold_size);
            tidx = decrypt(seed, tidx);
            uint64_t f = bytehamster::util::remix(hash.lo + seed + i);

            if (f < avail_thresholds[tidx]) {
                return offset + b;
            }

            offset += cur_buckets;
        }

        return gaps.select(phf(hash.hi ^ hash.lo));
    }

    size_t count_bits() const {
        size_t s = 0;
        for (auto &[nb,c] : layers) {
            s += c.count_bits() - sizeof(c) * 8;
        }
        return sizeof(*this) * 8
            - sizeof(avail_thresholds) * 8 // Lookup table independent of keys, could actually be static constexpr
            + layers.capacity() * sizeof(layers[0])
            + s
            + phf.getBits() - sizeof(phf) * 8
            + gaps.bitCount() - sizeof(gaps) * 8;
    }

private:
    static inline size_t calculateBucket(Hash128 key, size_t layer, size_t cur_buckets) {
        return bytehamster::util::fastrange64(bytehamster::util::remix(key.hi + layer), cur_buckets);
    }

    static inline uint64_t encrypt(uint64_t key, uint64_t tidx) {
        return (tidx ^ bytehamster::util::remix(key)) & ((1ul << threshold_size) - 1);
    }

    static inline uint64_t decrypt(uint64_t key, uint64_t tidx) {
        return (tidx ^ bytehamster::util::remix(key)) & ((1ul << threshold_size) - 1);
    }

    struct LayerBuilder {
        std::vector<Key> &keys;
        std::ranges::subrange<std::vector<Key>::iterator> current;
        uint64_t cur_bucket, total_buckets, layer, offset;
        std::vector<uint64_t> &emptySlots;
        std::vector<Hash128> bumped;

        const std::array<uint64_t, n_thresholds> &thresholds;
        AllowedImperfectBumping<k, threshold_size> &allowedImperfectBumping;

        LayerBuilder(const std::array<uint64_t, n_thresholds> &thresholds,
                    AllowedImperfectBumping<k, threshold_size> &allowedImperfectBumping,
                    std::vector<Key> &keys, uint64_t layer, uint64_t buckets,
                    uint64_t offset, std::vector<uint64_t> &emptySlots)
                : keys(keys), current(keys.begin(), keys.begin()),
                  cur_bucket(0), total_buckets(buckets), layer(layer), offset(offset),
                  emptySlots(emptySlots), thresholds(thresholds), allowedImperfectBumping(allowedImperfectBumping) {
        }

        std::optional<uint64_t> advance(uint64_t seed) {
            if (cur_bucket == total_buckets) {
                return {};
            }

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
            while (!emptySlots.empty() && emptySlots.back() == offset + cur_bucket - 1) {
                emptySlots.pop_back();
                cnt++;
            }
            bumped.resize(bumped.size() - (current.size() + cnt - k));

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
            return find(seed, decrypt(seed, prev));
        }

        std::optional<uint64_t> find(uint64_t seed, uint64_t prev) {
            uint64_t goal = std::min(k, current.size());
            AllowedImperfectBumpingEntry allowed = {};
            if (current.size() < allowedImperfectBumping.maxBucketSize()) {
                allowed = allowedImperfectBumping.getAllowed(current.size());
            } else {
                allowed = {k + 1, 0};
            }
            for (uint64_t tidx = prev; tidx > 0;) {
                tidx--;
                auto it = std::ranges::lower_bound(current, thresholds[tidx], {},
                    [](const Key &key) { return key.fingerprint; });
                uint64_t idx = it - current.begin();
                uint64_t error = goal - idx;
                if (error > allowed.allowedBumped) {
                    break;
                }
                if (error < allowed.allowedBumped || bytehamster::util::remix(seed + tidx) < allowed.probability) {
                    for (Key & key: current | std::views::drop(idx)) {
                        bumped.push_back(key.hash);
                    }
                    for (uint64_t x = idx; x < k; x++) {
                        emptySlots.push_back(offset + cur_bucket - 1);
                    }
                    return encrypt(seed, tidx);
                }
            }
            return {};
        }
    };

    public:
    ThresholdBasedBumpingConsensus() : n(0), gaps({}, 0) {
    }

    explicit ThresholdBasedBumpingConsensus(const std::vector<std::string> &keys, double overload)
            : ThresholdBasedBumpingConsensus(hashKeys(keys, overload), overload) {
    }

    static std::vector<Key> hashKeys(const std::vector<std::string> &keys, double overload) {
        uint64_t total_buckets = (keys.size() + k - 1) / k;
        double overload_bucket_size = k * overload;
        uint64_t cur_buckets = std::ceil(keys.size() / overload_bucket_size);
        if (cur_buckets >= total_buckets) {
            cur_buckets = total_buckets;
        } else if ((double)(keys.size() - k * cur_buckets) / (total_buckets - cur_buckets) > overload_bucket_size) {
            cur_buckets = total_buckets;
        }

        std::vector<Key> hashed_keys;
        hashed_keys.reserve(keys.size());
        for (const auto &key : keys) {
            Hash128 hash(key);
            hashed_keys.emplace_back(calculateBucket(hash, 0, cur_buckets), 0, hash);
        }
        return hashed_keys;
    }

    explicit ThresholdBasedBumpingConsensus(std::vector<Key> keys, double overload)
            : n(keys.size()), gaps({}, 0),
              avail_thresholds(ThresholdBasedBumping::compute_thresholds<n_thresholds>(k, k * overload)) {
        AllowedImperfectBumping<k, threshold_size> allowedImperfectBumping(overload);

        double overload_bucket_size = k * overload;
        uint64_t total_buckets = (n + k - 1) / k;

#ifdef STATS
        total_keys += n;
#endif

        uint64_t offset = 0;
        std::vector<uint64_t> emptySlots;
        std::vector<Hash128> bumped;
        for (uint64_t layer = 0; offset != total_buckets; layer++) {
            uint64_t remaining = total_buckets - offset;
            size_t cur_keys = layer > 0 ? bumped.size() : keys.size();
            uint64_t cur_buckets = std::ceil(cur_keys / overload_bucket_size);
            if (cur_buckets >= remaining) {
                cur_buckets = remaining;
            } else if ((double)(cur_keys - k * cur_buckets) / (remaining - cur_buckets) > overload_bucket_size) {
                cur_buckets = remaining;
            }

            if (layer > 0) {
                keys.clear();
                keys.reserve(bumped.size());
                for (const Hash128 &key : bumped) {
                    keys.emplace_back(calculateBucket(key, layer, cur_buckets), 0, key);
                }
                bumped.clear();
            }

            sort_buckets(keys);

            LayerBuilder builder(avail_thresholds, allowedImperfectBumping, keys, layer, cur_buckets, offset, emptySlots);
            Consensus consensus(builder);
            layers.emplace_back(cur_buckets, std::move(consensus));
            std::swap(bumped, builder.bumped);

            offset += cur_buckets;
        }
        layers.shrink_to_fit();

#ifdef STATS
        bumped_keys += bumped.size();
#endif

        std::vector<uint64_t> fallbackKeys;
        for (Hash128 key : bumped) {
            fallbackKeys.push_back(key.hi ^ key.lo);
        }
        phf = fips::FiPS<>(fallbackKeys);
        std::vector<uint64_t> actual_spots;
        for (uint64_t key : fallbackKeys) {
            uint64_t h = phf(key);
            if (h >= actual_spots.size()) {
                actual_spots.resize(h + 1);
            }
            actual_spots[h] = 1;
        }
        auto it = emptySlots.begin();
        for (uint64_t &x : actual_spots) {
            uint64_t v = *it;
            if (x) {
                ++it;
            }
            x = v;
        }
        gaps = sux::bits::EliasFano<>(actual_spots, total_buckets);
    }
};

}

