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
#include <Fips.h>

#include <hash128.hpp>
#include "common.hpp"
#include "consensus.hpp"

#include "optimalThresholds.hpp"

namespace kphf::ThresholdBasedBumpingConsensus {

inline double success_prob(double threshold, uint64_t n, uint64_t k) {
    if (k > n) {
        return 0;
    }
    double lbinom = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
    double l1 = log(threshold);
    double l2 = log(1 - threshold);
    if (std::isinf(l1)) {
        return k == 0;
    } else if (std::isinf(l2)) {
        return k == n;
    } else {
        return exp(lbinom + k * l1 + (n-k) * l2);
    }
}

void recomp(uint64_t max_n, const std::vector<double> &thresholds, std::vector<std::vector<double>> &success_exp) {
    for (uint64_t n = 0; n <= max_n; n++) {
        auto &v = success_exp[n];
        for (uint64_t i = 0; i < v.size(); i++) {
            double e = 0;
            for (double t : thresholds) {
                e += success_prob(t, n, i);
            }
            v[i] = e;
        }
    }
}

std::pair<double, double> best_error(uint64_t n, size_t k, const std::vector<std::vector<double>> &success_exp) {
    uint64_t goal = std::min(n, k);
    constexpr double TUNE = 1.05;
    double need = TUNE;
    need -= success_exp[n][goal];
    if (need <= 0) {
        return {0,0};
    }
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
}

template <uint64_t k, int threshold_size>
std::pair<std::vector<uint64_t>, std::vector<std::pair<uint64_t, uint64_t>>> compute_thresholds_and_error(double lamda) {
    uint64_t n_thresholds = 1ul << threshold_size;
    std::vector<double> thresholds = compute_thresholds(k, lamda*k, n_thresholds);

    uint64_t max_n = static_cast<uint64_t>(2ul * k * lamda);
    std::vector<std::vector<double>> success_exp(max_n+1, std::vector<double>(k+1));

    recomp(max_n, thresholds, success_exp);

    std::vector<uint64_t> thresholds_final(n_thresholds);
    for (uint64_t i = 0; i < n_thresholds; i++) {
        thresholds_final[i] = ThresholdBasedBumping::double_to_u64(thresholds[i]);
    }

    std::vector<std::pair<uint64_t, uint64_t>> errors(max_n+1);
    for (uint64_t n = 0; n <= max_n; n++) {
        double e = best_error(n, k, success_exp).first;
        uint64_t q = e;
        errors[n] = {q+1, ThresholdBasedBumping::double_to_u64(e-q)};
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

template <uint64_t k, int threshold_size>
class ThresholdBasedBumpingConsensus {
private:
    static_assert(threshold_size >= 1);
    static_assert(k > 0);
    uint64_t n;
    std::vector<uint64_t> thresholds;
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

            if (f < thresholds[tidx]) {
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

        const std::vector<uint64_t> &thresholds;
        const std::vector<std::pair<uint64_t, uint64_t>> &errors;

        LayerBuilder(const std::vector<uint64_t> &thresholds,
                    const std::vector<std::pair<uint64_t, uint64_t>> &errors,
                    std::vector<Key> &keys, uint64_t layer, uint64_t buckets,
                    uint64_t offset, std::vector<uint64_t> &emptySlots)
                : keys(keys), current(keys.begin(), keys.begin()),
                  cur_bucket(0), total_buckets(buckets), layer(layer), offset(offset),
                  emptySlots(emptySlots), thresholds(thresholds), errors(errors) {
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
            uint64_t error_bound, error_prob;
            if (current.size() < errors.size()) {
                std::tie(error_bound, error_prob) = errors[current.size()];
            } else {
                error_bound = k+1, error_prob = 0;
            }
            for (uint64_t tidx = prev; tidx > 0;) {
                tidx--;
                auto it = std::ranges::lower_bound(current, thresholds[tidx], {},
                    [](const Key &key) { return key.fingerprint; });
                uint64_t idx = it - current.begin();
                uint64_t error = goal - idx;
                if (error > error_bound) {
                    break;
                }
                if (error < error_bound || bytehamster::util::remix(seed + tidx) < error_prob) {
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
            : n(keys.size()), gaps({}, 0) {
        std::vector<std::pair<uint64_t, uint64_t>> errors;
        std::tie(thresholds, errors) = compute_thresholds_and_error<k, threshold_size>(overload);

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
            uint64_t cur_buckets = std::ceil(keys.size() / overload_bucket_size);
            if (cur_buckets >= remaining) {
                cur_buckets = remaining;
            } else if ((double)(keys.size() - k * cur_buckets) / (remaining - cur_buckets) > overload_bucket_size) {
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

            LayerBuilder builder(thresholds, errors, keys, layer, cur_buckets, offset, emptySlots);
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

