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

namespace kphf::ThresholdBasedBumpingConsensus {

inline double poisson(double lamda, uint64_t k) {
    return exp(k * log(lamda) - lamda - lgamma(k+1));
}

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

inline uint64_t double_to_u64(double x) {
    return x == 1.0 ? ~0ul : static_cast<uint64_t>(ldexp(x, 64));
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

void set_threshold(uint64_t max_n, uint64_t i, double value, std::vector<double> &thresholds, std::vector<std::vector<double>> &success_exp) {
    double old = std::exchange(thresholds[i], value);
    for (uint64_t n = 0; n <= max_n; n++) {
        auto &v = success_exp[n];
        for (uint64_t i = 0; i < v.size(); i++) {
            v[i] += success_prob(value, n, i) - success_prob(old, n, i);
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

double eval(uint64_t max_n, size_t k, double lamda, const std::vector<double> &thresholds, const std::vector<std::vector<double>> &success_exp) {
    if (!std::ranges::contains(thresholds, 0.0)) {
        return std::numeric_limits<double>::infinity();
    }
    double exp = k * lamda;
    double badness = 0;
    for (uint64_t n = 0; n <= max_n; n++) {
        badness += poisson(exp, n) * best_error(n, k, success_exp).second;
    }
    return badness;
}

std::pair<std::vector<uint64_t>, std::vector<std::pair<uint64_t, uint64_t>>>
        compute_thresholds(uint64_t k, double lamda, int threshold_size) {

    uint64_t n_thresholds = 1ul << threshold_size;
    std::vector<double> thresholds(n_thresholds);
    for (uint64_t i = 0; i < n_thresholds; i++) {
        thresholds[i] = static_cast<double>(i) / (n_thresholds-1);
    }

    uint64_t max_n = static_cast<uint64_t>(2ul * k * lamda);
    std::vector<std::vector<double>> success_exp(max_n+1, std::vector<double>(k+1));

    double temp = 0.8;
    double cur_badness;
    while (temp > 0.01 / k) {
        recomp(max_n, thresholds, success_exp);
        cur_badness = eval(max_n, k, lamda, thresholds, success_exp);
        for (uint64_t j = 0; j < n_thresholds; j++) {
            double b = thresholds[j];
            double a = std::max(b - temp, 0.0);
            double c = std::min(b + temp, 1.0);
            set_threshold(max_n, j, a, thresholds, success_exp);
            double a_badness = eval(max_n, k, lamda, thresholds, success_exp);
            set_threshold(max_n, j, c, thresholds, success_exp);
            double c_badness = eval(max_n, k, lamda, thresholds, success_exp);
            double best = std::min({cur_badness, a_badness, c_badness});
            if (best == c_badness) {
            } else if (best == a_badness) {
                set_threshold(max_n, j, a, thresholds, success_exp);
            } else if (best == cur_badness) {
                set_threshold(max_n, j, b, thresholds, success_exp);
            }
            cur_badness = best;
        }
        temp *= 0.8;
    }
    std::ranges::sort(thresholds);
    recomp(max_n, thresholds, success_exp);

    std::vector<uint64_t> thresholds_final(n_thresholds);
    for (uint64_t i = 0; i < n_thresholds; i++) {
        thresholds_final[i] = double_to_u64(thresholds[i]);
    }

    std::vector<std::pair<uint64_t, uint64_t>> errors(max_n+1);
    for (uint64_t n = 0; n <= max_n; n++) {
        double e = best_error(n, k, success_exp).first;
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

class ThresholdBasedBumpingConsensus {
private:
    uint64_t n;
    uint64_t k;
    double overload;
    int threshold_size;
    std::vector<uint64_t> thresholds;
    std::vector<std::pair<uint64_t, uint64_t>> errors;
    std::vector<std::pair<uint64_t, Consensus>> layers;
    fips::FiPS<> phf;
    mutable sux::bits::EliasFano<> gaps;
    using Key = kphf::ThresholdBasedBumping::Key;
public:

    uint64_t operator()(Hash128 key) const {
        uint64_t offset = 0;
        for (uint64_t i = 0; i < layers.size(); i++) {
            auto &[cur_buckets, consensus] = layers[i];
            uint64_t h = bytehamster::util::remix(key.hi + i);
            uint64_t b = bytehamster::util::fastrange64(h, cur_buckets);
            auto [seed,tidx] =
              consensus.get(b * threshold_size, threshold_size);
            tidx = decrypt(seed, tidx, threshold_size);
            uint64_t f = bytehamster::util::remix(key.lo + seed + i);

            if (f < thresholds[tidx]) return offset + b;

            offset += cur_buckets;
        }

        return gaps.select(phf(key.hi ^ key.lo));
    }

    size_t count_bits() const {
        size_t s = 0;
        for (auto &[nb,c]: layers) s += c.count_bits() - sizeof(c) * 8;
        return sizeof(*this) * 8
            + layers.capacity() * sizeof(layers[0])
            + s
            + phf.getBits() - sizeof(phf) * 8
            + gaps.bitCount() - sizeof(gaps) * 8;
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

        Builder(uint64_t k, uint64_t threshold_size,
            const std::span<uint64_t> &thresholds,
            const std::span<std::pair<uint64_t, uint64_t>> &errors,
            std::vector<Hash128> &hashes, uint64_t layer, uint64_t buckets,
            uint64_t offset, std::vector<uint64_t> *spots)
                : keys(std::move(prepare(hashes, layer, buckets))),
                  current(keys.begin(), keys.begin()),
                  cur_bucket(0), total_buckets(buckets), layer(layer), offset(offset),
                  spots(spots), bumped(&hashes),
                  k(k), threshold_size(threshold_size),
                  thresholds(thresholds), errors(errors) {
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
                auto it = std::ranges::lower_bound(current, thresholds[tidx], {},
                    [](Key k) -> std::uint64_t { return k.fingerprint; });
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
    ThresholdBasedBumpingConsensus() : n(0), gaps({}, 0) {
    }

    ThresholdBasedBumpingConsensus(size_t k, std::vector<Hash128> &keys, double overload, size_t threshold_size)
            : n(keys.size()), k(k), overload(overload), threshold_size(threshold_size), gaps({}, 0) {
        assert(overload > 1.0);
        assert(threshold_size >= 1);
        assert(k > 0);
        auto begin = std::chrono::high_resolution_clock::now();
        std::cout<<"Begin calculating thresholds"<<std::endl;
        std::tie(thresholds, errors) = compute_thresholds(k, overload, threshold_size);
        std::cout<<"Complete calculating thresholds" << std::endl;
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" << std::endl;

        double overload_bucket_size = k * overload;
        uint64_t total_buckets = (n + k - 1) / k;

#ifdef STATS
        total_keys += n;
#endif

        uint64_t offset = 0;
        std::vector<uint64_t> spots;
        for (uint64_t i = 0; offset != total_buckets; i++) {
            uint64_t remaining = total_buckets - offset;
            uint64_t cur_buckets = std::ceil(keys.size() / overload_bucket_size);
            if (cur_buckets >= remaining) {
                cur_buckets = remaining;
            } else if ((double)(keys.size() - k * cur_buckets) / (remaining - cur_buckets) > overload_bucket_size) {
                cur_buckets = remaining;
            }

            Builder builder(k, threshold_size, thresholds, errors, keys, i, cur_buckets, offset, &spots);
            Consensus consensus(builder);

            layers.emplace_back(cur_buckets, std::move(consensus));

            offset += cur_buckets;
        }
        layers.shrink_to_fit();

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
            if (h >= actual_spots.size()) {
                actual_spots.resize(h + 1);
            }
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

