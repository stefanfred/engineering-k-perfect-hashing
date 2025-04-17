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

#include "optimalThresholdsAsymp.hpp"

namespace kphf::ThresholdBasedBumping {

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

[[nodiscard]] inline uint64_t getBucket(Hash128 hash, size_t cur_buckets) {
    return bytehamster::util::fastrange64(hash.hi, cur_buckets);
}

inline void sort_buckets(std::vector<Hash128> &r, size_t cur_buckets) {
    ips2ra::sort(r.begin(), r.end(),
                 [&](const Hash128 &key) -> uint64_t { return getBucket(key, cur_buckets); });
}

inline void sort_fingerprints(std::vector<Hash128> &r) {
    ips2ra::sort(r.begin(), r.end(),
                 [](const Hash128 &key) -> uint64_t { return key.lo; });
}

template<typename RandomIt>
RandomIt partition(RandomIt first, RandomIt last) {
    auto pivot = *(last - 1); // choose the last element as pivot
    auto store = first;

    for (auto it = first; it < last - 1; ++it) {
        if (it->lo < pivot.lo) {
            std::iter_swap(it, store);
            ++store;
        }
    }
    std::iter_swap(store, last - 1); // move pivot to final place
    return store;
}

template<typename RandomIt>
RandomIt quickSelect(RandomIt first, RandomIt last, size_t k) {
    while (first < last) {
        auto pivotIt = partition(first, last);
        size_t index = pivotIt - first;

        if (index == k) {
            return pivotIt;
        } else if (index < k) {
            first = pivotIt + 1;
            k -= index + 1;
        } else {
            last = pivotIt;
        }
    }
    return last;
}

template<uint64_t K, int THRESHOLD_SIZE, bool packing>
class ThresholdBasedBumping {
    static_assert(THRESHOLD_SIZE > 1);
    static_assert(K > 0);
private:
    static constexpr uint64_t _k = K;
    static constexpr uint64_t threshold_size = THRESHOLD_SIZE;
    static constexpr uint64_t n_regions = uint64_t(uint64_t(1) << threshold_size);
    static constexpr uint64_t n_thresholds = n_regions - 1;

    std::array<uint64_t, n_thresholds> avail_thresholds;

    uint64_t n;
    std::vector<uint64_t> nbuckets;
    std::vector<uint64_t> layerSpotsBound;
    std::vector<uint8_t> thresholds;
    SimpleRibbon<1> retrieval;
    fips::FiPS<> phf;
    mutable sux::bits::EliasFano<> gaps;

public:
    uint64_t operator()(const std::string &key) const {
        return operator()(Hash128(key));
    }

    uint64_t operator()(Hash128 hash) const {
        uint64_t offset = 0;
        for (uint64_t i = 0; i < nbuckets.size(); i++) {
            uint64_t cur_buckets = nbuckets[i];
            uint64_t keyBucket = getBucket(hash, cur_buckets);
            keyBucket += offset;
            uint64_t tidx;
            uint64_t off = keyBucket * threshold_size;
            memcpy(&tidx, thresholds.data() + off/8, 8);
            tidx >>= off%8;
            tidx &= (uint64_t(1) << threshold_size) - 1;

            if constexpr (packing) {
                if (tidx != 0 && hash.lo < avail_thresholds[tidx-1]) {
                    return keyBucket;
                }
                if (tidx == n_thresholds || hash.lo < avail_thresholds[tidx]) {
                    if (retrieval.retrieve(hash.lo)) {
                        return keyBucket;
                    }
                }
            } else {
                if (tidx != 0 && hash.lo < avail_thresholds[tidx-1]) {
                    return keyBucket;
                }
            }

            if(hash.hi < layerSpotsBound[i]) {
                return gaps.select(phf(hash.hi ^ hash.lo));
            }

            hash = Hash128(bytehamster::util::remix(hash.hi),bytehamster::util::remix(hash.lo));

            offset += cur_buckets;
        }

        return gaps.select(phf(hash.hi ^ hash.lo));
    }

    size_t count_bits() const {
        return sizeof(*this) * 8
            - sizeof(avail_thresholds) * 8 // Lookup table independent of keys, could actually be static constexpr
            + nbuckets.capacity() * 64
            + layerSpotsBound.capacity() * 64
            + thresholds.capacity() * 8
            + (packing ? retrieval.sizeBytes() : 0) * 8
            + phf.getBits() - sizeof(phf) * 8
            + gaps.bitCount() - sizeof(gaps) * 8
        ;
    }
    ThresholdBasedBumping()
        : avail_thresholds(), n(0), gaps(std::vector<uint64_t>(), 0) {
    }

    explicit ThresholdBasedBumping(const std::vector<std::string> &keys, double overload)
        : ThresholdBasedBumping(hashKeys(keys), overload) {
    }

    static std::vector<Hash128> hashKeys(const std::vector<std::string> &keys) {
        std::vector<Hash128> hashed_keys;
        hashed_keys.reserve(keys.size());
        for (const auto &key : keys) {
            hashed_keys.emplace_back(Hash128(key));
        }
        return hashed_keys;
    }

    explicit ThresholdBasedBumping(std::vector<Hash128> keys, double overload)
            : avail_thresholds(compute_thresholds<n_thresholds>(_k, _k * overload)),
              n(keys.size()),
              gaps(std::vector<uint64_t>(), 0) {
        double overload_bucket_size = _k * overload;
        uint64_t total_buckets = (n + _k - 1) / _k;
        thresholds.resize((total_buckets * threshold_size + 7) / 8 + 7);
        thresholds.shrink_to_fit();

#ifdef STATS
        total_keys += n;
#endif

        uint64_t offset = 0;
        std::vector<std::pair<uint64_t, uint8_t>> retrievalItems;
        std::vector<uint64_t> fallbackKeys;
        std::vector<uint64_t> spots;
        std::vector<Hash128> bumped;
        for (uint64_t layer = 0; offset != total_buckets; layer++) {
            size_t nThisLayer = layer == 0 ? keys.size() : bumped.size();
            uint64_t cur_buckets = std::min(total_buckets - offset, std::max(uint64_t(1),uint64_t(std::ceil(nThisLayer / overload_bucket_size))));
            nbuckets.push_back(cur_buckets);

            if (layer > 0) {
                bumped.swap(keys);
                bumped.clear();
            }

            uint64_t spotsBeforeThisLayer = spots.size();
            std::vector<Hash128> bumpOrFill;
            sort_buckets(keys, cur_buckets);
            auto next = keys.begin();
            for (uint64_t j = 0; j < cur_buckets; j++) {
                auto start = next;
                while (next != keys.end() && getBucket(*next, cur_buckets) == j) ++next;
                auto bucket = std::ranges::subrange(start, next);
                uint64_t tidx;
                if (bucket.size() <= _k) {
                    tidx = n_thresholds;
                } else {
                    Hash128 kthKey = *quickSelect(bucket.begin(), bucket.end(), _k);
                    tidx = std::ranges::upper_bound(avail_thresholds, kthKey.lo) - avail_thresholds.begin();
                }

                uint64_t in_bucket=0;
                auto currKey = bucket.begin();
                std::vector<Hash128> bumpOrFillPossibly;
                while (currKey != bucket.end()) {
                    if constexpr (packing) {
                        if (currKey->lo < avail_thresholds[tidx - 1]) {
                            in_bucket++;
                        } else {
                            if (tidx == n_thresholds || currKey->lo < avail_thresholds[tidx]) {
                                bumpOrFillPossibly.emplace_back(*currKey);
                            } else {
                                bumpOrFill.emplace_back(*currKey);
                            }
                        }
                    } else {
                        if (currKey->lo < avail_thresholds[tidx - 1]) {
                            in_bucket++;
                        } else {
                            bumpOrFill.emplace_back(*currKey);
                        }
                    }
                    currKey++;
                }

                if constexpr (packing) {
                    for (Hash128 h: bumpOrFillPossibly) {
                        if (in_bucket < _k) {
                            in_bucket++;
                            retrievalItems.emplace_back(h.lo, true);
                        } else {
                            bumpOrFill.push_back(h);
                            retrievalItems.emplace_back(h.lo, false);
                        }
                    }
                }

#ifdef STATS
                if (in_bucket == _k) filled_buckets++;
#endif

                uint64_t b = offset + j;
                uint64_t off = b * threshold_size;
                uint64_t word;
                memcpy(&word, thresholds.data() + off/8, 8);
                word |= tidx << off%8;
                memcpy(thresholds.data() + off/8, &word, 8);

                for (uint64_t c = 0; c < _k - in_bucket; c++) {
                    spots.push_back(offset + j);
                }
            }

            uint64_t spotsInThisLayer = spots.size() - spotsBeforeThisLayer;
            uint64_t relFill = double_to_u64(double(spotsInThisLayer) / double(bumpOrFill.size()));
            layerSpotsBound.push_back(relFill);
            for (Hash128 key : bumpOrFill) {
                if (key.hi < relFill) {
                    fallbackKeys.push_back(key.hi ^ key.lo);
                } else {
                    bumped.emplace_back(bytehamster::util::remix(key.hi),bytehamster::util::remix(key.lo));
                }
            }
            offset += cur_buckets;
        }
        nbuckets.shrink_to_fit();

#ifdef STATS
        bumped_keys += bumped.size();
#endif

        for (Hash128 key : bumped) {
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
        if constexpr (packing) {
            retrieval = SimpleRibbon<1>(retrievalItems);
        }
    }
};

}
