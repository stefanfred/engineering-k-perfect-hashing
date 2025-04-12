#pragma once

#include <vector>
#include <ranges>

#ifdef STATS
#include <map>
#endif

#include <ips2ra.hpp>
#include <bytehamster/util/Function.h>

#include <hash128.hpp>

#include "../../extlib/simple-ribbon/extlib/ribbon/DySECT/include/bucket.h"

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

    template<size_t k, typename BucketFunction, typename Encoding>
    class HashDisplace {
    private:
        uint64_t nbuckets;
        uint64_t nbins;
        BucketFunction bucketFunction;
        Encoding seeds;
        std::array<uint32_t, 10> fastSeeds;
    public:
        uint64_t operator()(const std::string &item) const {
            return operator()(Hash128(item));
        }

        uint64_t operator()(Hash128 item) const {
            uint64_t bucket = bucketFunction(item.hi);
            uint64_t seed = bucket < fastSeeds.size() ? fastSeeds[bucket] : seeds[bucket];
            return bytehamster::util::fastrange64(item.lo * ~seed, nbins);
        }

        [[nodiscard]] size_t count_bits() const {
            return 8 * sizeof(*this)
                   + fastSeeds.size() * sizeof(uint32_t) * 8
                   + (seeds.count_bits() - 8 * sizeof(seeds));
        }

        HashDisplace()
                : nbuckets(0), nbins(0), bucketFunction(0, 1.0) {
        }

        HashDisplace(const std::vector<std::string> &items, uint64_t bucket_size, double load_factor = 1.0)
                : nbuckets((items.size() + bucket_size - 1) / bucket_size),
                  nbins(ceil(items.size() / load_factor / k)),
                  bucketFunction(nbuckets, load_factor) {
            assert(bucket_size > 0);
            assert(0.0 <= load_factor && load_factor <= 1.0);

            if (items.empty()) {
                return;
            }

            std::vector<std::pair<uint64_t, uint64_t>> sorted_items;
            sorted_items.reserve(items.size());
            for (const std::string &key: items) {
                Hash128 hash(key);
                sorted_items.emplace_back(bucketFunction(hash.hi), hash.lo);
            }
            ips2ra::sort(sorted_items.begin(), sorted_items.end(),
                         [](const std::pair<uint64_t, uint64_t> &key) { return key.first; });

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
                    while (it != sorted_items.end() && it->first == i) {
                        ++it;
                    }
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
                        uint64_t hash = bytehamster::util::fastrange64(bucket[j].second * ~seed, nbins);
                        if (counts[hash] == k) {
                            break;
                        }
                        counts[hash]++;
                    }
                    if (j == bucket.size()) {
                        break;
                    }
                    while (j > 0) {
                        uint64_t hash = bytehamster::util::fastrange64(bucket[--j].second * ~seed, nbins);
                        counts[hash]--;
                    }
                }
                seeds_vec[i] = seed;
#ifdef STATS
                *stats_it++ += seed;
#endif
            }
            seeds = Encoding(seeds_vec);
            for (size_t i = 0; i < fastSeeds.size(); ++i) {
                fastSeeds[i] = seeds_vec[i];
            }
        }
    };
}
