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
public:
    uint64_t operator()(const std::string &item) const {
        return operator()(Hash128(item));
    }

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

    PaCHash() : n_buckets(0) {
    }

    PaCHash(int k, double a, const std::vector<std::string> &keys) {
        n_buckets = std::ceil(a * keys.size() / k);
        std::vector<Hash128> keysHashed;
        keysHashed.reserve(keys.size());
        for (const std::string & key : keys) {
            Hash128 hash(key);
            hash.hi = bytehamster::util::fastrange64(hash.hi, n_buckets);
            keysHashed.push_back(hash);
        }
        build(k, a, keysHashed);
    }

private:
    void build(int k, double a, std::vector<Hash128> keys) {
        uint64_t n_bins = (keys.size() + k - 1) / k;
        n_buckets = std::ceil(a * keys.size() / k);

        ips2ra::sort(keys.begin(), keys.end(),[](const Hash128 &key) -> uint64_t { return key.hi; });

        std::vector<uint64_t> threshold(n_bins+1);
        threshold.front() = 0;
        for (size_t i = 1; i < n_bins; i++) {
            uint64_t b = keys[k*i].hi;
            uint64_t previousB = keys[k*i-1].hi;
            if (previousB < b - 1) { // Empty bin in between
                b--;
            } else if (previousB == b - 1) { // Store the smaller one
                size_t previousBucketSize = 0;
                while (k*i-1 >= previousBucketSize && keys[k*i-1 - previousBucketSize].hi == previousB) {
                    previousBucketSize++;
                }
                size_t nextBucketSize = 0;
                while (k*i + nextBucketSize < keys.size() && keys[k*i + nextBucketSize].hi == b) {
                    nextBucketSize++;
                }
                if (previousBucketSize < nextBucketSize) {
                    b--;
                }
            }
            threshold[i] = b + 1;
        }
        threshold.back() = n_buckets + 1;

        std::vector<std::pair<uint64_t, uint8_t>> ribbon_data;
        uint64_t offset = 0, next = 1;
        for (uint64_t i = 0; i < n_buckets; i++) {
            uint64_t start = next;
            while (threshold[next] == i+1) {
                next++;
            }
            uint64_t end = next;

            int bits = std::bit_width(end - start);
            uint64_t x = 0;
            for (Hash128 key; (key = keys[offset]).hi == i; offset++) {
                if (offset == (start+x) * k) {
                    x++;
                }
                for (int j = 0; j < bits; j++) {
                    ribbon_data.emplace_back(key.lo + uint64_t(j), (x >> j) & 1);
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

        ribbon = SimpleRibbon<1>(ribbon_data);
        std::cout<<"Num ribbon: "<<ribbon_data.size()<<std::endl;
        ef = EliasFano(threshold);
    }
};

}

