#pragma once

#include <ips2ra.hpp>
#include <bytehamster/util/Function.h>

namespace kphf::HashDisplace {
class UniformBucketFunction {
    private:
        size_t nbuckets;
    public:
        explicit UniformBucketFunction(size_t nbuckets, double loadFactor) : nbuckets(nbuckets) {
            (void) loadFactor;
        }

        uint64_t operator()(uint64_t hash) const {
            return bytehamster::util::fastrange64(hash, nbuckets);
        }

        static std::string name() {
            return "uniform";
        }
};
}
