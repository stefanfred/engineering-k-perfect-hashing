#pragma once

#include <ips2ra.hpp>
#include <bytehamster/util/Function.h>

#include <hash128.hpp>

namespace kphf::HashDisplace {

class UniformBucketFunction {
public:
    class Instance {
    private:
        uint64_t nbuckets;

        Instance(uint64_t nbuckets): nbuckets(nbuckets) {}

    public:
        uint64_t operator()(uint64_t z) const {
            return bytehamster::util::fastrange64(z, nbuckets);
        }

        size_t count_bits() const {
            return 8 * sizeof(*this);
        }

        friend class UniformBucketFunction;
    };

    UniformBucketFunction(uint64_t k) {
        (void) k;
    }

    Instance
      operator()(uint64_t n, uint64_t nbuckets, double load_factor) const {
        (void) n;
        (void) load_factor;
        return Instance(nbuckets);
    }

    static std::string name() { return "uniform"; }
};

}
