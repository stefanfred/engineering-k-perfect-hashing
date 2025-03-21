#pragma once

#include <cstdint>
#include <bytehamster/util/MurmurHash64.h>

struct Hash128 {
    uint64_t hi, lo;

    explicit Hash128(const std::string &string)
        : hi(bytehamster::util::MurmurHash64(string)),
          lo(bytehamster::util::MurmurHash64(string + "abc")) {
        // TODO: This is inefficient, use 128-bit hash function
    }
};
