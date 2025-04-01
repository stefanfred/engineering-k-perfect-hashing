#pragma once

#include <cstdint>
#include <bytehamster/util/MurmurHash64.h>
#include <sux/support/SpookyV2.hpp>

struct Hash128 {
    uint64_t hi, lo;

    explicit Hash128() : hi(0), lo(0) {
    }

    explicit Hash128(const std::string &string) : hi(0), lo(0) {
        SpookyHash::Short128(string.data(), string.length(), &hi, &lo);
    }
};
