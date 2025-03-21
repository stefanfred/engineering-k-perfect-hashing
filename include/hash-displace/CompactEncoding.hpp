#pragma once

#include <bit>
#include <vector>
#include <cstring>

namespace kphf::HashDisplace {
class CompactEncoding {
        uint64_t bits_per_bucket = 0;
        std::vector<uint8_t> encoding;
    public:
        CompactEncoding() = default;

        explicit CompactEncoding(const std::vector<uint64_t> &seeds) {
            bits_per_bucket = std::bit_width(*std::ranges::max_element(seeds));
            encoding.resize((bits_per_bucket * seeds.size() + 7) / 8 + 7);
            encoding.shrink_to_fit();
            for (size_t i = 0; i < seeds.size(); i++) {
                size_t off = i * bits_per_bucket;
                uint64_t word;
                memcpy(&word, encoding.data() + off / 8, 8);
                word |= seeds[i] << (off % 8);
                memcpy(encoding.data() + off / 8, &word, 8);
            }
        }

        uint64_t operator[](uint64_t i) const {
            size_t off = i * bits_per_bucket;
            uint64_t word;
            memcpy(&word, encoding.data() + off / 8u, 8);
            return (word >> (off % 8u)) & ((1ul << bits_per_bucket) - 1u);
        }

        [[nodiscard]] size_t count_bits() const {
            return 8 * sizeof(*this)
                   + encoding.capacity() * 8;
        }

        static std::string name() { return "compact"; }
};
}
