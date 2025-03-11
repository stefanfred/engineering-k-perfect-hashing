#pragma once

#include <bit>
#include <sux/bits/SimpleSelectZeroHalf.hpp>

namespace kphf {
namespace PaCHash {

class EliasFano {
private:
	using SimpleSelectZeroHalf = sux::bits::SimpleSelectZeroHalf<>;

	int l;
	std::vector<uint64_t> upper;
	mutable SimpleSelectZeroHalf upper_select;
	std::vector<uint8_t> lower;

	EliasFano(int l, std::vector<uint64_t> &&upper,
	  SimpleSelectZeroHalf &&upper_select, std::vector<uint8_t> &&lower):
	  l(l), upper(std::move(upper)), upper_select(std::move(upper_select)),
	  lower(std::move(lower)) {}

public:
	EliasFano(const std::vector<uint64_t> &items):
	  EliasFano(std::move(build(items))) {}

	std::pair<size_t, size_t> search(uint64_t v) const {
		uint64_t high = v >> l;
		uint64_t low = v & ((1 << l) - 1);
		uint64_t pos = upper_select.selectZero(high);
		pos--;
		while (get_upper_bit(pos) && get_lower(pos - high) > low) pos--;
		size_t end = pos - high + 1;
		while (get_upper_bit(pos) && get_lower(pos - high) == low) pos--;
		size_t start = pos - high + 1;
		return {start,end};
	}

	size_t count_bits() const {
		return sizeof(*this)
			+ upper.capacity() * sizeof(upper[0]) * 8
			+ upper_select.bitCount() - sizeof(upper_select) * 8
			+ lower.capacity() * sizeof(lower[0]) * 8
		;
	}

private:
	bool get_upper_bit(uint64_t pos) const {
		return (upper[pos/64] >> pos%64) & 1;
	}

	uint64_t get_lower(size_t idx) const {
		size_t off = idx * l;
		uint64_t word;
		memcpy(&word, lower.data() + off/8, 8);
		return (word >> off%8) & ((1 << l) - 1);
	}

	static EliasFano build(const std::vector<uint64_t> &items) {
		uint64_t u = items.back();
		uint64_t n = items.size();
		int l = std::bit_width((u - 1) / n);

		int upper_length = n + (u >> l) + 1;
		std::vector<uint64_t> upper((upper_length + 63) / 64);
		upper.shrink_to_fit();
		for (uint64_t i = 0; i < n; i++) {
			uint64_t j = i + (items[i] >> l);
			upper[j/64] |= uint64_t(1) << (j%64);
		}
		SimpleSelectZeroHalf upper_select(upper.data(), upper_length);

		std::vector<uint8_t> lower((l * n + 7) / 8 + 7);
		lower.shrink_to_fit();
		for (uint64_t i = 0; i < n; i++) {
			uint64_t off = i * l;
			uint64_t word;
			memcpy(&word, lower.data() + off/8, 8);
			word |= (items[i] & ((1 << l) - 1)) << off%8;
			memcpy(lower.data() + off/8, &word, 8);
		}

		return EliasFano(l, std::move(upper), std::move(upper_select),
		  std::move(lower));
	}
};

}
}
