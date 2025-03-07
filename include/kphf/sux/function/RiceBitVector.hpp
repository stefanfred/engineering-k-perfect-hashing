/*
 * Sux: Succinct data structures
 *
 * Copyright (C) 2019-2020 Emmanuel Esposito and Sebastiano Vigna
 *
 *  This library is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * Under Section 7 of GPL version 3, you are granted additional permissions
 * described in the GCC Runtime Library Exception, version 3.1, as published by
 * the Free Software Foundation.
 *
 * You should have received a copy of the GNU General Public License and a copy of
 * the GCC Runtime Library Exception along with this program; see the files
 * COPYING3 and COPYING.RUNTIME respectively.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "../support/common.hpp"
#include "../util/Vector.hpp"
#include <cstdint>
#include <cstdio>
#include <iostream>

namespace sux::function {

using namespace std;
using namespace sux;

/** Storage for Golomb-Rice codes of a RecSplit bucket.
 *
 * This class exists solely to implement RecSplit.
 * @tparam AT a type of memory allocation out of util::AllocType.
 */

template <util::AllocType AT = util::AllocType::MALLOC> class RiceBitVector {

  public:
	class Builder {
		util::Vector<uint64_t, AT> data;
		size_t bit_count = 0;

	  public:
		Builder() : Builder(16) {}

		Builder(const size_t alloc_words) : data(alloc_words) {}

		void appendFixed(const uint64_t v, const int log2golomb) {
			const uint64_t lower_bits = v & ((uint64_t(1) << log2golomb) - 1);

			data.resize((((bit_count + log2golomb + 7) / 8) + 7 + 7) / 8);

			appendFixedInner(lower_bits, log2golomb);
		}

		void appendUnaryAll(const std::vector<uint32_t> unary) {
			size_t bit_inc = 0;
			for (const auto u : unary) {
				bit_inc += u + 1;
			}

			data.resize((((bit_count + bit_inc + 7) / 8) + 7 + 7) / 8);

			for (const auto u : unary) appendUnaryInner(u);
		}

		void appendGolomb(const uint64_t v, const double p) {
			const uint64_t m = golomb_modulus(p);
			const uint32_t b = bit_width(m)-1;
			const uint64_t threshold = (uint64_t(2) << b) - m;

			const uint32_t unary = v / m;
			uint64_t fixed = v % m;
			const uint32_t fixed_len = b + (fixed >= threshold);
			if (fixed & (uint64_t(1) << b)) fixed += threshold;

			size_t bit_inc = unary + 1 + fixed_len;

			data.resize((((bit_count + bit_inc + 7) / 8) + 7 + 7) / 8);

			appendUnaryInner(unary);
			appendFixedInner(fixed, fixed_len);
		}

		uint64_t getBits() { return bit_count; }

		RiceBitVector<AT> build() {
			data.trimToFit();
			return RiceBitVector(std::move(data));
		}

	  private:
		inline void appendFixedInner(const uint64_t v, const int length) {
			int used_bits = bit_count & 63;
			uint64_t *append_ptr = &data + bit_count / 64;
			uint64_t cur_word = *append_ptr;

			cur_word |= v << used_bits;
			if (used_bits + length > 64) {
				*(append_ptr++) = cur_word;
				cur_word = v >> (64 - used_bits);
				used_bits += length - 64;
			}
			*append_ptr = cur_word;
			bit_count += length;
		}

		inline void appendUnaryInner(const uint32_t u) {
			bit_count += u;
			uint64_t *append_ptr = &data + bit_count / 64;
			*append_ptr |= uint64_t(1) << (bit_count & 63);
			++bit_count;
		}
	};

  private:
	util::Vector<uint64_t, AT> data;

	friend std::ostream &operator<<(std::ostream &os, const RiceBitVector<AT> &rbv) {
		os << rbv.data;
		return os;
	}

	friend std::istream &operator>>(std::istream &is, RiceBitVector<AT> &rbv) {
		is >> rbv.data;
		return is;
	}

	static inline uint64_t golomb_modulus(const double p) {
		return round(-log(2) / log1p(-p));
	}

  public:
	RiceBitVector() {}
	RiceBitVector(util::Vector<uint64_t, AT> data) : data(std::move(data)) {}

	size_t getBits() const { return data.size() * sizeof(uint64_t) * 8; }

	class Reader {
		size_t curr_fixed_offset = 0;
		uint64_t curr_window_unary = 0;
		uint64_t *curr_ptr_unary;
		int valid_lower_bits_unary = 0;
		util::Vector<uint64_t, AT> &data;

	  public:
		Reader(util::Vector<uint64_t, AT> &data) : data(data) {}

		uint64_t readNext(const int log2golomb) {
			uint64_t result = readUnary() << log2golomb;

			uint64_t fixed;
			memcpy(&fixed, (uint8_t *)&data + curr_fixed_offset / 8, 8);
			result |= (fixed >> curr_fixed_offset % 8) & ((uint64_t(1) << log2golomb) - 1);
			curr_fixed_offset += log2golomb;
			return result;
		}

		uint64_t readGolomb(const double p) {
			const uint32_t m = golomb_modulus(p);
			const int b = bit_width(m)-1;
			const uint64_t threshold = (uint64_t(2) << b) - m;

			uint64_t result = readUnary() * m;

			uint64_t fixed = curr_window_unary;
			if (valid_lower_bits_unary <= b) {
				int missing = b - valid_lower_bits_unary;
				curr_window_unary = *(curr_ptr_unary++);
				fixed |= curr_window_unary << valid_lower_bits_unary;
				curr_window_unary >>= missing;
				valid_lower_bits_unary = 64 - missing;
			} else {
				curr_window_unary >>= b;
				valid_lower_bits_unary -= b;
			}
			fixed &= (uint64_t(1) << b) - 1;

			if (fixed >= threshold) {
				if (curr_window_unary & 1) {
					fixed += (uint64_t(1) << b) - threshold;
				}
				curr_window_unary >>= 1;
				valid_lower_bits_unary--;
			}
			return result + fixed;
		}

		void skipSubtree(const size_t nodes, const size_t fixed_len) {
			assert(nodes > 0);
			size_t missing = nodes, cnt;
			while ((cnt = nu(curr_window_unary)) < missing) {
				curr_window_unary = *(curr_ptr_unary++);
				missing -= cnt;
				valid_lower_bits_unary = 64;
			}
			cnt = select64(curr_window_unary, missing - 1);
			curr_window_unary >>= cnt;
			curr_window_unary >>= 1;
			valid_lower_bits_unary -= cnt + 1;

			curr_fixed_offset += fixed_len;
		}

		void readReset(const size_t bit_pos, const size_t unary_offset) {
			// assert(bit_pos < bit_count);
			curr_fixed_offset = bit_pos;
			size_t unary_pos = bit_pos + unary_offset;
			curr_ptr_unary = &data + unary_pos / 64;
			curr_window_unary = *(curr_ptr_unary++) >> (unary_pos & 63);
			valid_lower_bits_unary = 64 - (unary_pos & 63);
		}

	  private:
		inline uint64_t readUnary() {
			uint64_t result = 0;

			if (curr_window_unary == 0) {
				result += valid_lower_bits_unary;
				curr_window_unary = *(curr_ptr_unary++);
				valid_lower_bits_unary = 64;
				while (__builtin_expect(curr_window_unary == 0, 0)) {
					result += 64;
					curr_window_unary = *(curr_ptr_unary++);
				}
			}

			const size_t pos = rho(curr_window_unary);

			curr_window_unary >>= pos;
			curr_window_unary >>= 1;
			valid_lower_bits_unary -= pos + 1;

			result += pos;
			return result;
		}
	};

	Reader reader() { return Reader(data); }
};

} // namespace sux::function
