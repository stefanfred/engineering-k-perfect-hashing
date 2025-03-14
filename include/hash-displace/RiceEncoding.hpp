#pragma once

#include <bit>
#include <vector>
#include <ranges>
#ifdef STATS
#include <map>
#endif

#include <ips2ra.hpp>
#include <bytehamster/util/Function.h>
#include <sux/bits/SimpleSelectHalf.hpp>

#include <hash128.hpp>

namespace kphf::HashDisplace {

class RiceEncoding {
private:
	uint64_t length;
	std::vector<uint64_t> unary;
	mutable sux::bits::SimpleSelectHalf<> unary_select;
	std::vector<uint8_t> fixed;

public:
	RiceEncoding() = default;

	RiceEncoding(const std::vector<uint64_t> &seeds) {
		{
			double avg =
			  (double) std::accumulate(seeds.begin(), seeds.end(), uint64_t(0))
					   / seeds.size();
			double p = 1 / (avg+1);
			if (p == 1) length = 0;
			else {
				double phi = (sqrt(5) + 1) / 2;
				length =
				  std::max(uint64_t(ceil(log2(- log(phi) / log1p(-p)))),
						   uint64_t(0));
			}
		}
		fixed.resize((length * seeds.size() + 7) / 8 + 7);
		fixed.shrink_to_fit();
		uint64_t pos = 0;
		for (size_t i = 0; i < seeds.size(); i++) {
			{
				uint64_t fixed_part = seeds[i] & ((uint64_t(1) << length) - 1);
				size_t fixed_off = i * length;
				uint64_t word;
				memcpy(&word, fixed.data() + fixed_off / 8, 8);
				word |= fixed_part << (fixed_off % 8);
				memcpy(fixed.data() + fixed_off / 8, &word, 8);
			}
			{
				uint64_t unary_part = seeds[i] >> length;
				pos += unary_part;
				unary.resize(pos/64 + 1);
				unary.back() |= uint64_t(1) << (pos % 64);
				pos++;
			}
		}
		unary.shrink_to_fit();
		unary_select = sux::bits::SimpleSelectHalf(unary.data(), pos);
	}

	uint64_t operator[](uint64_t i) const {
		uint64_t res;
		{
			size_t off = i * length;
			memcpy(&res, fixed.data() + off / 8, 8);
			res >>= (off%8);
			res &= (uint64_t(1) << length) - 1;
		}
		{
			size_t unary_start = i == 0 ? 0 : unary_select.select(i-1)+1;
			size_t pos = unary_start/64;
			uint64_t window =
			  unary[pos] & ~((uint64_t(1) << (unary_start%64)) - 1);
			while (window == 0) window = unary[++pos];
			uint64_t unary = 64 * pos + __builtin_ctzll(window) - unary_start;
			res |= unary << length;
		}
		return res;
	}

	size_t count_bits() const {
		return 8 * sizeof(*this)
		  + unary.capacity() * 64
		  + unary_select.bitCount() - sizeof(unary_select) * 8
		  + fixed.capacity() * 8;
	}

	static std::string name() { return "rice"; }
};

}
