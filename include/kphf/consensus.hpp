
#pragma once

#include <vector>
#include <optional>
#include <utility>
#include <cstdint>
#include <cstring>

uint64_t backtracks = 0;

class Consensus {
private:
	std::vector<uint8_t> data;

public:
	template<typename T>
	Consensus(T &&builder) {
		std::uint64_t size = 0;
		auto consensus = [&]() -> std::uint64_t {
			return get(size, 0).first;
		};
		auto push = [&](std::uint64_t value, std::uint64_t length) -> void {
			data.resize((64 + size + length + 7) / 8 + 8);
			std::uint64_t off = size + 64;
			std::uint64_t v;
			std::memcpy(&v, data.data() + off/8, 8);
			v &= ~(((uint64_t(1) << length) - 1) << off%8);
			v |= (value << off%8);
			std::memcpy(data.data() + off/8, &v, 8);
			size += length;
		};
		auto pop = [&](std::uint64_t length) -> std::pair<uint64_t, uint64_t> {
			size -= length;
			return get(size, length);
		};
		data.resize(16);

		uint64_t idx = 0;
		for (;;) {
			std::uint64_t cons = consensus();
			std::optional<std::uint64_t> s = builder.advance(cons);
			if (!s) break;
			idx++;
			std::uint64_t length = *s;

			std::optional<std::uint64_t> result;
			result = builder.find_first(cons);
			while (!result) {
				length = builder.backtrack();
				backtracks++;
				if (--idx == 0) break;
				auto [cons, prev] = pop(length);
				result = builder.find_next(cons, prev);
			}

			if (result) push(*result, length);
			else {
				std::uint64_t v;
				std::memcpy(&v, data.data(), 8);
				v++;
				std::memcpy(data.data(), &v, 8);
			}
		}
		data.resize((64 + size + 7) / 8 + 7);
		data.shrink_to_fit();
		
	}

	std::pair<std::uint64_t, std::uint64_t>
	  get(std::uint64_t off, std::uint64_t size) const {
		unsigned __int128 v;
		std::memcpy(&v, data.data() + off/8, 16);
		v >>= off%8;
		uint64_t consensus = v & uint64_t(-1);
		uint64_t value = (v >> 64) & ((uint64_t(1) << size) - 1);

		return {consensus, value};
	}

	std::size_t count_bits() const {
		return sizeof(*this) * 8 + data.capacity() * 8;
	}
};

#if 0
	auto advance(std::uint64_t consensus) -> std::optional<std::uint64_t>;
	auto backtrack() -> std::uint64_t;
	auto find_first(std::uint64_t consensus) -> std::optional<std::uint64_t>;
	auto find_next(std::uint64_t consensus, std::uint64_t prev) -> std::optional<std::uint64_t>;
#endif
