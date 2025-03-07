#pragma once

#include <cstdint>

namespace kphf {

uint64_t remix(uint64_t z) {
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}

uint64_t rescale(uint64_t z, uint64_t c) {
	return ((unsigned __int128) z * c) >> 64;
}

}
