#pragma once

namespace kphf::ThresholdBasedBumping {
	struct Key {
		uint64_t bucket, fingerprint;
		Hash128 hash;
	};

	template<typename R>
	void sort_buckets(R &&r) {
		ips2ra::sort(r.begin(), r.end(),
		  [](const Key &key) -> uint64_t { return key.bucket; });
	}

	template<typename R>
	void sort_fingerprints(R &&r) {
		ips2ra::sort(r.begin(), r.end(),
		  [](const Key &key) -> uint64_t { return key.fingerprint; });
	}

	inline uint64_t double_to_u64(double x) {
		return x == 1.0 ? ~0ul : static_cast<uint64_t>(ldexp(x, 64));
	}

}
