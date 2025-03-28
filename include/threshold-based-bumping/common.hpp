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

}
