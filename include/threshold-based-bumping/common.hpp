#pragma once

namespace kphf::ThresholdBasedBumping {
	struct Key {
		uint64_t bucket, fingerprint;
		Hash128 hash;
	};

	template<typename R>
	void sort_keys(R &&r) {
		std::sort(r.begin(), r.end(), [](const Key &a, const Key &b) {
			return a.bucket < b.bucket || (a.bucket == b.bucket && a.fingerprint < b.fingerprint);
		});
		// TODO: Use ips2ra again, which does not natively support uint128
		//ips2ra::sort(r.begin(), r.end(),
		//  [](const Key &key) -> unsigned __int128 {
		//	return (unsigned __int128)key.bucket << 64 | key.fingerprint;
		//});
	}

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
