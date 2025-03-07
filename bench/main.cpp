#include <map>
#include <ranges>
#include <string>

#include "bench.hpp"
#include "ThresholdBasedBumpingConsensusContender.hpp"
#include "HashDisplaceContender.hpp"
#include "PaCHashContender.hpp"

const std::map<std::string, void (*)(Benchmarks &)> contenders = {
	{"ThresholdBasedBumpingConsensus",
		&ThresholdBasedBumpingConsensusContender::benchmark},
	{"HashDisplace",
		&HashDisplaceContender::benchmark},
	{"PaCHash",
		&PaCHashContender::benchmark},
};

int main(int argc, char **argv) {
	Benchmarks bench;
	if (argc == 1) {
		for (auto f: contenders | std::views::values) (*f)(bench);
	} else {
		bool ok = true;
		for (int i = 1; i < argc; i++) {
			auto it = contenders.find(argv[i]);
			if (it == contenders.end()) {
				std::cerr << "Unknown contender: " << argv[i] << "\n";
				ok = false;
			} else (*it->second)(bench);
		}
		if (!ok) {
			std::cerr << "Available contenders:\n";
			for (auto name: contenders | std::views::keys) {
				std::cerr << "  " << name << "\n";
			}
			return 127;
		}
	}
	bench.run(10);
	return 0;
}
