#include <map>
#include <ranges>
#include <string>
#include <tlx/cmdline_parser.hpp>

#include "bench.hpp"
#include "ThresholdBasedBumpingConsensusContender.hpp"
#include "HashDisplaceContender.hpp"
#include "PaCHashContender.hpp"

// TODO: RecSplit?
// TODO: Add CMPH library, which also has a k-perfect (not minimal) hash function
const std::map<std::string, void (*)(Benchmarks &)> contenders = {
	{"ThresholdBasedBumpingConsensus",
		&ThresholdBasedBumpingConsensusContender::benchmark},
	{"HashDisplace",
		&HashDisplaceContender::benchmark},
	{"PaCHash",
		&PaCHashContender::benchmark},
};

int main(int argc, char **argv) {
	// TODO: Use these
    size_t numKeys = 1e6;
	size_t numQueries = 1e6;
	size_t k = 8;

    tlx::CmdlineParser cmd;
	cmd.add_bytes('n', "numKeys", numKeys, "Number of objects");
	cmd.add_bytes('q', "numQueries", numQueries, "Number of queries to perform");
	cmd.add_bytes('k', "k", k, "Maximum number of keys per bucket");
	if (!cmd.process(argc, argv)) {
		return 1;
	}

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
