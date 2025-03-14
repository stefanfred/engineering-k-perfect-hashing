#include <tlx/cmdline_parser.hpp>

#include "bench.hpp"
#include "ThresholdBasedBumpingConsensusContender.hpp"
#include "HashDisplaceContender.hpp"
#include "PaCHashContender.hpp"

int main(int argc, char **argv) {
	size_t k = 8;

    tlx::CmdlineParser cmd;
	cmd.add_bytes('n', "numKeys", Benchmark::numKeys, "Number of objects");
	cmd.add_bytes('q', "numQueries", Benchmark::numQueries, "Number of queries to perform");
	cmd.add_bytes('k', "k", k, "Maximum number of keys per bucket");
	if (!cmd.process(argc, argv)) {
		return 1;
	}

	PaCHashContender::benchmark(k);
	ThresholdBasedBumpingConsensusContender::benchmark(k);
	HashDisplaceContender::benchmark(k);
	// TODO: RecSplit, CMPH
	return 0;
}
