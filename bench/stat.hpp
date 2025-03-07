#pragma once

#include <numeric>
#include <cmath>
#include <iostream>

namespace stat {

struct Statistics {
	struct Builder {
		double minval =  std::numeric_limits<double>::infinity();
		double maxval = -std::numeric_limits<double>::infinity();
		double sum = 0;
		double sqsum = 0;
		int count = 0;

		Builder() {}

		void add(double val) {
			minval = std::min(minval, val);
			maxval = std::max(maxval, val);
			sum += val;
			sqsum += val * val;
			count++;
		}

		Statistics build() {
			Statistics res;
			res.min = minval;
			res.max = maxval;
			res.mean = sum / count;
			res.stddev = sqrt((sqsum - sum * sum / count) / (count-1));
			return res;
		}
	};

	double min, max, mean, stddev;
};

std::ostream &operator<<(std::ostream &o, const Statistics &stat) {
	return o << stat.min << " ... " << stat.max
	         << " (μ = " << stat.mean << ", σ = " << stat.stddev << ")";
}

}
