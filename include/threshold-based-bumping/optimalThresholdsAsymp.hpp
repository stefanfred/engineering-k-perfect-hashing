#pragma once
#include "optimalThresholds.hpp"

namespace kphf::ThresholdBasedBumping {
auto compute_thresholds_asymp(uint64_t _k, double bucket_size, uint64_t n_thresholds) {
    double shape = 0.5 * double(_k + 1);
    double rate = 0.5 * bucket_size;
    constexpr size_t samples = 10000;
    std::array<double, samples> vals{};
    vals[0] = 0;
    for (size_t i = 1; i < samples; ++i) {
        double x = (double(i) + 0.5) / (samples - 1);
        vals[i] = vals[i - 1] + gammapdf(x, rate, shape);
    }
    for (size_t i = 1; i < samples; ++i) {
        vals[i] /= vals[samples - 1];
    }


    std::vector<double> res;
    res.push_back(0.0);
    size_t arrayIndex = 0;
    for (size_t i = 1; i < n_thresholds - 1; ++i) {
        double area = double(i) / double(n_thresholds - 1);
        while (vals[arrayIndex] < area) {
            arrayIndex++;
        }
        res.push_back(double(arrayIndex) / double(samples));
    }
    res.push_back(1.0);

    return res;
}
}
