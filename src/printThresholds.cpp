#include "../include/threshold-based-bumping/optimalThresholdsAsymp.hpp"

#include <iostream>

int main(int argc, char **argv) {
    double gam = 1.2;
    size_t _k = 100;
    constexpr size_t n_thresholds = 32;

    auto res = kphf::ThresholdBasedBumping::compute_thresholds_normalized(_k, gam * _k, n_thresholds);
    auto resAsymp = kphf::ThresholdBasedBumping::compute_thresholds_asymp(_k, gam * _k, n_thresholds);
    for (size_t i = 0; i < n_thresholds; ++i) {
        std::cout << "RESULT m=exact i=" << (i + 1) << " t=" << res[i] << " k=" << _k << " g=" << gam << std::endl;
        std::cout << "RESULT m=asymp i=" << (i + 1) << " t=" << resAsymp[i] << " k=" << _k << " g=" << gam << std::endl;
    }
}