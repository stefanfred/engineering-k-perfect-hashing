#pragma once

#include <vector>
#include <limits>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <ranges>


double gammapdf(double x, double rate, double shape) {
    if (x < 0.000001) {
        return 0.0;
    }
    double expo = (shape - 1) * log(x) - rate * x + shape * log(rate) + shape - shape * log(shape);
    return exp(expo);
}

struct PrecomputedIntegral {
    static constexpr size_t granularity = 10000;
    double vals[granularity];

    PrecomputedIntegral(double rate, double shape) {
        double sum = 0;
        vals[0]=0;
        for (size_t i = 1; i < granularity; ++i) {
            double x = (double(i) + 0.5) / (granularity-1);
            sum += gammapdf(x, rate, shape) / x;
            vals[i] = sum / (granularity-1);
        }
    }

    double operator()(double x) {
        double index = x * double(granularity - 1);
        size_t indexN = size_t(index);
        if(indexN==granularity-1) {
            return vals[indexN];
        }
        double inter = index - double(indexN);
        return std::lerp(vals[indexN], vals[indexN + 1], inter);
    }
};

auto compute_thresholds(uint64_t _k, double bucket_size, uint64_t n_thresholds) {
    double shape = _k + 1;
    double rate = bucket_size;
    PrecomputedIntegral prec(rate, shape);
    std::vector<double> res(n_thresholds);
    double delta = 0.5;
    res[n_thresholds - 1] = 1.0; // last must be 1.0
    res[n_thresholds - 2] = 1.0 - delta; // initial guess for previous
    while (true) { //binary search
        bool fail = false;
        for (int64_t i = int64_t(n_thresholds) - 3; i >= 0; --i) {
            res[i] = res[i + 1] -
                     (res[i + 1] * (prec(res[i + 2]) - prec(res[i + 1]))) / gammapdf(res[i + 1], rate, shape);
            if (res[i] >= 1.0) {
                delta /= 2.0;
                res[n_thresholds - 2] -= delta;
                fail = true;
                break;
            }
            if (res[i] < 0.0) {
                delta /= 2.0;
                res[n_thresholds - 2] += delta;
                fail = true;
                break;
            }
        }
        if (fail) {
            continue;
        }
        if (res[0] > 0.001) {
            delta /= 2.0;
            res[n_thresholds - 2] -= delta;
            continue;
        }
        res[0] = 0.0;


        return res;
    }
}