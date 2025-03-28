#pragma once

#include <cstddef>

template <template <size_t _> class Runner>
void dispatchK(size_t k, size_t N) {
    if (k == 2) {
        Runner<2>()(N);
    } else if (k == 4) {
        Runner<4>()(N);
    } else if (k == 8) {
        Runner<8>()(N);
    } else if (k == 10) {
        Runner<10>()(N);
    } else if (k == 100) {
        Runner<100>()(N);
    } else if (k == 1000) {
        Runner<1000>()(N);
    } else {
        throw std::invalid_argument("Value of k not compiled into this binary");
    }
}