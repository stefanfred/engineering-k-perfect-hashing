#include "PaCHashContender.hpp"

void paCHashContenderRunner(size_t N, size_t k) {
    if (k == 0) {
        throw std::invalid_argument("k must be greater than 0");
    }
    for (double a = 0.6 * k; a <= 2.0 * k; a += 0.2 * k) {
        {PaCHashContender(N, k, a).run();}
    }
}
