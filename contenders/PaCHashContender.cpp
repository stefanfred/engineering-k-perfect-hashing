#include "PaCHashContender.hpp"

void paCHashContenderRunner(size_t N, size_t k) {
    if (k == 0) {
        throw std::invalid_argument("k must be greater than 0");
    }
    double aFrom = k / 4;
    double aTo = k * 4;
    for (double a = aFrom; a <= aTo; a *= 1.5) {
        {PaCHashContender(N, k, a).run();}
    }
}
