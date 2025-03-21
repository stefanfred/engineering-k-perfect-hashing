#include "PaCHashContender.hpp"

void paCHashContenderRunner(size_t N, size_t k) {
    if (k == 0) {
        throw std::invalid_argument("k must be greater than 0");
    }
    for (int i = 1; i < 40; i++) {
        {PaCHashContender(N, k, i / 20.0).run();}
    }
}
