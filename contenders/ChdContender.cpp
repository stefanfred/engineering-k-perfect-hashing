#include "ChdContender.hpp"

void chdContenderRunner(size_t N, size_t k, double loadFactor) {
    if (k == 0) {
        throw std::invalid_argument("k must be greater than 0");
    }
    // CHD clamps the bucket sizes to [1, 15]
    for (size_t bucketSize = 1; bucketSize <= 15; bucketSize++) {
        {ChdContender(N, k, loadFactor, bucketSize, false).run();}
        // {ChdContender(N, k, loadFactor, bucketSize, true).run();} // Enters an endless loop
    }
}
