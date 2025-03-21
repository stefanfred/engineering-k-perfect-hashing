#include "ChdContender.hpp"

void chdContenderRunner(size_t N, size_t k, double loadFactor) {
    if (k == 0) {
        throw std::invalid_argument("k must be greater than 0");
    }
    for (int keysPerBucket = 1; keysPerBucket < 8; keysPerBucket++) {
        {ChdContender(N, k, loadFactor, loadFactor, keysPerBucket, false).run();}
        {ChdContender(N, k, loadFactor, loadFactor, keysPerBucket, true).run();}
    }
}
