#include <KRecSplitContender.hpp>

#include "DispatchK.h"

template <size_t k, size_t leafSize>
void dispatchBucketSize(size_t N) {
    for (double bucketSizeFactor = 0.2; bucketSizeFactor < 1; bucketSizeFactor += 0.2) {
        size_t bucketSize = bucketSizeFactor * sux::function::krecsplit::MAX_BUCKET_SIZE<k, leafSize>;
        {KRecSplitContender<k, leafSize>(N, bucketSize).run();}
    }
}

template <size_t k>
struct KRecSplitContenderRunner {
    void operator() (size_t N) const {
        // TODO: Constexpr lookup tables get compile time error for k=1000
        if constexpr (k <= 100) {
            dispatchBucketSize<k, 2>(N);
            dispatchBucketSize<k, 4>(N);
            dispatchBucketSize<k, 6>(N);
        }
    }
};

void kRecSplitContenderRunner(const size_t N, const size_t k) {
    dispatchK<KRecSplitContenderRunner>(k, N);
}
