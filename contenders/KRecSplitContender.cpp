#include <KRecSplitContender.hpp>

#include "DispatchK.h"

template <size_t k, size_t leafSize>
void dispatchBucketSize(size_t N) {
    for (double bucketSizeFactor = 0.2; bucketSizeFactor < 1; bucketSizeFactor += 0.2) {
        size_t bucketSize = bucketSizeFactor * kphf::RecSplit::MAX_BUCKET_SIZE<k, leafSize>;
        {KRecSplitContender<k, leafSize>(N, bucketSize).run();}
    }
}

template <size_t k>
struct KRecSplitContenderRunner {
    void operator() (size_t N) const {
        dispatchBucketSize<k, 1>(N);
        dispatchBucketSize<k, 2>(N);
        if constexpr (k <= 100) {
            dispatchBucketSize<k, 4>(N);
        }
        if constexpr (k <= 50) {
            dispatchBucketSize<k, 6>(N);
        }
    }
};

void kRecSplitContenderRunner(const size_t N, const size_t k) {
    dispatchK<KRecSplitContenderRunner>(k, N);
}
