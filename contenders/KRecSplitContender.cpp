#include <KRecSplitContender.hpp>

#include "DispatchK.h"

template <size_t k>
struct KRecSplitContenderRunner {
    void operator() (size_t N) const {
        // TODO: Why does RecSplit overflow with k=10? Is this a wrong configuration?
        if constexpr (k <= 8) {
            {KRecSplitContender<k, k>(N, 2000).run();}
        }
    }
};

void kRecSplitContenderRunner(const size_t N, const size_t k) {
    dispatchK<KRecSplitContenderRunner>(k, N);
}
