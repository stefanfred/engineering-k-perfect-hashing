#include <KRecSplitContender.hpp>

#include "PaCHashContender.hpp"

void kRecSplitContenderRunner(size_t N, size_t k) {
    if (k == 0) {
        throw std::invalid_argument("k must be greater than 0");
    }
    // TODO: Make this nicer using templated dispatch function.
    if (k == 2) {
        {KRecSplitContender<2, 2>(N, 2000).run();}
    } else if (k == 4) {
        {KRecSplitContender<4, 4>(N, 2000).run();}
    } else if (k == 8) {
        // TODO: These are probably wrong configurations. They are super slow.
        {KRecSplitContender<8, 8>(N, 2000).run();}
    }
    // TODO: Larger k cause overflow in template evaluation
}
