#include <iostream>
#include <complex>
#include <vector>

double poisson_cdf(int x, double lambda) {
    double term = std::exp(-lambda);
    double sum = term;

    for (int k = 1; k <= x; ++k) {
        term *= lambda / k;
        sum += term;
    }

    return sum;
}

int main() {
    // note that this only works for small k because of numerical instabilities
    size_t k = 100;

    std::vector<double> res;
    double step = 0.0001;
    double x = 0;
    double y = 1;
    double integral = 0;

    while (y > 0 && x < 1) {
        y = poisson_cdf(k - 1, integral);
        integral += k*step / y;
        res.push_back(y);
        x += step;

        // probabilities
        //std::cout << "RESULT x=" << x << " y=" << y << std::endl;
    }


    std::vector<double> resInt;
    resInt.push_back(0);
    double prev = 0;
    for (auto v: res) {
        prev += std::log(v);
        resInt.push_back(prev);
    }
    for (size_t i = 0; i < resInt.size(); ++i) {
        // bucket assignment function
        std::cout << "RESULT x=" << (i * step) << " y=" << resInt[i]/resInt[resInt.size()-1] << std::endl;
    }

    return 0;
}

