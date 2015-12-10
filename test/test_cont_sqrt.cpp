#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/math/continuous_sqrt.hpp"


using namespace waveblocks;

int main()
{
    real_t PI = 3.14159265359;
    int N = 1000;

    math::ContinuousSqrt<real_t> sqrt;

    std::ofstream out("cont_sqrt.csv");
    out.precision(15);

    for (int i = 0; i < N; i++) {
        real_t angle = 4*PI * i / (N-1);
        complex_t y = std::polar(2.0+0.1*angle, angle);
        complex_t sy = sqrt(y);

        out << std::setw(20) << angle << ", " <<
            std::setw(20) << std::real(y) << ", " <<
            std::setw(20) << std::imag(y) << ", " <<
            std::setw(20) << std::real(sy) << ", " <<
            std::setw(20) << std::imag(sy) << std::endl;
    }

    return 0;
}
