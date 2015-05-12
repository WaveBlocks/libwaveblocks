#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/continuous_sqrt.hpp"

using namespace waveblocks;

int main()
{
    real_t PI = 3.14159265359;
    
    ContinuousSqrt<real_t> sqrt;
    
    //counter-clock-wise
    {
        std::ofstream out("cont_sqrt_ccw.csv");
        out << 
                std::setw(15) << "arg(in)" << " " << 
                std::setw(15) << "real(out)" << " " << 
                std::setw(15) << "imag(out)" << " " <<
                std::setw(15) << "abs(out)" << " " <<
                std::setw(15) << "arg(out)" << "\n";
        for (int i = 0; i <= 4*360; i += 10) {
            real_t angle = i*PI/real_t(180);
            
            complex_t result = sqrt( std::polar(1.0, angle) );
            
            out << 
                std::setw(15) << angle << " " << 
                std::setw(15) << std::real(result) << " " << 
                std::setw(15) << std::imag(result) << " " <<
                std::setw(15) << std::abs(result) << " " <<
                std::setw(15) << std::arg(result) << "\n";
        }
    }
    
    //clock-wise
    {
        std::ofstream out("cont_sqrt_cw.csv");
        out << 
                std::setw(15) << "arg(in)" << " " << 
                std::setw(15) << "real(out)" << " " << 
                std::setw(15) << "imag(out)" << " " <<
                std::setw(15) << "abs(out)" << " " <<
                std::setw(15) << "arg(out)" << "\n";
        for (int i = 4*360; i >= 0; i -= 10) {
            real_t angle = i*PI/real_t(180);
            
            complex_t result = sqrt( std::polar(1.0, angle) );
            
            out << 
                std::setw(15) << angle << " " << 
                std::setw(15) << std::real(result) << " " << 
                std::setw(15) << std::imag(result) << " " <<
                std::setw(15) << std::abs(result) << " " <<
                std::setw(15) << std::arg(result) << "\n";
        }
    }
    
    return 0;
}