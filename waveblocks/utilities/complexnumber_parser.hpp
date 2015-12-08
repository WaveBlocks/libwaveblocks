#pragma once

#include <complex>


namespace utilities
{
    /**
     *
     * This function basically invokes std::strtod until end of string and
     * accumulates those comverted numbers.
     *
     * If one of the characters 'i', 'j','I', 'J' immediately follows
     * after std::strtod, the value is interpreted as a complex number.
     *
     * Be aware, that this function does not check input properly.
     * Thus for pathological input you get unexpected results,
     * but this function is able to correctly interpret output of python/numpy.
     *
     * Valid input:
     *   "1.0+2.0j" yields (1.0,2.0)
     *   "0.7j+1.1+2.2+3.3" yields (6.6,0.7)
     *   "1e3-1.7e-3j" yield (1000,-0.0017)
     *
     * Invalid input:
     *   "-+2.0"
     *
     * \param text zero-terminated c-string
     */
    bool parse_complex(char const* text, std::complex<double> & result);
}
