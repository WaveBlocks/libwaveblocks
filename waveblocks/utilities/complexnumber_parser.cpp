#include "complexnumber_parser.hpp"

#include <cstdint>

namespace utilities {

bool parse_complex(char const* text, std::complex<double> & result)
{
    double real = 0.0;
    double imag = 0.0;

    do {
        char * end;
        double num = std::strtod(text, &end);

        if (end == text)
            return false; // require at least one processed character to avoid endless loop

            if (*end=='j' || *end=='i' || *end=='J' || *end=='I') {
                imag += num;
                ++end;
            } else {
                real += num;
            }

            text = end;
    } while (*text != 0);

    result = std::complex<double>(real, imag);

    return true;
}

}
