#pragma once

#include <string>
#include <sstream>


namespace waveblocks{
    namespace utilities {
        template<class T>
        std::string to_string(T in) {
            std::ostringstream strs;
            strs << in;
            return strs.str();
        }
    }
}
