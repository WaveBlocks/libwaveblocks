#ifndef YAML_COMPLEX_HPP
#define YAML_COMPLEX_HPP

#include <complex>

#include <yaml-cpp/yaml.h>

#include "../util/complexnumber_parser.hpp"

namespace YAML {
    template<>
    struct convert<std::complex<double> > {
        static Node encode(const std::complex<double> & rhs) {
            if (rhs.imag() == 0.0) {
                return Node(std::to_string(rhs.real()));
            }
            else if (rhs.real() == 0.0) {
                    return Node(""+std::to_string(rhs.imag())+"j");
            }
            else if (rhs.imag() < 0.0) {
                return Node(""+std::to_string(rhs.real())+std::to_string(rhs.imag())+"j");
            }
            else {
                return Node(""+std::to_string(rhs.real())+"+"+std::to_string(rhs.imag())+"j");
            }
        }

        static bool decode(const Node& node, std::complex<double> & rhs) {
            std::string str = node.as<std::string>();

            // deal with python's annoying habbit to add parentheses to complex numbers
            if (str.size() >= 2 && str.front() == '(' && str.back() == ')')
                str = str.substr(1, str.size()-2);

            return util::parse_complex(str.c_str(), rhs);
        }
    };
}

#endif
