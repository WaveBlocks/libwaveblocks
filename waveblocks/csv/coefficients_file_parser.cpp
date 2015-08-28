#include "coefficients_file_parser.hpp"

#include <cstdlib>
#include <string>
#include <iostream>

#include "../util/complexnumber_parser.hpp"

namespace waveblocks {
namespace csv {

CoefficientsFileParser::CoefficientsFileParser(std::string file, int dimensions, int components)
    : lattice_node(dimensions)
    , coefficients(components)
    , in_(file.c_str())
    , wp_dimensions_(dimensions)
    , wp_components_(components)
    , line_number_(0)
{
    if (!in_.good())
        throw std::runtime_error("cannot find/open file: "+file);
}

bool CoefficientsFileParser::next()
{
    ++line_number_;
    
    for (int i = 0; i < wp_dimensions_; i++) {
        in_ >> lattice_node[i];
        
        if (in_.eof())
            return false;
    }
    
    for (int i = 0; i < wp_components_; i++) {
        std::string token;
        in_ >> token;
        
        if (in_.eof())
            return false;
        
        // python encloses complex numbers in parentheses
        if (token.size() >= 2 && token.front() == '(' && token.back() == ')')
            token = token.substr(1, token.size()-2);
        
        if (!util::parse_complex(token.c_str(), coefficients[i]))
            throw std::runtime_error("parse error at line: "+std::to_string(line_number_));
    }
    
    if (in_.bad())
        throw std::runtime_error("parse error at line: "+std::to_string(line_number_));
    
    return in_.good();
}

}
}