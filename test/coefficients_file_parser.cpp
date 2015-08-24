#include "coefficients_file_parser.hpp"

#include <cstdlib>
#include <string>

using namespace waveblocks;

CoefficientsFileParser::CoefficientsFileParser(std::string file, int dimensions, int components)
    : lattice_node(dimensions)
    , coefficients(components)
    , in_(file.c_str())
    , wp_dimensions_(dimensions)
    , wp_components_(components)
{ }

bool CoefficientsFileParser::next()
{
    std::size_t line_number = 0;
    std::string line;
    
    while (in_.good()) {
        std::getline(in_, line);
        
        // filter empty lines
        if (line.size() == 0)
            continue;
        
        // filter comments
        if (line[0] == '#')
            continue;
        
        char const* pos = line.c_str();
        
        char * end;
        
        // read multi-indices
        for (int i = 0; i < wp_dimensions_; i++) {
            long int value = std::strtol(pos, &end, 10);
            if (end == pos)
                throw std::runtime_error("parse error at line: "+line_number);
            
            lattice_node[i] = value;
            
            pos = end;
        }
        
        // read coefficients
        for (int i = 0; i < wp_dimensions_; i++) {
            std::complex<double> value = std::strtod(pos, &end);
            if (end == pos)
                throw std::runtime_error("parse error at line: "+line_number);
            
            coefficients[i] = value;
            
            pos = end;
        }
        
        std::strtod(pos, &end);
        if (pos != end)
            throw std::runtime_error("parse error at line: "+line_number);
        
        ++line_number;
    }
    
    return in_.good();
}
