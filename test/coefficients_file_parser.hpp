#ifndef WAVEBLOCKS_COEFFICIENTS_FILE_PARSER_HPP
#define WAVEBLOCKS_COEFFICIENTS_FILE_PARSER_HPP

#include <vector>
#include <complex>
#include <fstream>

namespace waveblocks
{

class CoefficientsFileParser
{
public:
    std::vector<int> lattice_node;
    std::vector<std::complex<double> > coefficients;
    
    CoefficientsFileParser(std::string file, int wp_dimensions, int wp_components);
    
    bool next();
    
private:
    std::fstream in_;
    int wp_dimensions_;
    int wp_components_;
};

}

#endif