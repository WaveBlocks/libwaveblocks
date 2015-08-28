#ifndef WAVEBLOCKS_CSV_HAWP_COEFFS_LOADER_HPP
#define WAVEBLOCKS_CSV_HAWP_COEFFS_LOADER_HPP

#include <complex>
#include <fstream>
#include <memory>
#include <iostream>

#include "coefficients_file_parser.hpp"

#include "../shape_enum.hpp"

#include "boost/format.hpp"

namespace waveblocks {
namespace csv {
    
template<dim_t D, class MultiIndex>
class HaWpCoefficientsLoader
{
public:
    std::vector<complex_t> operator()(std::shared_ptr<ShapeEnum<D,MultiIndex> > shape, std::string filename)
    {
        std::vector<complex_t> coeffs(shape->n_entries());
        
        std::vector<int> occurrences(coeffs.size());
        
        CoefficientsFileParser parser(filename, D, 1);
        
        while (parser.next()) {
            MultiIndex index;
            
            int i_slice = 0;
            for (int i = 0; i < D; i++) {
                int entry = parser.lattice_node[i];
                if (entry < 0)
                    throw std::runtime_error((boost::format("negative index at line %i") % parser.line_number()).str());
                
                index[i] = entry;
                if (index[i] != entry)
                    throw std::runtime_error((boost::format("integer overflow at line %i") % parser.line_number()).str());
                
                i_slice += entry;
            }
            
            std::size_t ordinal;
            if (i_slice >= shape->n_slices() || !shape->slice(i_slice).try_find(index, ordinal))
                throw std::runtime_error((boost::format("multi-index at line %i is not element of shape") % parser.line_number()).str());
            
            ordinal += shape->slice(i_slice).offset();
            
            if (occurrences[ordinal] >= 1)
                throw std::runtime_error((boost::format("duplicate entry at line %i") % parser.line_number()).str());
            
            coeffs[ordinal] = parser.coefficients[0];
            occurrences[ordinal] += 1;
        }
        
        for (std::size_t i = 0; i < occurrences.size(); i++) {
            if (occurrences[i] == 0) {
                throw std::runtime_error((boost::format("missing entry for lattice node: %s") % shape->at(i)).str());
            }
        }
        
        return coeffs;
    }
};

}
}

#endif