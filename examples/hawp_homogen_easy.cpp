#include <Eigen/Core>

#include <iostream>

#include "waveblocks/tiny_multi_index.hpp"

#include "waveblocks/shape_commons.hpp"

#include "waveblocks/shape_enum.hpp"
#include "waveblocks/shape_enumerator.hpp"

#include "waveblocks/hawp_commons.hpp"

using namespace waveblocks;

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    
    const dim_t D = 5;
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;
    
    // (1) Define shapes
    HyperCubicShape<D> shape1({4,4,4,2,2});
    HyperCubicShape<D> shape2({2,2,4,4,4});
    LimitedHyperbolicCutShape<D> shape3(9.0, {4,4,4,4,4});
    
    // (2) Enumerate shapes
    ShapeEnumerator<D,MultiIndex> shape_enumerator;
    std::vector< ShapeEnumSharedPtr<D,MultiIndex> > shape_enums;
    shape_enums.push_back( std::make_shared< ShapeEnum<D,MultiIndex> >(shape_enumerator.enumerate(shape1)) );
    shape_enums.push_back( std::make_shared< ShapeEnum<D,MultiIndex> >(shape_enumerator.enumerate(shape2)) );
    shape_enums.push_back( std::make_shared< ShapeEnum<D,MultiIndex> >(shape_enumerator.enumerate(shape3)) );
    
    // (3) Initialize wavepacket-components
    HomogeneousHaWp<D,MultiIndex> wavepacket(3); // 3 = number of components
    wavepacket.eps() = 0.9;
    for (std::size_t c = 0; c < wavepacket.n_components(); c++) {
        wavepacket[c].shape() = shape_enums[c];
        
        wavepacket[c].coefficients().resize(shape_enums[c]->n_entries());
    }
    
    // (4) Define quadrature points
    int npts = 5; // number of quadrature points
    const int N = Eigen::Dynamic; // number of quadrature points (template parameter)
    
    Eigen::Matrix<complex_t, D, N> grid(D,npts);
    
    for (int n = 0; n < npts; n++) {
        
    }
    
    // (5) Evaluate wavepacket-components
    for (std::size_t c = 0; c < wavepacket.n_components(); c++) {
        Eigen::Matrix<complex_t, 1, N> value = wavepacket[c].evaluate(grid);
        
        std::cout << "component [" << c << "]: ";
        for (int n = 0; n < npts; n++) {
            std::cout << value(0,n) << " ";
        }
        std::cout << std::endl;
    }
    
    return 0;
}