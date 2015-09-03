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
    
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n   ", "[", "]");
    
    const dim_t D = 5;
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;
    
    // (1) Define shapes
    HyperCubicShape<D> shape1({2,2,4,4,4});
    LimitedHyperbolicCutShape<D> shape23(9.0, {4,4,4,4,4});
    
    // (2) Enumerate shapes
    ShapeEnumerator<D,MultiIndex> shape_enumerator;
    
    std::vector< ShapeEnumSharedPtr<D,MultiIndex> > shape_enums;
    shape_enums.push_back( shape_enumerator.enumerate(shape1) );
    shape_enums.push_back( shape_enumerator.enumerate(shape23) );
    shape_enums.push_back( shape_enumerator.enumerate(shape23) );
    
    // (3) Initialize wavepacket-components
    HomogeneousHaWp<D,MultiIndex> wavepacket(3); // 3 = number of components
    wavepacket.eps() = 0.9;
    wavepacket.parameters() = HaWpParamSet<D>{};
    for (std::size_t c = 0; c < wavepacket.n_components(); c++) {
        wavepacket[c].shape() = shape_enums[c];
        
        std::size_t n_basis_shapes = shape_enums[c]->n_entries();
        
        // initialize wavepacket coefficients (with some bogus-values)
        wavepacket[c].coefficients() = std::vector<complex_t>(n_basis_shapes, 
                                                              complex_t(1.0/n_basis_shapes,1.0/n_basis_shapes));
    }
    
    // (4) Define quadrature points
    RMatrix<D, Eigen::Dynamic> grid(D,1);
    for (int i = 0; i < 1; i++) {
        grid(i,0) = -1.0 + 2.0*(i-1)/D;
    }
    
    // (5) Evaluate wavepacket-components
    
    std::cout << "Evaluate each component one by one ... " << std::endl;
    for (std::size_t c = 0; c < wavepacket.n_components(); c++) {
        CMatrix<1, Eigen::Dynamic> result = wavepacket[c].evaluate(grid);
        
        std::cout << "   " << result.format(CleanFmt) << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "Evaluate all components at once ... " << std::endl;
    CMatrix<Eigen::Dynamic, Eigen::Dynamic> result = wavepacket.evaluate(grid);
    std::cout << "   " << result.format(CleanFmt) << std::endl;
    
    std::cout << "Done" << std::endl;
    
    return 0;
}