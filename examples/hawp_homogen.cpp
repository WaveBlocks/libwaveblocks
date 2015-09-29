#include <Eigen/Core>

#include <iostream>

#include "waveblocks/tiny_multi_index.hpp"

#include "waveblocks/shape_commons.hpp"

#include "waveblocks/shape_enum.hpp"
#include "waveblocks/shape_enumerator.hpp"

#include "waveblocks/hawp_commons.hpp"

#include "waveblocks/util/timer.hpp"

using namespace waveblocks;

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    
    Eigen::IOFormat CleanFmt(8, 0, ", ", "\n   ", "[", "]");
    
    const dim_t D = 8;
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;
    
    // (1) Define shapes
    HyperCubicShape<D> shape1({2,2,4,4,4,5,3,5});
    LimitedHyperbolicCutShape<D> shape23(1 << D, 4);
    
    // (2) Enumerate shapes
    ShapeEnumerator<D,MultiIndex> shape_enumerator;
    
    int n_components = 10;
    
    std::vector< ShapeEnumSharedPtr<D,MultiIndex> > shape_enums(n_components);
    shape_enums[0] = shape_enumerator.enumerate(shape1);
    for (int i = 1; i < n_components; i++) {
        shape_enums[i] = shape_enumerator.enumerate(shape23);
    }
    
    // (3) Initialize wavepacket-components
    HomogeneousHaWp<D,MultiIndex> wavepacket(n_components);
    
    wavepacket.eps() = 0.9;
    wavepacket.parameters() = HaWpParamSet<D>{};
    wavepacket.parameters().p = RMatrix<D,1>::Random();
    wavepacket.parameters().q = RMatrix<D,1>::Random();
    wavepacket.parameters().P = CMatrix<D,D>::Random();
    wavepacket.parameters().Q = CMatrix<D,D>::Random();
    for (std::size_t c = 0; c < wavepacket.n_components(); c++) {
        wavepacket[c].shape() = shape_enums[c];
        
        std::size_t n_basis_shapes = shape_enums[c]->n_entries();

        Coefficients coeffs(n_basis_shapes);
        for (std::size_t i = 0; i < n_basis_shapes; i++) {
            coeffs[i] = std::exp(complex_t(0,i))/std::sqrt(n_basis_shapes);
        }
        
        wavepacket[c].coefficients() = std::move(coeffs);
    }
    
    // (4) Define quadrature points
    const int numQ = 1;
    RMatrix<D, numQ> grid(D,1);
    for (int i = 0; i < D; i++) {
        grid(i,0) = -1.0 + 2.0*(i-1)/D;
    }
    
    // (5) Evaluate wavepacket-components
    Timer timer;
    
    std::cout << "Evaluate each component one by one ... " << std::endl;
    {
        double cumm_time = 0.0;
        for (std::size_t c = 0; c < wavepacket.n_components(); c++) {
            timer.start();
            CMatrix<1,numQ> result = wavepacket[c].evaluate(grid);
            timer.stop();
            cumm_time += timer.millis();
            std::cout << "   " << result.format(CleanFmt) << std::endl;
            //std::cout << "   time: " << timer.millis() << " [ms]" << std::endl;
        }
        std::cout << "   time: " << cumm_time << " [ms]" << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "Evaluate all components at once ... " << std::endl;
    {
        timer.start();
        CMatrix<Eigen::Dynamic,numQ> result = wavepacket.evaluate(grid);
        timer.stop();
        std::cout << "   " << result.format(CleanFmt) << std::endl;
        std::cout << "   time: " << timer.millis() << " [ms]" << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "Done" << std::endl;
    
    return 0;
}