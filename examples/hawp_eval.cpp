
#include <iostream>

#include <Eigen/Core>

#include <waveblocks/util/timer.hpp>

#include <waveblocks/hawp_commons.hpp>
#include <waveblocks/hawp_gradient_operator.hpp>
#include <waveblocks/tiny_multi_index.hpp>

#include <waveblocks/shape_commons.hpp>
#include <waveblocks/shape_enumerator.hpp>



#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

using namespace waveblocks;

int pow(int m, int e)
{
    int res = 1;
    for (int i = 0; i < e; i++) {
        res *= m;
    }
    return res;
}

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    
    if (argc != 3) {
        std::cout << "2 arguments required: " << std::endl;
        std::cout << "   sparsity: integer value \\in {0,1, ... ,D}; higher value => more basis shape nodes" << std::endl;
        std::cout << "   basis shape limit: integer value <= 64/D" << std::endl;
        return -1;
    }
    
    int shape_sparsity = boost::lexical_cast<int>(argv[1]);
    int shape_limit = boost::lexical_cast<int>(argv[2]);
    
    // (1) Define dimensionality
    const dim_t D = 7;
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;
    
    // (2) Define shapes
    LimitedHyperbolicCutShape<D> shape(pow(2,D-shape_sparsity)*pow(3,shape_sparsity), shape_limit);
    
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "     CREATE WAVEPACKET" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    
    ScalarHaWp<D,MultiIndex> wp;
    wp.eps() = 0.9;
    
    // (3) Define parameters
    wp.parameters() = HaWpParamSet<D>{};
    wp.parameters().p = RMatrix<D,1>::Random()*0.1;
    wp.parameters().q = RMatrix<D,1>::Random()*0.1;
    wp.parameters().P = CMatrix<D,D>::Random();
    wp.parameters().Q = CMatrix<D,D>::Random();
    
    ShapeEnumerator<D,MultiIndex> enumerator;
    wp.shape() = enumerator.enumerate(shape);
    
    std::size_t bsize = wp.shape()->n_entries();
    wp.coefficients() = std::vector<complex_t>(bsize, complex_t(1.0/std::sqrt(bsize), 1.0/std::sqrt(bsize)));
    
    std::cout << wp.parameters() << std::endl;
    std::cout << std::endl;
    std::cout << "dimensionality: " << D << std::endl;
    std::cout << "basis shape:   " << shape << std::endl;
    std::cout << "    #entries:  " << wp.shape()->n_entries() << std::endl;
    std::cout << "    #slices:   " << wp.shape()->n_slices() << std::endl;
    
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "     EVALUATE WAVEPACKET" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n   ", "[", "]");
    
    const int numQ = 1;
    RMatrix<D,numQ> grid(D,1);
    std::cout << boost::format("Define %i quadrature points (columns)") % grid.cols() << std::endl;
    for (int i = 0; i < D; i++) {
        grid(i,0) = -1.0 + 2.0*(i)/(D-1);
    }
    
    std::cout << "   " << grid.format(CleanFmt) << std::endl;
    std::cout << std::endl;
    
    std::cout << boost::format("Evaluate wavepacket on %i quadrature points") % grid.cols() << std::endl;
    Timer timer;
    {
        CMatrix<1,numQ> result;
        
        timer.start();
        result = wp.evaluate(grid);
        timer.stop();
        std::cout << "   " << result.format(CleanFmt) << std::endl;
        std::cout << "   time: " << timer.millis() << " [ms] " << std::endl;
    }
    std::cout << std::endl;
}