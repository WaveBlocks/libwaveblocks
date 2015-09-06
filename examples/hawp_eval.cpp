
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

#include <algorithm>

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
        std::cout << "   basis shape limit: integer value <= 2^floor(64/D)" << std::endl;
        return -1;
    }
    
    int shape_sparsity = boost::lexical_cast<int>(argv[1]);
    int shape_limit = boost::lexical_cast<int>(argv[2]);
    
    // (1) Define dimensionality
    const dim_t D = 6;
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;
    
    // (2) Define shapes
    LimitedHyperbolicCutShape<D> shape(pow(2,D-shape_sparsity)*pow(3,shape_sparsity), shape_limit);
    
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "     CREATE WAVEPACKET" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    
    ScalarHaWp<D,MultiIndex> wp;
    wp.eps() = 0.9;
    
    // (3) Define "random" parameters
    wp.parameters() = HaWpParamSet<D>{};
    for (int i = 0; i < D; i++) {
        wp.parameters().p(i,0) += 0.1*std::sin(0.9*i);
        wp.parameters().q(i,0) += 0.1*std::cos(0.8*i);
        for (int j = 0; j < D; j++) {
            wp.parameters().P(i,j) += 0.1*std::exp(complex_t(0,0.7*(i*D+j)));
            wp.parameters().Q(i,j) += 0.1*std::exp(complex_t(0,0.6*(i*D+j)));
        }
    }
    std::cout << wp.parameters() << std::endl;
    ShapeEnumerator<D,MultiIndex> enumerator;
    wp.shape() = enumerator.enumerate(shape);
    
    std::size_t bsize = wp.shape()->n_entries();
    wp.coefficients() = std::vector<complex_t>(bsize);
    
    // (4) Define "random" coefficients
    std::vector<MultiIndex> order(bsize);
    for (int islice = 0; islice < wp.shape()->n_slices(); islice++) {
        auto const& slice = wp.shape()->slice(islice);
        for (std::size_t i = 0; i < slice.size(); i++) {
            order[slice.offset()+i] = slice[i];
        }
    }
    
    std::less<MultiIndex> cmp;
    std::sort(order.begin(), order.end(), cmp);
    for (std::size_t i = 0; i < bsize; i++) {
        int islice = 0;
        for (int d = 0; d < D; d++)
            islice += order[i][d];
        
        auto const& slice = wp.shape()->slice(islice);
        
        std::size_t ordinal = slice.offset() + slice.find(order[i]);
        
        wp.coefficients()[ordinal] = std::exp(complex_t(0,i))/std::sqrt(bsize);
    
        //std::cout << order[i] << wp.coefficients()[i] << std::endl;
    }
    
    //std::cout << wp.parameters() << std::endl;
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
    RMatrix<D,numQ> grid(D,numQ);
    std::cout << boost::format("Define %i quadrature points (columns)") % grid.cols() << std::endl;
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < numQ; j++) {
            grid(i,j) = (-1.0 + 2.0*(i)/(D-1))*(1+j)/numQ;
        }
    }
    
    std::cout << "   " << grid.format(CleanFmt) << std::endl;
    std::cout << std::endl;
    
//     std::cout << boost::format("Evaluate basis functions on %i quadrature points") % grid.cols() << std::endl;
//     {
//         CMatrix<Eigen::Dynamic, numQ> result; // size = ( number of basis functions X number of quadrature points)
//         
//         result = wp.evaluate_basis(grid);
//         
//         for (int i = 0; i < wp.shape()->n_slices(); i++) {
//             auto const& slice = wp.shape()->slice(i);
//             
//             for (std::size_t j = 0; j < slice.size(); j++) {
//                 std::cout << "   " << slice[j] << ": " << result.row(slice.offset() + j) << std::endl;
//             }
//         }
//     }
//     std::cout << std::endl;
    
    
    std::cout << boost::format("Evaluate wavepacket on %i quadrature points") % grid.cols() << std::endl;
    Timer timer;
    {
        CMatrix<1,numQ> result(1,numQ);
        
        timer.start();
        //for (int i = 0; i < 1000; i++) {
            result = wp.evaluate(grid);
        //}
        timer.stop();
        std::cout << "   " << result.format(CleanFmt) << std::endl;
        std::cout << "   time: " << timer.millis() << " [ms] " << std::endl;
    }
    std::cout << std::endl;
}