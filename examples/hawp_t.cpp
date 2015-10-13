#include <iostream>
#include <algorithm>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include <omp.h>

#include <waveblocks/util/timer.hpp>
#include <waveblocks/hawp_commons.hpp>
#include <waveblocks/hawp_gradient_operator.hpp>
#include <waveblocks/tiny_multi_index.hpp>
#include <waveblocks/shape_commons.hpp>
#include <waveblocks/shape_enumerator.hpp>


using namespace waveblocks;


template<dim_t D, class MultiIndex>
void print_coefficients(AbstractScalarHaWp<D,MultiIndex> const& wp)
{
    {
        for (int i = 0; i < wp.shape()->n_slices(); i++) {
            auto const& slice = wp.shape()->slice(i);

            for (std::size_t j = 0; j < slice.size(); j++) {
                std::cout << "   " << slice[j] << ": " << wp.coefficients()[slice.offset() + j] << std::endl;
            }
        }
    }
    std::cout << std::endl;
}



int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;

    // (1) Define dimensionality
    const dim_t D = 1;
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;

    // (2) Create Wavepacket
    ScalarHaWp<D,MultiIndex> wp;
    wp.eps() = 0.1;

    // (3) Define standard parameters
    wp.parameters() = HaWpParamSet<D>{};

    wp.parameters().Q(wp.parameters().Q() / 2);
    wp.parameters().P(wp.parameters().P() * 2);
    wp.parameters().S(complex_t(1,0));

    // (4) Define basis shape
    HyperCubicShape<D> shape(12);
    ShapeEnumerator<D,MultiIndex> enumerator;
    wp.shape() = enumerator.enumerate(shape);

    // (5) Define coefficients
    std::size_t bsize = wp.shape()->n_entries();
    wp.coefficients() = Coefficients(bsize);
    wp.coefficients()[0] = 1.0;


    // // wtf??
    // std::vector<MultiIndex> order(bsize);
    // for (int s = 0; s < wp.shape()->n_slices(); s++) {
    //     auto const& slice = wp.shape()->slice(s);
    //     for (std::size_t i = 0; i < slice.size(); i++) {
    //         order[slice.offset()+i] = slice[i];
    //     }
    // }

    // std::less<MultiIndex> cmp;
    // std::sort(order.begin(), order.end(), cmp);
    // for (std::size_t i = 0; i < bsize; i++) {
    //     int islice = 0;
    //     for (int d = 0; d < D; d++)
    //         islice += order[i][d];
    //     auto const& slice = wp.shape()->slice(islice);
    //     std::size_t ordinal = slice.offset() + slice.find(order[i]);
    //     wp.coefficients()[ordinal] = 1.0;
    // }

    complex_t prefactor = wp.prefactor();
    complex_t phase = wp.phasefactor();

    // The final wavepacket
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "     WAVEPACKET CONFIGURATION" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "dimensionality: " << D << std::endl;
    std::cout << "basis shape:    " << shape << std::endl;
    std::cout << "    #entries:   " << wp.shape()->n_entries() << std::endl;
    std::cout << "    #slices:    " << wp.shape()->n_slices() << std::endl;
    std::cout << "scaling eps:    " << wp.eps() << std::endl;
    std::cout << "parameterset:   " << wp.parameters();
    std::cout << "prefactor:      " << prefactor << std::endl;
    std::cout << "global phase:   " << phase << std::endl;
    std::cout << "coefficients:" << std::endl;
    print_coefficients(wp);

    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "     EVALUATE WAVEPACKET" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;

    Eigen::IOFormat CleanFmt(16, 0, ", ", "\n   ", "[", "]");

    // Construct a grid of points
    const int numQ = 5;
    RMatrix<D,numQ> grid(D,numQ);
    std::cout << boost::format("Define %i quadrature points (columns)") % grid.cols() << std::endl;
    for (int d = 0; d < D; d++) {
        for (int n = 0; n < numQ; n++) {
            grid(d,n) = -1.0 + 2.0*n/float(numQ-1);
        }
    }
    std::cout << "   " << grid.format(CleanFmt) << std::endl;
    std::cout << std::endl;

    // Evaluate Wavepacket
    std::cout << boost::format("Evaluate wavepacket on %i quadrature points") % grid.cols() << std::endl;

    CMatrix<1,numQ> result(1,numQ);
    result = prefactor * phase * wp.evaluate(grid);
    std::cout << "   " << result.format(CleanFmt) << std::endl;

    std::cout << std::endl;
}
