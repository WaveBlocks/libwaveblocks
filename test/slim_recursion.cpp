#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include "util/time.hpp"

#include "waveblocks/hypercubic_shape.hpp"
#include "waveblocks/hagedorn_wavepacket.hpp"
#include "sample_wavepacket.hpp"

#include "test_shape_enumeration.hpp"

using namespace waveblocks;

template<dim_t D>
MultiIndex<D> createFilledMultiIndex(int entry)
{
    MultiIndex<D> index;
    for (dim_t i = 0; i < D; i++)
        index[i] = entry;
    return index;
}

int main(int argc, char *argv[])
{
    const dim_t D = 3;
    typedef HyperCubicShape<D> S;
    
    S shape(createFilledMultiIndex<D>(5));
    
    std::shared_ptr< const HagedornParameterSet<D> > parameters = createSampleParameters<D>();
    
    std::shared_ptr< const ShapeEnumeration<D> > enumeration( new SlicedShapeEnumeration<D,S>(shape) );
    
    testSlicedShapeEnumeration(*enumeration);
    
    std::shared_ptr< const std::valarray<complex_t> > coefficients = createSampleCoefficients<D>(enumeration);
    
    HagedornWavepacket<D> wavepacket(0.9, parameters, enumeration, {coefficients});
    
    // evaluate wavepacket at a chosen location
    {
        std::cout << "chosen evaluation {" << std::endl;
        
        double start = getRealTime();
        Eigen::Matrix<real_t,D,1> x;
        for (dim_t d = 0; d < D; d++)
            x(d,0) = (d+1)/real_t(2*D);
        double stop = getRealTime();
        
        auto psi = wavepacket(x);
        
        std::cout << "   psi: " << psi << '\n';
        std::cout << "   time: " << (stop - start) << '\n';
        std::cout << "}" << std::endl;
    }
    
    /**
     * read space delimited csv files with columns x_1 x_2 .. x_d real(psi[x]) imag(psi[x])
     */
    std::cout << "compare results to reference files {" << std::endl;
    for (int iarg = 1; iarg < argc; iarg++) {
        std::cout << "   [FILE] " << argv[iarg] << std::endl;
        std::ifstream in(argv[iarg]);
        
        std::size_t lines = 0;
        while (in.good()) {
            ++lines;
            
            //read position
            Eigen::Matrix<real_t,D,1> x;
            for (dim_t d = 0; d < D; d++)
                in >> x(d,0);
            
            //read reference value
            real_t ref_real, ref_imag;
            in >> ref_real;
            in >> ref_imag;
            complex_t ref(ref_real, ref_imag);
            
            //compute wavepacket value
            complex_t psi = wavepacket(x)(0,0);
            
            auto error = std::norm(psi - ref)/std::norm(ref);
            
            if (error > 10e-10) {
                std::cout << "      [FAILURE] mismatch at line " << lines << ". error = " << error << std::endl;
            }
        }

        std::cout << "      [INFO] processed " << lines << " lines" << std::endl;
    }
    std::cout << "}" << std::endl;

    return 0;
}
