#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Core>

#include "util/time.hpp"

#include "waveblocks/shape_hyperbolic.hpp"
#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/hagedorn_wavepacket.hpp"
#include "waveblocks/tiny_multi_index.hpp"
#include "sample_wavepacket.hpp"

#include "check_shape_enumeration.hpp"
#include "check_wavepacket.hpp"

using namespace waveblocks;

int main()
{
    const dim_t D = 3;
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;
    typedef HyperCubicShape<D> S;
    
    S shape({7,7,7,7,7, 7,7,7,7,7});
    
    std::shared_ptr< HagedornParameterSet<D> > parameters = createSampleParameters<D>();
    
    std::shared_ptr< ShapeEnumeration<D> > enumeration( new DefaultShapeEnumeration<D,MultiIndex,S>(shape) );
    
    std::shared_ptr< std::valarray<complex_t> > coefficients = createSampleCoefficients<D>(enumeration);
    
    HagedornWavepacket<D> wavepacket(0.9, parameters, enumeration, {coefficients});

    // evaluate wavepacket at a chosen location
    {
        std::cout << "chosen evaluation {" << std::endl;
        
        Eigen::Matrix<complex_t,D,1> x;
        for (dim_t d = 0; d < D; d++) {
            x(d,0) = complex_t( (d+1)/real_t(2*D), (D-d)/real_t(2*D));
        }
        
        double start = getRealTime();
        
        auto psi = wavepacket(x);
        
        double stop = getRealTime();
        
        std::cout << "   x: " << x.transpose() << '\n';
        std::cout << "   psi: " << psi.transpose() << '\n';
        std::cout << "   psi (with prefactor): " << psi.transpose()*wavepacket.prefactor() << std::endl;
        std::cout << "   time: " << (stop - start) << '\n';
        std::cout << "}" << std::endl;
    }

    return 0;
}
