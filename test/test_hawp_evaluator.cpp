#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Core>

#include "util/time.hpp"

#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/shape_hyperbolic.hpp"

#include "waveblocks/tiny_multi_index.hpp"
#include "waveblocks/shape_enum.hpp"
#include "waveblocks/shape_enumerator.hpp"

#include "waveblocks/hawp.hpp"

#include "sample_wavepacket.hpp"

#include "check_shape_enumeration.hpp"
#include "check_wavepacket.hpp"

using namespace waveblocks;

int main(int argc, char *argv[])
{
    const dim_t D = 3;
    typedef TinyMultiIndex<std::size_t,D> MultiIndex;
    typedef LimitedHyperbolicCutShape<D> S;
    
    S shape(7.0, {5});
    
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(shape);
    
    HaWpParamSet<D> parameters = createSampleParameters<D>();
    
    std::vector<complex_t> coefficients = createSampleCoefficients(shape_enum);
    
    double eps = 0.9;
    
    checkShapeEnumeration(shape_enum, "wavepacket shape enumeration");
    
    // evaluate wavepacket at a chosen location
    {
        std::cout << "chosen evaluation {" << std::endl;
        
        Eigen::Matrix<complex_t,D,1> x;
        for (dim_t d = 0; d < D; d++) {
            x(d,0) = complex_t( (d+1)/real_t(2*D), (D-d)/real_t(2*D) );
        }
        
        double start = getRealTime();
        
        Eigen::Array<complex_t,1,1> psi = hawp::basis(eps, &parameters, &shape_enum).at(x).reduce(coefficients);
        
        double stop = getRealTime();
        
        std::cout << "   x: " << x.transpose() << '\n';
        std::cout << "   psi: " << psi.transpose() << '\n';
        std::cout << "   psi (with prefactor): " << psi.transpose()*hawp::prefactor(parameters) << std::endl;
        std::cout << "   time: " << (stop - start) << '\n';
        std::cout << "}" << std::endl;
    }
    
    if (argc == 2)
        compareWavepacketToReferenceFile(eps, parameters, shape_enum, coefficients, argv[1]);

    return 0;
}
