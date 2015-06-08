#include <iomanip>
#include <fstream>
#include <unordered_set>

#include "util/time.hpp"

#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/shape_hyperbolic.hpp"
#include "waveblocks/shape_extended.hpp"

#include "waveblocks/shape_enumeration_default.hpp"

#include "waveblocks/hagedorn_wavepacket.hpp"
#include "waveblocks/hagedorn_gradient_operator.hpp"

#include "waveblocks/tiny_multi_index.hpp"

#include "sample_wavepacket.hpp"

#include "check_shape_enumeration.hpp"
#include "check_wavepacket.hpp"
#include "check_coefficients.hpp"

using namespace waveblocks;

int main(int argc, char* argv[])
{
    double start, stop;

    const dim_t D = 5;

    typedef TinyMultiIndex<std::size_t,D> MultiIndex;
    typedef LimitedHyperbolicCutShape<D> S;
    
    S shape(9,{4,4,4,4,4});
    
    auto parameters = createSampleParameters<D>();

    std::shared_ptr< ShapeEnumeration<D> > wave_enum( new DefaultShapeEnumeration<D,MultiIndex,S>(shape) );
    
    std::shared_ptr< ShapeEnumeration<D> > grad_enum( new DefaultShapeEnumeration<D,MultiIndex,ExtendedShape<D,S> >( ExtendedShape<D,S> {shape} ));
    
    checkShapeEnumeration(*wave_enum, "wavepacket enumeration");
    checkShapeEnumeration(*grad_enum, "gradient enumeration");
    
    auto wave_coeffs = createSampleCoefficients<D>(wave_enum);
    
    HagedornWavepacket<D> wavepacket(0.9, parameters, wave_enum, {{wave_coeffs}});
    
    GradientOperator<D> nabla(grad_enum);
    
    start = getRealTime();
    std::array< HagedornWavepacket<D>, D> gradient = nabla(wavepacket);
    stop = getRealTime();
    std::cout << "[TIME] apply gradient: " << (stop - start) << std::endl;
    
    // evaluate wavepacket at a chosen location
    {
        std::cout << "evaluate gradient {" << std::endl;
        
        Eigen::Matrix<real_t,D,1> x;
        for (dim_t d = 0; d < D; d++)
            x(d,0) = (d+1)/real_t(2*D);
        
        start = getRealTime();
        for (auto & dir : gradient)
            std::cout << "   " << dir(x) << std::endl;
        stop = getRealTime();
        std::cout << "   time: " << (stop - start) << '\n';
        std::cout << "}" << std::endl;
    }
    
    //compare coefficients to csv file
    if (argc == 3) {
        compareCoefficientsToReferenceFile(gradient, argv[1]);
        //compareWavepacketToReferenceFile(gradient, argv[2]);
    }

    return 0;
}
