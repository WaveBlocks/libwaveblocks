#include <iomanip>
#include <fstream>
#include <unordered_set>

#include "util/time.hpp"

#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/shape_hyperbolic.hpp"
#include "waveblocks/shape_extended.hpp"

#include "waveblocks/shape_enum.hpp"
#include "waveblocks/shape_enumerator.hpp"

#include "waveblocks/hawp.hpp"

#include "waveblocks/tiny_multi_index.hpp"

#include "sample_wavepacket.hpp"

#include "check_shape_enumeration.hpp"
//#include "check_wavepacket.hpp"
#include "check_coefficients.hpp"

using namespace waveblocks;

int main(int argc, char* argv[])
{
    double start, stop;

    const dim_t D = 5;

    typedef TinyMultiIndex<std::size_t,D> MultiIndex;
    typedef LimitedHyperbolicCutShape<D> S;
    
    S shape(9,{4,4,4,4,4});
    
    HagedornParameterSet<D> parameters = createSampleParameters<D>();
    
    ShapeEnumerator<D,MultiIndex> enumerator;
    
    ShapeEnum<D,MultiIndex> wave_enum = enumerator.generate(shape);
    ShapeEnum<D,MultiIndex> grad_enum = enumerator.generate( ExtendedShape<D,S>{shape} );
    
    checkShapeEnumeration(wave_enum, "wavepacket enumeration");
    checkShapeEnumeration(grad_enum, "gradient enumeration");
    
    std::vector<complex_t> wave_coeffs = createSampleCoefficients<D>(wave_enum);
    
    double eps = 0.9;
    
    auto basis = hawp::basis(eps, &parameters, &wave_enum);
    
    start = getRealTime();
    std::array< std::vector<complex_t>, D> grad_coeffs = hawp::nabla(basis, &grad_enum).apply(wave_coeffs);
    stop = getRealTime();
    std::cout << "[TIME] apply gradient: " << (stop - start) << std::endl;
    
    // evaluate wavepacket at a chosen location
    {
        std::cout << "evaluate gradient {" << std::endl;
        
        Eigen::Matrix<real_t,D,1> x;
        for (dim_t d = 0; d < D; d++)
            x(d,0) = (d+1)/real_t(2*D);
        
        start = getRealTime();
        
        auto evaluator = basis.at(x);
        
        for (auto & dir : grad_coeffs)
            std::cout << "   " << evaluator.reduce(dir) << std::endl;
        
        stop = getRealTime();
        std::cout << "   time: " << (stop - start) << '\n';
        std::cout << "}" << std::endl;
    }
    
    //compare coefficients to csv file
    if (argc == 3) {
        compareCoefficientsToReferenceFile(argv[1], &grad_enum, grad_coeffs);
        //compareWavepacketToReferenceFile(gradient, argv[2]);
    }

    return 0;
}
