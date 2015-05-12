#include <iomanip>
#include <fstream>
#include <unordered_set>

#include "util/time.hpp"

#include "waveblocks/hypercubic_shape.hpp"
#include "waveblocks/hyperbolic_shape.hpp"
#include "waveblocks/shape_extension.hpp"

#include "waveblocks/sliced_shape_enumeration.hpp"

#include "waveblocks/hagedorn_wavepacket.hpp"
#include "waveblocks/hagedorn_gradient_operator.hpp"

#include "sample_wavepacket.hpp"

#include "check_shape_enumeration.hpp"
#include "check_wavepacket.hpp"
#include "check_coefficients.hpp"

using namespace waveblocks;

int main(int argc, char* argv[])
{
    double start;

    const dim_t D = 4;

    typedef HyperbolicCutShape<D> S;

    S shape(7.0);

    auto parameters = createSampleParameters<D>();

    std::shared_ptr< ShapeEnumeration<D> > wave_enum( new SlicedShapeEnumeration<D,S>(shape) );
    std::shared_ptr< ShapeEnumeration<D> > grad_enum( new SlicedShapeEnumeration<D,ExtendedShape<D,S> >( ExtendedShape<D,S> {shape} ));

    checkShapeEnumeration(*wave_enum);
    checkShapeEnumeration(*grad_enum);

    auto wave_coeffs = createSampleCoefficients<D>(wave_enum);

    HagedornWavepacket<D> wavepacket(0.9, parameters, wave_enum, {{wave_coeffs}});

    GradientOperator<D> nabla(wave_enum, grad_enum);


    start = getRealTime();
    HagedornWavepacket<D,D> gradient = nabla(wavepacket);
    std::cout << "[TIME] apply gradient: " << (getRealTime() - start) << std::endl;

    // evaluate wavepacket at a chosen location
    {
        std::cout << "chosen evaluation {" << std::endl;

        double start = getRealTime();
        Eigen::Matrix<real_t,D,1> x;
        for (dim_t d = 0; d < D; d++)
            x(d,0) = (d+1)/real_t(2*D);
        double stop = getRealTime();

        auto psi = gradient(x);

        std::cout << "   psi: " << psi.transpose() << '\n';
        std::cout << "   time: " << (stop - start) << '\n';
        std::cout << "}" << std::endl;
    }

    //compare coefficients to csv file
    if (argc == 3) {
        compareCoefficientsToReferenceFile(gradient, argv[1]);
        compareWavepacketToReferenceFile(gradient, argv[2]);
    }

    return 0;
}
