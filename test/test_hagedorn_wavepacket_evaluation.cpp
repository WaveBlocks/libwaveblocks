#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Core>

#include "util/time.hpp"

#include "waveblocks/hypercubic_shape.hpp"
#include "waveblocks/hagedorn_wavepacket.hpp"
#include "sample_wavepacket.hpp"

#include "check_shape_enumeration.hpp"
#include "check_wavepacket.hpp"

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

    std::shared_ptr< HagedornParameterSet<D> > parameters = createSampleParameters<D>();

    std::shared_ptr< ShapeEnumeration<D> > enumeration( new SlicedShapeEnumeration<D,S>(shape) );

    checkShapeEnumeration(*enumeration);

    std::shared_ptr< std::valarray<complex_t> > coefficients = createSampleCoefficients<D>(enumeration);

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

        std::cout << "   x: " << x.transpose() << '\n';
        std::cout << "   psi: " << psi.transpose() << '\n';
        std::cout << "   time: " << (stop - start) << '\n';
        std::cout << "}" << std::endl;
    }

    if (argc == 2)
        compareWavepacketToReferenceFile(wavepacket, argv[1]);

    return 0;
}
