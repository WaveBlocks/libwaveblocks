#include <iomanip>
#include <fstream>
#include <unordered_set>

#include "util/time.hpp"

#include "waveblocks/hypercubic_shape.hpp"
#include "waveblocks/hyperbolic_shape.hpp"
#include "waveblocks/shape_extension.hpp"

#include "waveblocks/sliced_shape_enumeration.hpp"

#include "waveblocks/hagedorn_wavepacket.hpp"
#include "waveblocks/gradient_operator.hpp"


#include "sample_wavepacket.hpp"

using namespace waveblocks;

int main(int argc, char* argv[])
{
    double start;

    const dim_t D = 3;
    typedef HyperCubicShape<D> S;
    
    S shape(MultiIndex<D> {{5,5,5}});
    
    auto parameters = createSampleParameters<D>();
    
    std::shared_ptr< ShapeEnumeration<D> > wave_enum( new SlicedShapeEnumeration<D,S>(shape) );
    std::shared_ptr< ShapeEnumeration<D> > grad_enum( new SlicedShapeEnumeration<D,ExtendedShape<D,S> >( ExtendedShape<D,S>{shape} ));
    
    auto wave_coeffs = createSampleCoefficients<D>(wave_enum);
    
    HagedornWavepacket<D> wavepacket(0.9, parameters, wave_enum, {{wave_coeffs}});
    
    GradientOperator<D> nabla(wave_enum, grad_enum);
    
    start = getRealTime();
    HagedornWavepacket<D,D> gradient = nabla(wavepacket);
    std::cout << "[TIME] apply gradient: " << (getRealTime() - start) << std::endl;
    
    start = getRealTime();
    std::cout << gradient(Eigen::Matrix<real_t,D,1>::Ones()) << std::endl;
    std::cout << "[TIME] evaluate gradient: " << (getRealTime() - start) << std::endl;
    
    //compare coefficients to csv file
    if (argc == 2) {
        std::cout << "compare gradient coefficients to reference file {" << std::endl;
        std::ifstream csv(argv[1]);
        if (csv.fail())
            std::cout << "   [ERROR] File not found!" << std::endl;
        
        std::size_t lines = 0;
        while (csv.good()) {
            ++lines;
            real_t tol = 1e-10;

            Eigen::Matrix<complex_t,D,1> solution;

            //read index
            MultiIndex<D> index;
            for (dim_t d = 0; d < D; d++) {
                int entry;
                csv >> entry;
                index[d] = entry;
            }

            //read real gradient part
            real_t temp[D];
            for (dim_t d = 0; d < D; d++) {
                real_t value;
                csv >> value;
                temp[d] = value;
            }

            //read imaginary gradient part
            for (dim_t d = 0; d < D; d++) {
                real_t value;
                csv >> value;
                solution(d,0) = complex_t(temp[d],value);
            }
            
            if (!grad_enum->contains(index)) {
                std::cout << "   [FAILURE] missing a node in our extended enumeration: " << index << std::endl;
                continue;
            }
            
            Eigen::Matrix<complex_t,D,1> result = gradient.coefficient( gradient.enumeration()->find(index) );
            
            real_t error = (result - solution).norm()/solution.norm();
            
            if ( error > tol)
                std::cout << "   [FAILURE] gradient at node " << index << " does not match entry in csv file. error = " << error << std::endl;
        }
        std::cout << "   [INFO] lines processed: " << lines << std::endl;
        std::cout << "}" << std::endl;
    }

    return 0;
}
