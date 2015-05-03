#include <iomanip>
#include <fstream>
#include <unordered_set>

#include "waveblocks.hpp"
#include "util/time.hpp"

#include "sample_wavefunc.hpp"

using namespace waveblocks;

int main()
{
    double start;

    const dim_t D = 3;
    typedef HyperCubicShape<D> S;

    S shape(MultiIndex<D> {{10,10,10}});

    auto parameters = createSampleParameters<D>();
    auto enumeration = std::make_shared< SlicedShapeEnumeration<D,S> >(shape);
    auto coefficients = createSampleCoefficients<D,S>(enumeration);
    
    HagedornWavepacket<D,S> wavepacket(parameters, coefficients);
    
    GradientOperator<D,S> nabla(shape, enumeration);
    
    //gather-type
    start = getRealTime();
    auto gradient = nabla(wavepacket, 0);
    std::cout << "[TIME] evaluate gradient (gather type): " << (getRealTime() - start) << std::endl;
    
//     //compare result to csv file
//     std::cout << "compare results to reference file {" << std::endl;
//     std::ifstream csv("../../test/gradient.csv");
//     if (csv.fail())
//         std::cout << "   [ERROR] File not found!" << std::endl;
// 
//     std::size_t lines = 0;
//     while (csv.good()) {
//         ++lines;
//         real_t tol = 1e-10;
// 
//         Eigen::Matrix<complex_t,D,1> solution;
// 
//         //read index
//         MultiIndex<D> index;
//         for (dim_t d = 0; d < D; d++) {
//             int entry;
//             csv >> entry;
//             index[d] = entry;
//         }
// 
//         //read real gradient part
//         real_t temp[D];
//         for (dim_t d = 0; d < D; d++) {
//             real_t value;
//             csv >> value;
//             temp[d] = value;
//         }
// 
//         //read imaginary gradient part
//         for (dim_t d = 0; d < D; d++) {
//             real_t value;
//             csv >> value;
//             solution(d,0) = complex_t(temp[d],value);
//         }
// 
//         //compare
//         std::size_t ordinal = enumeration.find(index);
//         if (ordinal >= enumeration.size()) {
//             ordinal = extension.find(index);
//             if (ordinal >= extension.size()) {
//                 //check that solution is 0
//                 if (solution.norm() > tol)
//                     std::cout << "   [FAILURE] value of foreign node " << index << " is non-zero in csv file" << std::endl;
//                 continue;
//             }
// 
//             ordinal += enumeration.size();
//         }
// 
//         Eigen::Matrix<complex_t,D,1> result;
//         for (dim_t d = 0; d < D; d++)
//             result[d] = gradient[d][ordinal];
// 
//         real_t error = (result - solution).norm()/solution.norm();
// 
//         if ( error > tol)
//             std::cout << "   [FAILURE] gradient at node " << index << " does not match entry in csv file. error = " << error << std::endl;
//     }
//     std::cout << "   [INFO] lines processed: " << lines << std::endl;
//     std::cout << "}" << std::endl;

    return 0;
}
