#include <iomanip>
#include <fstream>
#include <unordered_set>

#include "waveblocks.hpp"
#include "util/time.hpp"

using namespace waveblocks;

int main()
{
    double start;
    
    const dim_t D = 6;
    typedef HyperCubicShape<D> Shape;
    
    Shape shape(MultiIndex<D>{{10,10,10,10,10,10}});
    
    start = getRealTime();
    SlicedShapeEnumeration<D,Shape> enumeration(shape);
    std::cout << "create enumeration: " << (getRealTime() - start) << std::endl;
    
    start = getRealTime();
    ShapeExtensionEnumeration<D,Shape> extension(shape);
    std::cout << "create extension: " << (getRealTime() - start) << std::endl;
    
    HagedornParameterSet<D> parameters;
    
    std::vector<complex_t> coefficients(enumeration.size());
    
    //std::cout << "coefficients: " << std::endl;
    start = getRealTime();
    {
        std::size_t ordinal = 0;
        for (auto index : enumeration) {
            int sum = 0;
            for (dim_t d = 0; d < D; d++)
                sum += index[d];
            
            real_t falloff = 0.1;
            
            coefficients[ordinal] = complex_t(
                std::exp(-falloff*sum),
                std::exp(-falloff*sum));
            
            //std::cout << index << ": " << coefficients[ordinal] << std::endl;
            
            ordinal++;
        }
    }
    std::cout << "create coefficients: " << (getRealTime() - start) << std::endl;
    
    std::vector<complex_t> gradient[D];
    Eigen::Matrix<real_t,D,1> x;
    
    //scatter-type
//     start = getRealTime();
//     for (dim_t d = 0; d < D; d++)
//         evaluateWavepacketGradient(coefficients, parameters, enumeration, extension, x, d, gradient[d]);
//     std::cout << "evaluate gradient (scatter type): " << (getRealTime() - start) << std::endl;
    
    //gather-type
    start = getRealTime();
    GradientOperator<D,Shape> nabla;
    for (dim_t d = 0; d < D; d++)
        nabla.evaluateWavepacketGradient(coefficients, parameters, enumeration, extension, x, d, gradient[d]);
    std::cout << "evaluate gradient (gather type): " << (getRealTime() - start) << std::endl;
    
//     std::size_t ordinal = 0;
//     {
//         for (auto index : enumeration) {
//             std::cout << index << "  ";
//             for (dim_t d = 0; d < D; d++) {
//                 complex_t value = gradient[d][ordinal];
//                 std::cout << std::setw(11) << value.real() << " ";
//             }
//             for (dim_t d = 0; d < D; d++) {
//                 complex_t value = gradient[d][ordinal];
//                 std::cout << std::setw(11) << value.imag() << " ";
//             }
//             std::cout << "\n";
//             ++ordinal;
//         }
//         
//         for (auto index : extension) {
//             std::cout << index << "  ";
//             for (dim_t d = 0; d < D; d++) {
//                 complex_t value = gradient[d][ordinal];
//                 std::cout << std::setw(11) << value.real() << " ";
//             }
//             for (dim_t d = 0; d < D; d++) {
//                 complex_t value = gradient[d][ordinal];
//                 std::cout << std::setw(11) << value.imag() << " ";
//             }
//             std::cout << "\n";
//             ++ordinal;
//         }
//     }
    
    return 0;
    
    //compare result to csv file 
    std::cout << "compare results to csv file {" << std::endl;
    std::ifstream csv("../../test/gradient.csv");
    if (csv.fail())
        std::cout << "   [ERROR] File not found!" << std::endl;
    
    std::size_t lines = 0;
    while (csv.good()) {
        ++lines;
        real_t tol = 1e-10;
        
        Eigen::Matrix<complex_t,D,1> solution;
        
        //read index
        MultiIndex<D> index;
        for (dim_t d = 0; d < D; d++)
            csv >> index[d];
        
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
        
        //compare
        std::size_t ordinal = enumeration.find(index);
        if (ordinal >= enumeration.size()) {
            ordinal = extension.find(index);
            if (ordinal >= extension.size()) {
                //check that solution is 0
                if (solution.norm() > tol)
                    std::cout << "   [FAILURE] value of foreign node " << index << " is non-zero in csv file" << std::endl;
                continue;
            }
            
            ordinal += enumeration.size();
        }
        
        Eigen::Matrix<complex_t,D,1> result;
        for (dim_t d = 0; d < D; d++)
            result[d] = gradient[d][ordinal];
        
        real_t error = (result - solution).norm()/solution.norm();
        
        if ( error > tol)
            std::cout << "   [FAILURE] gradient at node " << index << " does not match entry in csv file. error = " << error << std::endl;
    }
    std::cout << "   [INFO] lines processed: " << lines << std::endl;
    std::cout << "}" << std::endl;
    
    return 0;
}