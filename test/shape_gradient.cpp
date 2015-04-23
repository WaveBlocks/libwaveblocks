#include <iomanip>
#include <fstream>
#include <unordered_set>

#include "waveblocks.hpp"

using namespace waveblocks;

int main()
{
    const dim_t D = 5;
    typedef HyperCubicShape<D> Shape;
    
    Shape shape(MultiIndex<D>{{10,10,10,10,10}});
    
    SlicedShapeEnumeration<D,Shape> enumeration(shape);
    ShapeExtensionEnumeration<D,Shape> extension(shape);
    
    HagedornParameterSet<D> parameters;
    
    std::vector<complex_t> coefficients(enumeration.size());
    
    //std::cout << "coefficients: " << std::endl;
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
    
    std::vector<complex_t> gradient[D];
    
    Eigen::Matrix<real_t,D,1> x;
    for (dim_t d = 0; d < D; d++)
        evaluateWavepacketGradient(coefficients, parameters, enumeration, extension, x, d, gradient[d]);
    
    if (false) {
    std::size_t ordinal = 0;
    {
        for (auto index : enumeration) {
            complex_t value = gradient[0][ordinal++];
            std::cout << index << "  " << value.real() << " " << std::setw(10) << value.imag() << "\n";
        }
        
        for (auto index : extension) {
            complex_t value = gradient[0][ordinal++];
            std::cout << index << "  " << std::setw(10) << value.real() << " " << std::setw(10) << value.imag() << "\n";
        }
    }
    }
    
    //compare result to csv file 
    std::ifstream csv("../../test/gradient.csv");
    while (csv.good()) {
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
    
    return 0;
}