#include "waveblocks/hypercubic_shape.hpp"
#include "waveblocks/slim_basis_recursion.hpp"

#include <vector>
#include <iostream>
#include <fstream>

using namespace waveblocks;

template<std::size_t D>
void test(std::size_t limit)
{
    MultiIndex<D> limits;
    
    for (std::size_t d = 0; d < D; d++)
        limits[d] = limit;
    
    HyperCubicShape<D> shape(limits);
    
    SlicedShapeEnumeration<D,HyperCubicShape<D>> slices(shape);
    
    HagedornParameterSet<D> parameters;
    
    std::vector<complex_t> coefficients(slices[slices.count()-1].offset() + slices[slices.count()-1].size());
    
    Eigen::Matrix<real_t,D,1> x;
    
    std::cout << evaluateWavepacket(coefficients, parameters, slices, x) << std::endl;
}

int main()
{
    // test 2 dimensional
    const std::size_t D = 2;
    
    typedef Eigen::Matrix<complex_t,D,D> CMatrix;
    typedef Eigen::Matrix<complex_t,D,1> CVector;
    
    typedef Eigen::Matrix<real_t,D,D> RMatrix;
    typedef Eigen::Matrix<real_t,D,1> RVector;
    
    HyperCubicShape<D> shape(MultiIndex<D>{{4,4}});
    
    LexicalShapeEnumeration<D,HyperCubicShape<D>> enumeration(shape, 5);
    
    SlicedShapeEnumeration<D,HyperCubicShape<D>> slices(shape);
    
    HagedornParameterSet<D> parameters;
    parameters.eps = 0.9;
    //parameters.q << -1, 0;
    parameters.p << 2, 0;
    //parameters.Q << complex_t(2,2), complex_t(-1,-1), complex_t(1,1), complex_t(1,-3);
    //parameters.P << complex_t(1,2), complex_t(1,1), complex_t(-1,1), complex_t(-2,3);
    
    std::cout << parameters << std::endl;
    
    std::vector<complex_t> coefficients(slices[slices.count()-1].offset() + slices[slices.count()-1].size());
    
    //set coefficients
    //std::cout << "COEFFICIENTS: " << std::endl;
    for (auto slice : slices) {
        std::size_t i = slice.offset();
        for (auto index : slice) {
            coefficients[i++] = complex_t(std::exp(-double(0.5*index[0])),std::exp(-double(0.5*index[1])));
        }
    }
    
    ////DEBUG
    std::cout << "Test One Evaluation" << std::endl;
    Eigen::Matrix<real_t,D,1> x;
    x << 0.5, -0.3;
    std::cout << "x: \n" << x << std::endl;
    std::cout << "psi: " << evaluateWavepacket(coefficients, parameters, slices, x) << '\n';
    
    std::size_t n1 = 20, n2 = 20;
    double a1 = -5.0, b1 = 5.0;
    double a2 = -5.0, b2 = 5.0;
    
    std::ofstream out("slim_recursion_cpp.csv");
    for (std::size_t i1 = 0; i1 <= n1; i1++) {
        for (std::size_t i2 = 0; i2 <= n2; i2++) {
            Eigen::Matrix<real_t,D,1> x;
            
            x[0] = a1 + i1*(b1-a1)/n1;
            x[1] = a2 + i2*(b2-a2)/n2;
            
            out << x[0] << ' ';
            out << x[1] << ' ';
            
            complex_t psi = evaluateWavepacket(coefficients, parameters, slices, x);
            out << psi.real() << ' ';
            out << psi.imag() << '\n';
        }
    }
    out.close();
    
    return 0;
}