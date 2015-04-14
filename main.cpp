#include <iostream>

#include <Eigen/Core>

#include "waveblocks.hpp"
#include <fstream>

//#define BOOST_TEST_MODULE MyTest
//#include <boost/test/unit_test.hpp>

using namespace waveblocks;

void example()
{
    /*
    HyperbolicCutShape<3> shape(3.0);
    
    LexicalShapeIterator<3,HyperbolicCutShape<3>> it(shape);
    
    for (unsigned x = 0; x <= shape.getSurface(0, MultiIndex<3>{{0,0,0}}); x++) {
        for (unsigned y = 0; y <= shape.getSurface(1, MultiIndex<3>{{x,0,0}}); y++) {
            for (unsigned z = 0; z <= shape.getSurface(2, MultiIndex<3>{{x,y,0}}); z++) {
                MultiIndex<3> index = {x,y,z};
                BOOST_CHECK(index == it.getMultiIndex());
                it.advance();
            }
        }
    }*/
}

int main() {
    example();
    
    /*
    const std::size_t D = 30;
    
    LimitedHyperbolicCutShape<D> shape(13.0, MultiIndex<D>{{13,13,13,13,13, 13,13,13,13,13, 13,13,13,13,13, 13,13,13,13,13, 13,13,13,13,13, 13,13,13,13,13}});
    
    LexicalShapeIterator<D,LimitedHyperbolicCutShape<D>> it(shape);
    
    LexicalShapeEnumeration<D,LimitedHyperbolicCutShape<D>> enumeration(shape, 5);
    
    std::cout << "size: " << enumeration.size() << std::endl;
    do {
        auto comp = enumeration.find(it.getMultiIndex());
        if (it != comp)
            std::cout << it << " vs " << comp << std::endl;
    } while (it.advance());*/
    
    const std::size_t D = 2;
    
    typedef Eigen::Matrix<complex_t,D,D> CMatrix;
    typedef Eigen::Matrix<complex_t,D,1> CVector;
    
    typedef Eigen::Matrix<real_t,D,D> RMatrix;
    typedef Eigen::Matrix<real_t,D,1> RVector;
    
    HyperCubicShape<D> shape(MultiIndex<D>{{3,3}});
    
    LexicalShapeEnumeration<D,HyperCubicShape<D>> enumeration(shape, 5);
    
    HagedornWavepacket<D,HyperCubicShape<D>> packet(enumeration);
    
    SlicedShapeEnumeration<D,HyperCubicShape<D>> slices(shape);
    
    /*
    for (auto slice : slices) {
        std::cout << "---" << std::endl;
        for (auto index : slice) {
            std::cout << index << std::endl;
        }
    }
    
    std::cout << slices[4].size() << std::endl;
    std::cout << slices[4].at(2) << std::endl;
    std::cout << slices[4].find(MultiIndex<D>{{2,3}}) << std::endl;
    */
    
    LexicalShapeIterator<D,HyperCubicShape<D>> it(shape);
    
    do {
        std::cout << it.getMultiIndex() << std::endl;
    } while (it.advance());
    
    HagedornParameterSet<D> parameters;
    //parameters.eps = 0.9;
    //parameters.q << -1, 0;
    parameters.p << 2, 0;
    //parameters.Q << complex_t(2,2), complex_t(-1,-1), complex_t(1,1), complex_t(1,-3);
    //parameters.P << complex_t(1,2), complex_t(1,1), complex_t(-1,1), complex_t(-2,3);
    
    std::cout << parameters << std::endl;
    
    std::vector<complex_t> coefficients(slices[slices.count()-1].offset() + slices[slices.count()-1].size());
    
    coefficients[0] = complex_t(1.0,0.0);
    
    std::size_t n1 = 20, n2 = 20;
    double a1 = -5.0, b1 = 5.0;
    double a2 = -5.0, b2 = 5.0;
    
    std::ofstream out("wavepacket.csv");
    for (std::size_t i = 0; i < n1; i++) {
        for (std::size_t j = 0; j < n2; j++) {
            Eigen::Matrix<real_t,D,1> x;
            
            x[0] = a1 + i*(b1-a1)/n1;
            x[1] = a2 + j*(b2-a2)/n2;
            
            out << x[0] << ' ';
            out << x[1] << ' ';
            out << evaluateWavepacket(coefficients, parameters, slices, x).real() << '\n';
        }
    }
    out.close();
    
    //CMatrix A, B;
    
    //A << complex_t(1,2), complex_t(2,1), complex_t(-1,-1), complex_t(0,1);
    //B << complex_t(-1,0), complex_t(1,1), complex_t(1,2), complex_t(-1,-1);
    
    //std::cout << A*B << std::endl;
    
    return 0;
}