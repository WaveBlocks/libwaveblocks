#include <iostream>

#include <Eigen/Core>

#include "waveblocks.hpp"

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
    
    /*
    packet.parameters.eps = 0.9;
    packet.parameters.q << -1, 0;
    packet.parameters.p << 1, 1;
    packet.parameters.Q << complex_t(2,2), complex_t(-1,-1), complex_t(1,1), complex_t(1,-3);
    packet.parameters.P << complex_t(1,2), complex_t(1,1), complex_t(-1,1), complex_t(-2,3);
    
    std::cout << packet.parameters.q << std::endl;
    std::cout << packet.parameters.p << std::endl;
    std::cout << packet.parameters.Q << std::endl;
    std::cout << packet.parameters.P << std::endl;*/
    
    //CMatrix A, B;
    
    //A << complex_t(1,2), complex_t(2,1), complex_t(-1,-1), complex_t(0,1);
    //B << complex_t(-1,0), complex_t(1,1), complex_t(1,2), complex_t(-1,-1);
    
    //std::cout << A*B << std::endl;
    
    return 0;
}