#include <vector>
#include <iostream>
#include <fstream>

#include "waveblocks.hpp"
#include "sample_wavefunc.hpp"

using namespace waveblocks;

template<dim_t D, class S>
void testSlicedShapeEnumeration(const SlicedShapeEnumeration<D,S> &enumeration)
{
    std::cout << "check enumeration {" << std::endl;
    
    {
        std::size_t ordinal = 0;
        for (auto index : enumeration) {
            if (ordinal != enumeration.find(index)) {
                std::cout << "   [FAILURE] find("<<index<<") != "<<ordinal << std::endl;
            }
            
            if (index != enumeration[ordinal]) {
                std::cout << "   [FAILURE] at("<<ordinal<<") != "<<index << std::endl;
            }
            ordinal++;
        }
        
        if (ordinal != enumeration.size()) {
            std::cout << "   [FAILURE] size() != "<<ordinal << std::endl;
        }
    }
    
    {
        std::size_t ordinal = 0;
        for (std::size_t islice = 0; islice < enumeration.slices().count(); islice++) {
            auto slice = enumeration.slice(islice);
            
            if (ordinal != slice.offset()) {
                std::cout << "   [FAILURE] slice_"<<islice<<".offset() != "<< ordinal << std::endl;
            }
            
            std::size_t ientry = 0;
            for (auto index : slice) {
                if (ientry != slice.find(index)) {
                    std::cout << "   [FAILURE] slice_"<<islice<<".find("<<index<< ") != "<<ientry << std::endl;
                }
                
                if (index != slice[ientry]) {
                    std::cout << "   [FAILURE] slice_"<<islice<<".at("<<ientry<< ") != "<<index << std::endl;
                }
                
                if (index != enumeration[ordinal]) {
                    std::cout << "   [FAILURE] order of slice_"<<islice<<" iteration != order of full iteration" << std::endl;
                }
                
                ++ientry;
                ++ordinal;
            }
            
            if (ientry != slice.size()) {
                std::cout << "   [FAILURE] slice_"<<islice<<".size() != "<<ientry << std::endl;
            }
        }
    }
    
    std::cout << "}" << std::endl;
}

int main()
{
    // test 2 dimensional
    const dim_t D = 5;
    
    HyperCubicShape<D> shape(MultiIndex<D>{{4,4,4,4,4}});
    
    SlicedShapeEnumeration<D,HyperCubicShape<D>> enumeration(shape);
    
    testSlicedShapeEnumeration(enumeration);
    
    HagedornParameterSet<D> parameters;
    createSampleParameters(parameters);
    
    std::cout << parameters << std::endl;
    
    std::vector<complex_t> coefficients(enumeration.size());
    
    //set coefficients
    //std::cout << "COEFFICIENTS: " << std::endl;
    createSampleCoefficients(enumeration, coefficients);
    
    // evaluate wavepacket at a chosen location
    {
        std::cout << "Test One Evaluation" << std::endl;
        
        Eigen::Matrix<real_t,D,1> x;
        for (dim_t d = 0; d < D; d++)
            x(d,0) = (d+1)/real_t(2*D);
        
        std::cout << "x: \n" << x << std::endl;
        std::cout << "psi: " << evaluateWavepacket(coefficients, parameters, enumeration, x) << '\n';
    }
    
    // plot wavepacket
//     std::size_t n1 = 20, n2 = 20;
//     double a1 = -5.0, b1 = 5.0;
//     double a2 = -5.0, b2 = 5.0;
//     
//     std::ofstream out("wavepacket.csv");
//     for (std::size_t i1 = 0; i1 <= n1; i1++) {
//         for (std::size_t i2 = 0; i2 <= n2; i2++) {
//             Eigen::Matrix<real_t,D,1> x;
//             
//             x[0] = a1 + i1*(b1-a1)/n1;
//             x[1] = a2 + i2*(b2-a2)/n2;
//             
//             out << x[0] << ' ';
//             out << x[1] << ' ';
//             
//             complex_t psi = evaluateWavepacket(coefficients, parameters, enumeration, x);
//             out << psi.real() << ' ';
//             out << psi.imag() << '\n';
//         }
//     }
//     out.close();
    
    std::cout << "compare results to reference file {" << std::endl;
    {
        std::ifstream in("../../test/wavepacket.csv");
        
        std::size_t lines = 0;
        while (in.good()) {
            ++lines;
            
            //read position
            Eigen::Matrix<real_t,D,1> x;
            for (dim_t d = 0; d < D; d++)
                in >> x(d,0);
            
            //read reference value
            real_t ref_real, ref_imag;
            in >> ref_real;
            in >> ref_imag;
            complex_t ref(ref_real, ref_imag);
            
            //compute wavepacket value
            complex_t psi = evaluateWavepacket(coefficients, parameters, enumeration, x);
            
            auto error = std::norm(psi - ref)/std::norm(ref);
            
            if (error > 10e-10) {
                std::cout << "   [FAILURE] mismatch at line " << lines << ". error = " << error << std::endl;
            }
        }
        
        std::cout << "   [INFO] processed " << lines << " lines" << std::endl;
    }
    std::cout << "}" << std::endl;
    
    return 0;
}