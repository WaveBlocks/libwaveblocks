#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Core>

#include "util/time.hpp"

#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/shape_hyperbolic.hpp"

#include "waveblocks/tiny_multi_index.hpp"
#include "waveblocks/shape_enum.hpp"
#include "waveblocks/shape_enumerator.hpp"

#include "waveblocks/hawp.hpp"

#include "sample_wavepacket.hpp"

#include "check_shape_enumeration.hpp"
#include "check_wavepacket.hpp"

using namespace waveblocks;

template<dim_t D, class MultiIndex>
struct Test
{
    double eps;
    HaWpParamSet<D> const& parameters;
    ShapeEnum<D, MultiIndex> const& shape_enum;
    std::vector<complex_t> const& coefficients;
    
    template<int N>
    void run(int npts)
    {
        Eigen::Matrix<complex_t,D,N> x(D,npts);
        for (dim_t d = 0; d < D; d++) {
            complex_t xd = { (d+1)/real_t(2*D), (D-d)/real_t(2*D) };
            
            for (int i = 0; i < npts; i++)
                x(d,i) = xd;
        }
        
        // (1) Evaluate wavepacket using reduce()
        
        Eigen::Array<complex_t,1,N> psi1 = hawp::basis(eps, &parameters, &shape_enum).at(x).reduce(coefficients);
        
        // (2) Evaluate wavepacket using all()
        
        {
            HaWpBasisVector<N> basis_vector = hawp::basis(eps, &parameters, &shape_enum).at(x).all();
            
            Eigen::Array<complex_t,1,N> psi2(1,npts);
            psi2.setZero();
            
            std::size_t j = 0;
            for (auto cj : coefficients) {
                psi2 += cj*basis_vector.row(j++);
            }
            
            std::cout << "   both lines should be identical" << std::endl;
            std::cout << "1. " << psi1 << std::endl;
            std::cout << "2. " << psi2 << std::endl;
        }
        
        //     std::cout << "   x: " << x.transpose() << '\n';
        //     std::cout << "   psi: " << psi1.transpose() << '\n';
        //     std::cout << "   psi (with prefactor): " << psi1.transpose()*hawp::prefactor(parameters) << std::endl;
        //     std::cout << "}" << std::endl;
    }
};

int main(int argc, char *argv[])
{
    const dim_t D = 3;
    typedef TinyMultiIndex<std::size_t,D> MultiIndex;
    typedef LimitedHyperbolicCutShape<D> S;
    
    S shape(7.0, {5,5,5});
    
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(shape);
    
    HaWpParamSet<D> parameters = createSampleParameters<D>();
    
    std::vector<complex_t> coefficients = createSampleCoefficients(shape_enum);
    
    double eps = 0.9;
    
    checkShapeEnumeration(shape_enum, "wavepacket shape enumeration");
    
    // evaluate wavepacket at a chosen location
    {
        std::cout << "chosen evaluation {" << std::endl;
        
        // Number of quadrature points can be set either at compile time or at runtime
        // One should get the same result either way!
        Test<D,MultiIndex> test = {eps, parameters, shape_enum, coefficients};
        
        std::cout << "\n# ---- static $N (number of quadrature points) ---- \n";
        std::cout << "\n#  $N = 1 \n";
        test.run<1>(1);
        std::cout << "\n#  $N = 5 \n";
        test.run<5>(5);
        
        std::cout << "\n# ---- dynamic $N (number of quadrature points) ---- \n";
        std::cout << "\n#  $N = 1 \n";
        test.run<Eigen::Dynamic>(1);
        std::cout << "\n#  $N = 5 \n";
        test.run<Eigen::Dynamic>(5);
        
        std::cout << "}\n\n" << std::endl;
    }
    
    if (argc == 2)
        compareWavepacketToReferenceFile(eps, parameters, shape_enum, coefficients, argv[1]);

    return 0;
}
