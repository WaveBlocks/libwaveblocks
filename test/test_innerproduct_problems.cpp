#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/gauss_hermite_qr.hpp"
#include "waveblocks/genz_keister_qr.hpp"
#include "waveblocks/hawp_commons.hpp"
#include "waveblocks/hawp_paramset.hpp"
#include "waveblocks/homogeneous_inner_product.hpp"
#include "waveblocks/inhomogeneous_inner_product.hpp"
#include "waveblocks/vector_inner_product.hpp"
#include "waveblocks/shape_enumerator.hpp"
#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/stdarray2stream.hpp"
#include "waveblocks/tensor_product_qr.hpp"
#include "waveblocks/tiny_multi_index.hpp"
#include "waveblocks/utilities/timer.hpp"

using namespace waveblocks;

// Problem with TinyMultiIndex.
// Using "unsigned short" leads to
// "Assertion `value >= 0 && UINT(value) <= mask()' failed."
void runMultiD1()
{
    const dim_t D = 8;
    const real_t eps = 0.2;
    const dim_t n_coeffs = 5;
    const dim_t order = 5;
    using QR = GaussHermiteQR<order>;
    using TQR = TensorProductQR<QR,QR,QR,QR,QR,QR,QR,QR>;

    std::cout << D << "-D homogeneous quadrature of order " <<
        order << " with " << n_coeffs << " coefficients per dimension.\n";

    // Note the unsigned short here.
    using MultiIndex = TinyMultiIndex<unsigned short, D>;
    using IP = HomogeneousInnerProduct<D, MultiIndex, TQR>;

    // Set up sample wavepacket.
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum =
        enumerator.generate(HyperCubicShape<D>(n_coeffs));
    HaWpParamSet<D> param_set;
    Coefficients coeffs = Coefficients::Ones(std::pow(n_coeffs, D), 1);

    ScalarHaWp<D, MultiIndex> packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    // Calculate quadrature, volatile to force evaluation.
    volatile complex_t result = IP::quadrature(packet);
    (void)result;
}

// bad_alloc.
void runMultiD2()
{
    const dim_t D = 6;
    const real_t eps = 0.2;
    const dim_t n_coeffs = 6;
    const dim_t order = 6;
    using QR = GaussHermiteQR<order>;
    using TQR = TensorProductQR<QR,QR,QR,QR,QR,QR>;

    std::cout << D << "-D homogeneous quadrature of order " <<
        order << " with " << n_coeffs << " coefficients per dimension.\n";

    using MultiIndex = TinyMultiIndex<unsigned long, D>;
    using IP = HomogeneousInnerProduct<D, MultiIndex, TQR>;

    // Set up sample wavepacket.
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum =
        enumerator.generate(HyperCubicShape<D>(n_coeffs));
    HaWpParamSet<D> param_set;
    Coefficients coeffs = Coefficients::Ones(std::pow(n_coeffs, D), 1);

    ScalarHaWp<D, MultiIndex> packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    // Calculate quadrature, volatile to force evaluation.
    volatile complex_t result = IP::quadrature(packet);
    (void)result;
}

// TinyMultiIndex compiler warning.
// Using "unsigned long" leads to "left shift count >= width of type".
// No warning if using unsigned short.
// No runtime errors in either case.
void runMultiD3()
{
    const dim_t D = 1;
    const real_t eps = 0.2;
    const dim_t n_coeffs = 8;
    const dim_t order = 8;
    using QR = GaussHermiteQR<order>;
    using TQR = TensorProductQR<QR>;

    std::cout << D << "-D homogeneous quadrature of order " <<
        order << " with " << n_coeffs << " coefficients per dimension.\n";

    // Note the unsigned long here.
    using MultiIndex = TinyMultiIndex<unsigned long, D>;
    using IP = HomogeneousInnerProduct<D, MultiIndex, TQR>;

    // Set up sample wavepacket.
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum =
        enumerator.generate(HyperCubicShape<D>(n_coeffs));
    HaWpParamSet<D> param_set;
    Coefficients coeffs = Coefficients::Ones(std::pow(n_coeffs, D), 1);

    ScalarHaWp<D, MultiIndex> packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    // Calculate quadrature, volatile to force evaluation.
    volatile complex_t result = IP::quadrature(packet);
    (void)result;
}

int main()
{
    runMultiD1();
    runMultiD2();
    runMultiD3();

    return 0;
}
