#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "waveblocks/types.hpp"
#include "waveblocks/wavepackets/hawp_commons.hpp"
#include "waveblocks/wavepackets/hawp_paramset.hpp"
#include "waveblocks/wavepackets/shapes/tiny_multi_index.hpp"
#include "waveblocks/wavepackets/shapes/shape_enumerator.hpp"
#include "waveblocks/wavepackets/shapes/shape_hypercubic.hpp"
#include "waveblocks/innerproducts/gauss_hermite_qr.hpp"
#include "waveblocks/innerproducts/tensor_product_qr.hpp"
#include "waveblocks/innerproducts/genz_keister_qr.hpp"
#include "waveblocks/innerproducts/homogeneous_inner_product.hpp"
#include "waveblocks/innerproducts/inhomogeneous_inner_product.hpp"
#include "waveblocks/innerproducts/vector_inner_product.hpp"
#include "waveblocks/utilities/stdarray2stream.hpp"


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
    using QR = innerproducts::GaussHermiteQR<order>;
    using TQR = innerproducts::TensorProductQR<QR,QR,QR,QR,QR,QR,QR,QR>;

    std::cout << D << "-D homogeneous quadrature of order " <<
        order << " with " << n_coeffs << " coefficients per dimension.\n";

    // Note the unsigned short here.
    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;
    using IP = innerproducts::HomogeneousInnerProduct<D, MultiIndex, TQR>;

    // Set up sample wavepacket.
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(n_coeffs));
    wavepackets::HaWpParamSet<D> param_set;
    Coefficients coeffs = Coefficients::Ones(std::pow(n_coeffs, D), 1);

    wavepackets::ScalarHaWp<D, MultiIndex> packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
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
    using QR = innerproducts::GaussHermiteQR<order>;
    using TQR = innerproducts::TensorProductQR<QR,QR,QR,QR,QR,QR>;

    std::cout << D << "-D homogeneous quadrature of order " <<
        order << " with " << n_coeffs << " coefficients per dimension.\n";

    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned long, D>;
    using IP = innerproducts::HomogeneousInnerProduct<D, MultiIndex, TQR>;

    // Set up sample wavepacket.
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(n_coeffs));
    wavepackets::HaWpParamSet<D> param_set;
    Coefficients coeffs = Coefficients::Ones(std::pow(n_coeffs, D), 1);

    wavepackets::ScalarHaWp<D, MultiIndex> packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
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
    using QR = innerproducts::GaussHermiteQR<order>;
    using TQR = innerproducts::TensorProductQR<QR>;

    std::cout << D << "-D homogeneous quadrature of order " <<
        order << " with " << n_coeffs << " coefficients per dimension.\n";

    // Note the unsigned long here.
    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned long, D>;
    using IP = innerproducts::HomogeneousInnerProduct<D, MultiIndex, TQR>;

    // Set up sample wavepacket.
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(n_coeffs));
    wavepackets::HaWpParamSet<D> param_set;
    Coefficients coeffs = Coefficients::Ones(std::pow(n_coeffs, D), 1);

    wavepackets::ScalarHaWp<D, MultiIndex> packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
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
