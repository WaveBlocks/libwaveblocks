#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/gauss_hermite_qr.hpp"
#include "waveblocks/hawp_commons.hpp"
#include "waveblocks/hawp_paramset.hpp"
#include "waveblocks/homogeneous_inner_product.hpp"
#include "waveblocks/inhomogeneous_inner_product.hpp"
#include "waveblocks/shape_enumerator.hpp"
#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/stdarray2stream.hpp"
#include "waveblocks/tensor_product_qr.hpp"
#include "waveblocks/tiny_multi_index.hpp"

using namespace waveblocks;

void test1DGaussHermite()
{
    std::cout << "\n1D Homogeneous Gauss-Hermite Test\n";
    std::cout <<   "---------------------------------\n";

    const real_t eps = 0.2;
    const dim_t D = 1;
    const dim_t N = 10;
    const dim_t order = 8;
    using MultiIndex = TinyMultiIndex<unsigned short, D>;
    using QR = GaussHermiteQR<order>;

    // Set up sample 1D wavepacket.
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum =
        enumerator.generate(HyperCubicShape<D>(N));
    HaWpParamSet<D> param_set;
    std::cout << param_set << std::endl;
    std::vector<complex_t> coeffs(N, 1.0);

    // Print QR nodes and weights.
    std::cout << "nodes: {";
    for (auto x : QR::nodes()) std::cout << " " << x;
    std::cout << " }\n";

    std::cout << "weights: {";
    for (auto x : QR::weights()) std::cout << " " << x;
    std::cout << " }\n";

    ScalarHaWp<D, MultiIndex> packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    // Calculate inner product matrix, print it.
    HomogeneousInnerProduct<D, MultiIndex, QR> ip;
    CMatrix<Eigen::Dynamic, Eigen::Dynamic> mat =
        ip.build_matrix(packet);

    std::cout << "IP matrix:\n" << mat << std::endl;

    // Calculate quadrature.
    std::cout << "Quadrature: " << ip.quadrature(packet) << "\n";
}

void test1DGaussHermiteOperator()
{
    std::cout << "\n1D Homogeneous Gauss-Hermite Custom Operator Test\n";
    std::cout <<   "-------------------------------------------------\n";

    const real_t eps = 0.2;
    const dim_t D = 1;
    const dim_t N = 10;
    const dim_t order = 8;
    using MultiIndex = TinyMultiIndex<unsigned short, D>;
    using QR = GaussHermiteQR<order>;
    using CMatrix1N = CMatrix<1, Eigen::Dynamic>;
    using CMatrixD1 = CMatrix<D, 1>;
    using CMatrixDN = CMatrix<D, Eigen::Dynamic>;

    // Set up sample 1D wavepacket.
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum =
        enumerator.generate(HyperCubicShape<D>(N));
    HaWpParamSet<D> param_set;
    std::cout << param_set << std::endl;
    std::vector<complex_t> coeffs(N, 1.0);

    // Print QR nodes and weights.
    std::cout << "nodes: {";
    for (auto x : QR::nodes()) std::cout << " " << x;
    std::cout << " }\n";

    std::cout << "weights: {";
    for (auto x : QR::weights()) std::cout << " " << x;
    std::cout << " }\n";

    ScalarHaWp<D, MultiIndex> packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    // Calculate inner product matrix, print it.
    // Use an operator that returns a sequence 1, 2, ..., number_nodes.
    HomogeneousInnerProduct<D, MultiIndex, QR> ip;
    auto op =
        [] (CMatrixDN nodes, CMatrixD1 pos) -> CMatrix1N
    {
        (void)pos;
        const dim_t n_nodes = nodes.cols();
        CMatrix1N result(1, n_nodes);
        for(int i = 0; i < n_nodes; ++i) result(0, i) = i+1;
        return result;
    };
    CMatrix<Eigen::Dynamic, Eigen::Dynamic> mat =
        ip.build_matrix(packet, op);

    std::cout << "IP matrix:\n" << mat << std::endl;

    // Calculate quadrature.
    std::cout << "Quadrature: " << ip.quadrature(packet, op) << "\n";
}

void test3DGaussHermite()
{
    std::cout << "\n3D Homogeneous Tensor-Product Gauss-Hermite Test\n";
    std::cout <<   "------------------------------------------------\n";

    const real_t eps = 0.2;
    const dim_t D = 3;
    const dim_t N = 5;
    using MultiIndex = TinyMultiIndex<unsigned short, D>;

    // Set up sample 3D wavepacket.
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum =
        enumerator.generate(HyperCubicShape<D>(N));
    HaWpParamSet<D> param_set;
    std::cout << param_set << std::endl;
    std::vector<complex_t> coeffs(std::pow(N, D), 1.0);
    ScalarHaWp<D, MultiIndex> packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    //const dim_t order = 4;
    //using QR = GaussHermiteQR<order>;
    //using TQR = TensorProductQR<QR, QR, QR>;
    using TQR = TensorProductQR<GaussHermiteQR<3>, GaussHermiteQR<4>, GaussHermiteQR<5>>;
    TQR::NodeMatrix nodes;
    TQR::WeightVector weights;
    std::tie(nodes, weights) = TQR::nodes_and_weights();
    std::cout << "number of nodes: " << TQR::number_nodes() << "\n";
    //std::cout << "node matrix:\n" << nodes << "\n";
    //std::cout << "weight vector:\n" << weights << "\n";

    // Calculate inner product matrix, print it.
    HomogeneousInnerProduct<D, MultiIndex, TQR> ip;
    CMatrix<Eigen::Dynamic, Eigen::Dynamic> mat =
        ip.build_matrix(packet);

    std::cout << "IP matrix diagonal:\n";
    for(int i = 0; i < mat.rows(); ++i)
    {
        std::cout << " " << mat(i, i);
        if(i%3 == 2) std::cout << "\n";
    }
    std::cout << "\n";

    // Calculate quadrature.
    std::cout << "Quadrature: " << ip.quadrature(packet) << "\n";
}

void test1DInhomog()
{
    std::cout << "\n1D Inhomogeneous Gauss-Hermite Test\n";
    std::cout <<   "-----------------------------------\n";

    const real_t eps = 0.2;
    const dim_t D = 1;
    const dim_t N1 = 10, N2 = 12;
    const dim_t order = 8;
    using MultiIndex = TinyMultiIndex<unsigned short, D>;
    using QR = GaussHermiteQR<order>;

    // Set up sample 1D wavepacket.
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum1 =
        enumerator.generate(HyperCubicShape<D>(N1));
    ShapeEnum<D, MultiIndex> shape_enum2 =
        enumerator.generate(HyperCubicShape<D>(N2));
    HaWpParamSet<D> param_set1, param_set2;
    std::vector<complex_t> coeffs1(N1, 1.0);
    std::vector<complex_t> coeffs2(N2, 1.4);

    ScalarHaWp<D, MultiIndex> packet1;
    packet1.eps() = eps;
    packet1.parameters() = param_set1;
    packet1.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum1);
    packet1.coefficients() = coeffs1;

    ScalarHaWp<D, MultiIndex> packet2;
    packet2.eps() = eps;
    packet2.parameters() = param_set2;
    packet2.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum2);
    packet2.coefficients() = coeffs2;

    // Calculate inner product matrix, print it.
    InhomogeneousInnerProduct<D, MultiIndex, QR> ip;
    CMatrix<Eigen::Dynamic, Eigen::Dynamic> mat =
        ip.build_matrix(packet1, packet2);

    std::cout << "IP matrix:\n" << mat << std::endl;

    // Calculate quadrature.
    std::cout << "Quadrature: " << ip.quadrature(packet1, packet2) << "\n";
}

int main()
{
    //test1DGaussHermite();
    //test1DGaussHermiteOperator();
    //test3DGaussHermite();
    test1DInhomog();

    return 0;
}
