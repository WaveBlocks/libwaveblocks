#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/gauss_hermite_qr.hpp"
#include "waveblocks/hawp.hpp"
#include "waveblocks/hawp_paramset.hpp"
#include "waveblocks/homogeneous_inner_product.hpp"
#include "waveblocks/shape_enumerator.hpp"
#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/stdarray2stream.hpp"
#include "waveblocks/tiny_multi_index.hpp"

using namespace waveblocks;

int main()
{
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
    std::vector<complex_t> coeffs(N, 0);

    // Print QR nodes and weights.
    std::cout << "nodes: {";
    for (auto x : QR::nodes()) std::cout << " " << x;
    std::cout << " }\n";

    std::cout << "weights: {";
    for (auto x : QR::weights()) std::cout << " " << x;
    std::cout << " }\n";

    HaWp<D, MultiIndex> packet(eps, &param_set, &shape_enum, &coeffs);

    // Calculate inner product matrix, print it.
    HomogeneousInnerProduct<D, MultiIndex, QR> ip;
    CMatrix<Eigen::Dynamic, Eigen::Dynamic> mat =
        ip.build_matrix(packet);

    std::cout << "IP matrix:\n" << mat << std::endl;

    return 0;
}
