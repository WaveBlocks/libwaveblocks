#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/hawp.hpp"
#include "waveblocks/hawp_paramset.hpp"
#include "waveblocks/homogeneous_inner_product.hpp"
#include "waveblocks/shape_enumerator.hpp"
#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/tiny_multi_index.hpp"

using namespace waveblocks;

int main()
{
    const dim_t D = 1;
    const dim_t N = 10;
    const dim_t order = 8;
    using MultiIndex = TinyMultiIndex<unsigned short, D>;

    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum =
        enumerator.generate(HyperCubicShape<D>(N));
    HaWpParamSet<D> param_set;
    std::cout << param_set << std::endl;
    std::vector<complex_t> coeffs(N, 0);

    HaWp<D, MultiIndex> packet(0.6, &param_set, &shape_enum, &coeffs);

    HomogeneousInnerProduct<D, MultiIndex> ip;
    CMatrix<Eigen::Dynamic, Eigen::Dynamic> mat =
        ip.build_matrix(packet, order);

    std::cout << "IP matrix:\n" << mat << std::endl;

    return 0;
}
