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
#include "waveblocks/util/timer.hpp"

using namespace waveblocks;

void run1D()
{
    const real_t eps = 0.2;
    const dim_t D = 1;
    const dim_t order = 8;
    using MultiIndex = TinyMultiIndex<unsigned short, D>;
    using QR = GaussHermiteQR<order>;
    using IP = HomogeneousInnerProduct<D, MultiIndex, QR>;
    const int n_runs = 1000;

    std::cout << "\n1D homogeneous quadrature of order " << order << ".\n";
    std::cout << "number of coefficients vs. evaluation time [ms] " <<
        "averaged over " << n_runs << " runs:\n";

    // Do runs for differently sized wavepackets.
    for (dim_t n_coeffs = 8 * order; n_coeffs <= 256 * order; n_coeffs *= 2)
    {
        // Set up sample 1D wavepacket.
        ShapeEnumerator<D, MultiIndex> enumerator;
        ShapeEnum<D, MultiIndex> shape_enum =
            enumerator.generate(HyperCubicShape<D>(n_coeffs));
        HaWpParamSet<D> param_set;
        Coefficients coeffs = Coefficients::Ones(n_coeffs, 1);

        ScalarHaWp<D, MultiIndex> packet;
        packet.eps() = eps;
        packet.parameters() = param_set;
        packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
        packet.coefficients() = coeffs;

        // Time many quadrature calculations.
        Timer timer;
        timer.start();
        for (int i = 0; i < n_runs; ++i)
        {
            // Calculate quadrature, volatile to force evaluation.
            volatile complex_t result = IP::quadrature(packet);
            (void)result;
        }
        timer.stop();

        std::cout << n_coeffs << "\t" << (timer.millis() / n_runs) << "\n";
    }
}

// Helper struct to avoid code duplication for different dimensionalities.
template<dim_t D, class TQR>
struct MultiDHelper
{
    static void run(real_t eps, dim_t n_coeffs, int n_runs)
    {
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

        // Time many quadrature calculations.
        Timer timer;
        timer.start();
        for (int i = 0; i < n_runs; ++i)
        {
            // Calculate quadrature, volatile to force evaluation.
            volatile complex_t result = IP::quadrature(packet);
            (void)result;
        }
        timer.stop();

        std::cout << D << "\t" << (timer.millis() / n_runs) << "\n";
    }
};

void runMultiD()
{
    const real_t eps = 0.2;
    const dim_t n_coeffs = 10;
    const dim_t order = 8;
    using QR = GaussHermiteQR<order>;

    std::cout << "\nHomogeneous quadrature of order " <<
        order << " with " << n_coeffs << " coefficients per dimension.\n";
    std::cout << "number of dimensions vs. evaluation time [ms]:\n";

    // Do runs for different numbers of dimensions.
    MultiDHelper<1, TensorProductQR<QR>>::run(eps, n_coeffs, 50000);
    MultiDHelper<2, TensorProductQR<QR,QR>>::run(eps, n_coeffs, 5000);
    MultiDHelper<3, TensorProductQR<QR,QR,QR>>::run(eps, n_coeffs, 100);
    MultiDHelper<4, TensorProductQR<QR,QR,QR,QR>>::run(eps, n_coeffs, 2);
}

template<dim_t LEVEL>
void runGenzKeister()
{
    const real_t eps = 0.2;
    const dim_t n_coeffs = 10;

    std::cout << "\nHomogeneous Genz-Keister quadrature of level " <<
        LEVEL << " with " << n_coeffs << " coefficients per dimension.\n";
    std::cout << "number of dimensions vs. evaluation time [ms]:\n";

    // Do runs for different numbers of dimensions.
    MultiDHelper<1, GenzKeisterQR<1, LEVEL>>::run(eps, n_coeffs, 50000);
    MultiDHelper<2, GenzKeisterQR<2, LEVEL>>::run(eps, n_coeffs, 5000);
    MultiDHelper<3, GenzKeisterQR<3, LEVEL>>::run(eps, n_coeffs, 100);
    MultiDHelper<4, GenzKeisterQR<4, LEVEL>>::run(eps, n_coeffs, 4);
}

void runMultiD2()
{
    const real_t eps = 0.2;
    const dim_t n_coeffs = 3;
    const dim_t order = 3;
    using QR = GaussHermiteQR<order>;

    std::cout << "\nHomogeneous quadrature of order " <<
        order << " with " << n_coeffs << " coefficients per dimension.\n";
    std::cout << "number of dimensions vs. evaluation time [ms]:\n";

    // Do runs for different numbers of dimensions.
    MultiDHelper<1, TensorProductQR<QR>>::run(eps, n_coeffs, 1000000);
    MultiDHelper<2, TensorProductQR<QR,QR>>::run(eps, n_coeffs, 500000);
    MultiDHelper<3, TensorProductQR<QR,QR,QR>>::run(eps, n_coeffs, 200000);
    MultiDHelper<4, TensorProductQR<QR,QR,QR,QR>>::run(eps, n_coeffs, 20000);
    MultiDHelper<5, TensorProductQR<QR,QR,QR,QR,QR>>::run(eps, n_coeffs, 1500);
    MultiDHelper<6, TensorProductQR<QR,QR,QR,QR,QR,QR>>::run(eps, n_coeffs, 100);
    MultiDHelper<7, TensorProductQR<QR,QR,QR,QR,QR,QR,QR>>::run(eps, n_coeffs, 4);
    MultiDHelper<8, TensorProductQR<QR,QR,QR,QR,QR,QR,QR,QR>>::run(eps, n_coeffs, 1);
}

template<dim_t LEVEL>
void runGenzKeister2()
{
    const real_t eps = 0.2;
    const dim_t n_coeffs = 3;

    std::cout << "\nHomogeneous Genz-Keister quadrature of level " <<
        LEVEL << " with " << n_coeffs << " coefficients per dimension.\n";
    std::cout << "number of dimensions vs. evaluation time [ms]:\n";

    // Do runs for different numbers of dimensions.
    MultiDHelper<1, GenzKeisterQR<1, LEVEL>>::run(eps, n_coeffs, 1000000);
    MultiDHelper<2, GenzKeisterQR<2, LEVEL>>::run(eps, n_coeffs, 500000);
    MultiDHelper<3, GenzKeisterQR<3, LEVEL>>::run(eps, n_coeffs, 200000);
    MultiDHelper<4, GenzKeisterQR<4, LEVEL>>::run(eps, n_coeffs, 50000);
    MultiDHelper<5, GenzKeisterQR<5, LEVEL>>::run(eps, n_coeffs, 8000);
    MultiDHelper<6, GenzKeisterQR<6, LEVEL>>::run(eps, n_coeffs, 1000);
    MultiDHelper<7, GenzKeisterQR<7, LEVEL>>::run(eps, n_coeffs, 50);
    MultiDHelper<8, GenzKeisterQR<8, LEVEL>>::run(eps, n_coeffs, 5);
    MultiDHelper<9, GenzKeisterQR<9, LEVEL>>::run(eps, n_coeffs, 1);
}

void runMultiComponent()
{
    const real_t eps = 0.2;
    const dim_t D = 2;
    const dim_t order = 8;
    const dim_t n_coeffs = 10;
    const int n_runs = 200;
    using MultiIndex = TinyMultiIndex<unsigned short, D>;
    using QR = GaussHermiteQR<order>;
    using TQR = TensorProductQR<QR,QR>;
    using IP = VectorInnerProduct<D, MultiIndex, TQR>;


    std::cout << "\nMulti-component " << D <<
        "D homogeneous quadrature of order " <<
        order << " with " << n_coeffs << " coefficients per dimension.\n";
    std::cout << "number of components vs. evaluation time [ms] " <<
        "averaged over " << n_runs << " runs:\n";

    for (dim_t n_components = 1; n_components <= 8; ++n_components)
    {
        // Set up sample 1D wavepacket.
        ShapeEnumerator<D, MultiIndex> enumerator;
        HaWpParamSet<D> param_set;

        HomogeneousHaWp<D, MultiIndex> packet(n_components);
        packet.eps() = eps;
        packet.parameters() = param_set;
        for (int comp = 0; comp < n_components; ++comp)
        {
            packet.component(comp).shape() =
                std::make_shared<ShapeEnum<D,MultiIndex>>(
                        enumerator.generate(HyperCubicShape<D>(n_coeffs)));
            packet.component(comp).coefficients() = Coefficients::Constant(
                    std::pow(n_coeffs,D), 1, 1);
        }

        // Time many quadrature calculations.
        Timer timer;
        timer.start();
        for (int i = 0; i < n_runs; ++i)
        {
            // Calculate quadrature, volatile to force evaluation.
            volatile auto result = IP::quadrature(packet);
            (void)result;
        }
        timer.stop();

        std::cout << n_components << "\t" << (timer.millis() / n_runs) << "\n";
    }
}

int main()
{
    run1D();
    runMultiD();
    runGenzKeister<6>();
    runMultiComponent();
    runMultiD2();
    runGenzKeister2<3>();

    return 0;
}
