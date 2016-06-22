#include <iostream>
#include <fstream>

#include "waveblocks/types.hpp"
#include "waveblocks/wavepackets/shapes/tiny_multi_index.hpp"
#include "waveblocks/potentials/potentials.hpp"
#include "waveblocks/potentials/bases.hpp"
#include "waveblocks/wavepackets/hawp_paramset.hpp"
#include "waveblocks/wavepackets/hawp_commons.hpp"
#include "waveblocks/wavepackets/shapes/shape_enumerator.hpp"
#include "waveblocks/wavepackets/shapes/shape_hypercubic.hpp"
#include "waveblocks/innerproducts/gauss_hermite_qr.hpp"
#include "waveblocks/innerproducts/tensor_product_qr.hpp"
#include "waveblocks/propagators/Hagedorn.hpp"
#include "waveblocks/observables/energy.hpp"
#include "waveblocks/io/hdf5writer.hpp"


using namespace waveblocks;

int main() {
    const int N = 1;
    const int D = 2;
    const int K = 4;

    const real_t sigma_x = 0.5;
    const real_t sigma_y = 0.5;

    const real_t tol = 1e-14;

    const real_t T = 12;
    const real_t dt = 0.01;

    const real_t eps = 0.1;

    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned long, D>;

    // The parameter set of the initial wavepacket
    CMatrix<D,D> Q = CMatrix<D,D>::Identity();
    CMatrix<D,D> P = complex_t(0,1) * CMatrix<D,D>::Identity();
    RVector<D> q = {-3.0, 0.0};
    RVector<D> p = { 0.0, 0.5};
    complex_t S = 0.0;
    wavepackets::HaWpParamSet<D> param_set(q,p,Q,P,S);

    // Basis shape
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K));

    // Gaussian Wavepacket phi_00 with c_00 = 1
    Coefficients coeffs = Coefficients::Zero(std::pow(K, D), 1);
    coeffs[0] = 1.0;
    Coefficients coefforig = Coefficients(coeffs);

    // Assemble packet
    wavepackets::ScalarHaWp<D,MultiIndex> packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    // Defining the potential
    typename CanonicalBasis<N,D>::potential_type potential =
        [sigma_x,sigma_y](CVector<D> x) {
        return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
    };
    typename ScalarLeadingLevel<D>::potential_type leading_level = potential;
    typename ScalarLeadingLevel<D>::jacobian_type leading_jac =
        [sigma_x,sigma_y](CVector<D> x) {
        return CVector<D>{sigma_x*x[0], sigma_y*x[1]};
    };
    typename ScalarLeadingLevel<D>::hessian_type leading_hess =
        [sigma_x,sigma_y](CVector<D> /*x*/) {
        CMatrix<D,D> res;
        res(0,0) = sigma_x;
        res(1,1) = sigma_y;
        return res;
    };

    ScalarMatrixPotential<D> V(potential,leading_level,leading_jac,leading_hess);

    // Quadrature rules
    using TQR = innerproducts::TensorProductQR<innerproducts::GaussHermiteQR<K+4>,
                                               innerproducts::GaussHermiteQR<K+4>>;

    // Defining the propagator
    propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;


    io::hdf5writer<D> mywriter2("harmonic_2D_cpp.hdf5");
    mywriter2.set_write_energy(true);
    mywriter2.prestructuring<MultiIndex>(packet,dt);

    //write time = 0
    real_t ekin = observables::kinetic_energy<D,MultiIndex>(packet);
    real_t epot = observables::potential_energy<ScalarMatrixPotential<D>,D,MultiIndex,TQR>(packet,V);
    mywriter2.store_packet(packet);
    mywriter2.store_norms(packet);
    mywriter2.store_energies(epot,ekin);


    std::cout << "Time: " << 0 << std::endl;
    std::cout << packet.parameters() << std::endl;

    // Propagation
    for (real_t t = dt; t < T; t += dt) {
        std::cout << "Time: " << t << std::endl;

        // Propagate
        propagator.propagate(packet,dt,V);
        std::cout << packet.parameters() << std::endl;

        real_t ekin = observables::kinetic_energy<D,MultiIndex>(packet);
        real_t epot = observables::potential_energy<ScalarMatrixPotential<D>,D,MultiIndex,TQR>(packet,V);

        mywriter2.store_packet(packet);
        mywriter2.store_energies(epot,ekin);
        mywriter2.store_norms(packet);

        std::cout << "E: (p,k,t) " << epot << ", " << ekin << ", " << ekin+epot << std::endl;

        // Assure constant coefficients
        auto diff = (packet.coefficients() - coefforig).array().abs();
        auto norm = diff.matrix().template lpNorm<Eigen::Infinity>();
        bool flag = norm > tol ? false : true;
        std::cout << "Coefficients constant? " << (flag ? "yes" : "no") << std::endl;


    }
    mywriter2.poststructuring();
}
