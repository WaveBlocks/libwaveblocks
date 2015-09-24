#include <iostream>
#include <Eigen/Core>
#include <waveblocks/hawp_paramset.hpp>


using namespace waveblocks;


int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;

    /*
    // Test 1
    const dim_t D = 1;
    HaWpParamSet<D> PIbra = HaWpParamSet<D>();
    HaWpParamSet<D> PIket = HaWpParamSet<D>();

    real_t t = 0.9;
    real_t u = 1.2;
    real_t a = std::sin(t);
    real_t b = std::cos(t);
    real_t c = std::sin(u);
    real_t d = std::cos(u);

    RMatrix<D,1> qb;
    qb << 0.8;
    PIbra.q(qb);
    CMatrix<D,D> Pb;
    Pb << complex_t(-a,b);
    CMatrix<D,D> Qb;
    Qb << complex_t(b,a);
    PIbra.Q(Qb);
    PIbra.P(Pb);

    RMatrix<D,1> qk;
    qk << 1.1;
    PIket.q(qk);
    CMatrix<D,D> Pk;
    Pk << complex_t(-c,d);
    CMatrix<D,D> Qk;
    Qk << complex_t(d,c);
    PIket.Q(Qk);
    PIket.P(Pk);
    */


    // Test 2
    const dim_t D = 2;
    HaWpParamSet<D> PIbra = HaWpParamSet<D>();
    HaWpParamSet<D> PIket = HaWpParamSet<D>();

    RMatrix<D,1> qb;
    qb <<  0.1, -0.2;
    PIbra.q(qb);
    PIbra.Q(PIbra.Q() / 2);
    PIbra.P(PIbra.P() * 2);

    RMatrix<D,1> qk;
    qk << 0.5, 1.0;
    PIket.q(qk);
    PIket.Q(PIket.Q() * 3);
    PIket.P(PIket.P() / 3);


    // Mix it
    auto PImix = PIbra.mix(PIket);

    Eigen::IOFormat CleanFmt(16, 0, ", ", "\n   ", "[", "]");

    std::cout << "Parameter set from bra part:" << std::endl;
    std::cout << "q: " << PIbra.q().format(CleanFmt) << std::endl;
    std::cout << "p: " << PIbra.p().format(CleanFmt) << std::endl;
    std::cout << "Q: " << PIbra.Q().format(CleanFmt) << std::endl;
    std::cout << "P: " << PIbra.P().format(CleanFmt) << std::endl;
    std::cout << "compatible: " << (PIbra.compatible() ? "yes" : "no") << std::endl;

    std::cout << "Parameter set from ket part:" << std::endl;
    std::cout << "q: " << PIket.q().format(CleanFmt) << std::endl;
    std::cout << "p: " << PIket.p().format(CleanFmt) << std::endl;
    std::cout << "Q: " << PIket.Q().format(CleanFmt) << std::endl;
    std::cout << "P: " << PIket.P().format(CleanFmt) << std::endl;
    std::cout << "compatible: " << (PIket.compatible() ? "yes" : "no") << std::endl;

    std::cout << "Mixed Parameters q and Q:" << std::endl;
    std::cout << "q: " << std::get<0>(PImix).format(CleanFmt) << std::endl;
    std::cout << "Q: " << std::get<1>(PImix).format(CleanFmt) << std::endl;

    std::cout << "Done" << std::endl;
}
