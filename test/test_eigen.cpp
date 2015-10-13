#include <iostream>

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

typedef double real_t;
typedef std::complex<real_t> complex_t;

template<int R, int C>
using CMatrix = Eigen::Matrix<complex_t,R,C>;

// gcc flags
/*
  -std=c++11
  -Wall
  -Wno-deprecated-declarations
  -Wextra
  -pedantic
  -Ofast
*/

// Output with 'auto'
/*
(-1.22474,0)        (0,0)  (1.22474,0)
Qs: [(0.707107,0.707107)]
Qs: [(0.146029,-0.296525)]
(-0.0666025,-0.0866025)                (0.02,0)    (0.106603,0.0866025)
*/

// Output with CMatrix<D,D>
/*
(-1.22474,0)        (0,0)  (1.22474,0)
Qs: [(1.00005,0)]
Qs: [(1.00005,0)]
(-0.102481,0)      (0.02,0)  (0.142481,0)
*/

void foo() {
    const int D = 1;

    CMatrix<1,3> nodes;
    nodes << -1.22474487139, 0.0, 1.22474487139;
    const CMatrix<1,3> cnodes = complex_t(1,0) * nodes;

    CMatrix<D,1> q;
    q << 0.02;

    CMatrix<D,D> Q;
    Q << complex_t(1, 0.01);


    std::cout << cnodes << std::endl;


    // Compute affine transformation.
    //CMatrix<D,D> Qs = (Q * Q.adjoint()).sqrt();
    auto Qs = (Q * Q.adjoint()).sqrt();

    std::cout << "Qs: [" << Qs << "]" << std::endl;

    CMatrix<D,3> transformed_nodes = q.replicate(1, 3) + 0.1 * (Qs * cnodes);

    std::cout << "Qs: [" << Qs << "]" << std::endl;


    std::cout << transformed_nodes << std::endl;
}


int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;

    foo();
}
