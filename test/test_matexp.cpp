#include <iostream>

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

typedef double real_t;
typedef std::complex<real_t> complex_t;

template<int R, int C>
using RMatrix = Eigen::Matrix<complex_t,R,C>;

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;

    RMatrix<2,2> M;
    M << 1, -1, -1, 1;

    // M = [[ 1, -1],
    //      [-1,  1]]

    // exp(A) = [[ 4.19452805, -3.19452805],
    //           [-3.19452805,  4.19452805]]

    std::cout << M << std::endl;
    std::cout << M.exp() << std::endl;
}
