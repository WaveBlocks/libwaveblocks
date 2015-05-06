#ifndef WAVEBLOCKS_HAGEDORN_PARAMETER_SET
#define WAVEBLOCKS_HAGEDORN_PARAMETER_SET

#include <cmath>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "basic_types.hpp"
#include "multi_index.hpp"


namespace waveblocks {

template<dim_t D>
struct HagedornParameterSet
{
    typedef Eigen::Matrix<real_t,D,D> RMatrix;
    typedef Eigen::Matrix<complex_t,D,D> CMatrix;
    
    typedef Eigen::Matrix<real_t,D,1> RVector;
    typedef Eigen::Matrix<complex_t,D,1> CVector;
    
    real_t eps;
    RVector q, p;
    CMatrix Q, P;
    
    HagedornParameterSet()
        : eps(1.0)
        , q(RVector::Zero())
        , p(RVector::Zero())
        , Q(CMatrix::Identity())
        , P(CMatrix::Identity()*complex_t(0.0, 1.0))
    { }
    
    HagedornParameterSet(real_t eps, const RVector &q, const RVector &p, const CMatrix &Q, const CMatrix &P)
        : eps(eps)
        , q(q)
        , p(p)
        , Q(Q)
        , P(P)
    { }
    
    complex_t evaluateGroundState(const Eigen::Matrix<real_t,D,1> &x) const
    {
        const real_t pi = 3.14159265359;
        
        Eigen::Matrix<real_t,D,1> dx = x - q;
        
        complex_t pr1 = dx.transpose()*P*Q.inverse()*dx;
        real_t pr2 = p.transpose()*dx;
        
        complex_t exponent = complex_t(0.0, 1.0)/(eps*eps) * (0.5*pr1 + pr2);
        
        return 1.0/std::pow(pi*eps*eps, D/4.0) * std::exp(exponent);
    }
};

template<dim_t D>
std::ostream &operator<<(std::ostream &out, const HagedornParameterSet<D> &parameters)
{
    Eigen::IOFormat matrix_format(Eigen::StreamPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
    
    out << "--------------------\n";
    out << "HagedornParameterSet\n";
    out << "#eps: " << parameters.eps << '\n';
    out << "#q: \n" << parameters.q.format(matrix_format) << '\n';
    out << "#p: \n" << parameters.p.format(matrix_format) << '\n';
    out << "#Q: \n" << parameters.Q.format(matrix_format) << '\n';
    out << "#P: \n" << parameters.P.format(matrix_format) << '\n';
    out << "--------------------\n";
    return out;
}

}

#endif