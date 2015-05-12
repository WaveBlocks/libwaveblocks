#ifndef WAVEBLOCKS_HAGEDORN_PARAMETER_SET
#define WAVEBLOCKS_HAGEDORN_PARAMETER_SET

#include <cmath>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "basic_types.hpp"
#include "multi_index.hpp"
#include "continuous_sqrt.hpp"


namespace waveblocks {

template<dim_t D>
struct HagedornParameterSet
{
    typedef Eigen::Matrix<real_t,D,D> RMatrix;
    typedef Eigen::Matrix<complex_t,D,D> CMatrix;
    
    typedef Eigen::Matrix<real_t,D,1> RVector;
    typedef Eigen::Matrix<complex_t,D,1> CVector;
    
    RVector q, p;
    CMatrix Q, P;
    ContinuousSqrt<real_t> sqrt_detQ;
    
    HagedornParameterSet()
        : q(RVector::Zero())
        , p(RVector::Zero())
        , Q(CMatrix::Identity())
        , P(CMatrix::Identity()*complex_t(0.0, 1.0))
        , sqrt_detQ() //detQ = 1.0 => sqrt(detQ) = 1.0
    { }
    
    HagedornParameterSet(const HagedornParameterSet &that)
        : q(that.q)
        , p(that.p)
        , Q(that.Q)
        , P(that.P)
        , sqrt_detQ(that.sqrt_detQ)
    { }
    
    HagedornParameterSet(const RVector &q, const RVector &p, const CMatrix &Q, const CMatrix &P)
        : q(q)
        , p(p)
        , Q(Q)
        , P(P)
        , sqrt_detQ( std::sqrt(Q.determinant()) ) //choose 1st root
    { }
    
    HagedornParameterSet(const RVector &q, const RVector &p, const CMatrix &Q, const CMatrix &P, ContinuousSqrt<real_t> sqrt_detQ)
        : q(q)
        , p(p)
        , Q(Q)
        , P(P)
        , sqrt_detQ(sqrt_detQ)
    { }
    
    HagedornParameterSet &operator=(const HagedornParameterSet &that)
    {
        q = that.q;
        p = that.p;
        Q = that.Q;
        P = that.P;
        sqrt_detQ = that.sqrt_detQ;
        return *this;
    }
};

template<dim_t D>
std::ostream &operator<<(std::ostream &out, const HagedornParameterSet<D> &parameters)
{
    Eigen::IOFormat matrix_format(Eigen::StreamPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
    
    out << "--------------------\n";
    out << "HagedornParameterSet\n";
    out << "#q: \n" << parameters.q.format(matrix_format) << '\n';
    out << "#p: \n" << parameters.p.format(matrix_format) << '\n';
    out << "#Q: \n" << parameters.Q.format(matrix_format) << '\n';
    out << "#P: \n" << parameters.P.format(matrix_format) << '\n';
    out << "#sqrt(detQ): " << std::abs(parameters.sqrt_detQ()) << "*exp(" << std::arg(parameters.sqrt_detQ()) << "*i)" << '\n';
    out << "--------------------\n";
    return out;
}

}

#endif