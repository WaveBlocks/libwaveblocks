#ifndef WAVEBLOCKS_HAGEDORN_PARAMETER_SET
#define WAVEBLOCKS_HAGEDORN_PARAMETER_SET

#include <cmath>
#include <iostream>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "basic_types.hpp"
#include "continuous_sqrt.hpp"


namespace waveblocks {

template<dim_t D>
struct HagedornParameterSet
{
    RMatrix<D,1> q, p;
    CMatrix<D,D> Q, P;
    ContinuousSqrt<real_t> sqrt_detQ;
    
    HagedornParameterSet()
        : q(RMatrix<D,1>::Zero())
        , p(RMatrix<D,1>::Zero())
        , Q(CMatrix<D,D>::Identity())
        , P(CMatrix<D,D>::Identity()*complex_t(0.0, 1.0))
        , sqrt_detQ() //detQ = 1.0 => sqrt(detQ) = 1.0
    { }
    
    HagedornParameterSet(const HagedornParameterSet &that)
        : q(that.q)
        , p(that.p)
        , Q(that.Q)
        , P(that.P)
        , sqrt_detQ(that.sqrt_detQ)
    { }
    
    HagedornParameterSet(const RMatrix<D,1> &q, const RMatrix<D,1> &p, const CMatrix<D,D> &Q, const CMatrix<D,D> &P)
        : q(q)
        , p(p)
        , Q(Q)
        , P(P)
        , sqrt_detQ( std::sqrt(Q.determinant()) ) //choose 1st root
    { }
    
    HagedornParameterSet(const RMatrix<D,1> &q, const RMatrix<D,1> &p, const CMatrix<D,D> &Q, const CMatrix<D,D> &P, ContinuousSqrt<real_t> sqrt_detQ)
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
    
    /**
     * Checks for compatibility relations
     * For details see master thesis 3.11
     */
    bool compatible() const
    {
        CMatrix<D,D> C1 = Q.adjoint()*P - P.adjoint()*Q - CMatrix<D,D>::Identity()*complex_t(0,2);
        CMatrix<D,D> C2 = P.transpose()*Q - Q.transpose()*P;
        
        double tol = 1e-10;
        
        return C1.template lpNorm<1>() < tol && C2.template lpNorm<1>() < tol;
    }
    
    /**
     * 
     */
    std::pair< RMatrix<D,1>, RMatrix<D,D> > mix(const HagedornParameterSet<D>& other) const
    {
        // Mix the parameters
        CMatrix<D,D> Gr = P * Q.inverse();
        CMatrix<D,D> Gc = other.P * other.Q.inverse();
        
        RMatrix<D,D> G = (Gc - Gr.adjoint()).imag();
        RMatrix<D,1> g = (Gc*other.q - Gr.adjoint()*q).imag();
        RMatrix<D,1> q0 = G.inverse() * g;
        RMatrix<D,D> Q0 = 0.5 * G;
        
        // We can not avoid the matrix root by using svd
        RMatrix<D,D> Qs = Q0.sqrt().inverse();
        
        // Assign (q0, Qs)
        return {q0,Qs};
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
    out << "#compatible?: " << (parameters.compatible() ? "yes" : "no") << '\n';
    out << "--------------------\n";
    out << std::endl; //flush
    return out;
}

}

#endif