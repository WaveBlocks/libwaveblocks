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
struct HaWpParamSet
{
    RMatrix<D,1> q, p;
    CMatrix<D,D> Q, P;
    complex_t S;
    ContinuousSqrt<real_t> sqrt_detQ;
    
    HaWpParamSet()
        : q(RMatrix<D,1>::Zero())
        , p(RMatrix<D,1>::Zero())
        , Q(CMatrix<D,D>::Identity())
        , P(CMatrix<D,D>::Identity()*complex_t(0.0, 1.0))
        , S(complex_t(0.0, 0.0))
        , sqrt_detQ() //detQ = 1.0 => sqrt(detQ) = 1.0
    { }
    
    HaWpParamSet(const HaWpParamSet &that)
        : q(that.q)
        , p(that.p)
        , Q(that.Q)
        , P(that.P)
        , S(that.S)
        , sqrt_detQ(that.sqrt_detQ)
    { }

    HaWpParamSet(const RMatrix<D,1> &q,
                 const RMatrix<D,1> &p,
                 const CMatrix<D,D> &Q,
                 const CMatrix<D,D> &P,
                 const complex_t &S)
        : q(q)
        , p(p)
        , Q(Q)
        , P(P)
        , S(S)
        , sqrt_detQ( std::sqrt(Q.determinant()) ) //choose 1st root
    { }

    HaWpParamSet(const RMatrix<D,1> &q,
                 const RMatrix<D,1> &p,
                 const CMatrix<D,D> &Q,
                 const CMatrix<D,D> &P,
                 const complex_t &S,
                 ContinuousSqrt<real_t> sqrt_detQ)
        : q(q)
        , p(p)
        , Q(Q)
        , P(P)
        , sqrt_detQ(sqrt_detQ)
        , S(S)
    { }
    
    HaWpParamSet &operator=(const HaWpParamSet &that)
    {
        q = that.q;
        p = that.p;
        Q = that.Q;
        P = that.P;
        S = that.S;
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
    std::pair< RMatrix<D,1>, RMatrix<D,D> > mix(const HaWpParamSet<D>& other) const
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
std::ostream &operator<<(std::ostream &out, const HaWpParamSet<D> &parameters)
{
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n     ", "[", "]");
    
    out << "HaWpParamSet {\n";
    out << "  q: " << parameters.q.format(CleanFmt) << '\n';
    out << "  p: " << parameters.p.format(CleanFmt) << '\n';
    out << "  Q: " << parameters.Q.format(CleanFmt) << '\n';
    out << "  P: " << parameters.P.format(CleanFmt) << '\n';
    out << "  S: " << parameters.S << '\n';
    out << "  sqrt(detQ): " << std::abs(parameters.sqrt_detQ()) << "*exp(" << std::arg(parameters.sqrt_detQ()) << "*i)" << '\n';
    out << "  compatible(): " << (parameters.compatible() ? "yes" : "no") << '\n';
    out << "}" << std::endl;
    return out;
}

}

#endif
