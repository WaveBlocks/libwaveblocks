#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>

#include "basic_types.hpp"
#include "lexical_shape_enumerator.hpp"

namespace waveblocks {

template<std::size_t D>
struct HagedornParameters
{
public:
    typedef Eigen::Matrix<real_t,D,D> RMatrix;
    typedef Eigen::Matrix<complex_t,D,D> CMatrix;
    
    typedef Eigen::Matrix<real_t,D,1> RVector;
    typedef Eigen::Matrix<complex_t,D,1> CVector;
    
    real_t eps;
    RVector q, p;
    CMatrix Q, P;
    
    HagedornParameters()
        : eps(1.0)
        , q(RVector::Zero())
        , p(RVector::Zero())
        , Q(CMatrix::Identity())
        , P(CMatrix::Identity()*complex_t(0.0, 1.0))
    { }
    
    HagedornParameters(real_t eps, const RVector &q, const RVector &p, const CMatrix &Q, const CMatrix &P)
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
    
    /**
     * 
     */
    complex_t evaluateBasis(const RVector &x, 
                            std::size_t axis,
                            MultiIndex<D> k, 
                            complex_t cur_v, 
                            const CVector &prev_v) const
    {
        CMatrix Qinv = Q.inverse();
        CMatrix QhQinvt = Q.adjoint()*Qinv.transpose();
        
        //compute {sqrt(k[i])*phi[k-e[i]]}
        //  e[i]: unit vector aligned to i-th axis
        CVector prev_v_scaled = prev_v;
        for (std::size_t i = 0; i < D; i++)
            prev_v_scaled(i,0) *= std::sqrt(k[i]);
        
        complex_t pr1 = std::sqrt(2.0)/eps * Qinv.row(axis).dot(x-q)*cur_v;
        complex_t pr2 = QhQinvt.row(axis).dot(prev_v_scaled);
        
        return (pr1 - pr2) / std::sqrt(k[axis]+1);
    }
};

template<std::size_t D, class S>
class HagedornWavepacket
{
private:
    LexicalShapeEnumeration<D,S> enumeration_;
    std::vector<complex_t> coefficients_;
    
public:
    HagedornParameters<D> parameters;
    
    HagedornWavepacket(LexicalShapeEnumeration<D,S> enumeration)
        : enumeration_(enumeration)
        , coefficients_(enumeration.size())
        , parameters()
    { }
    
    HagedornWavepacket(LexicalShapeEnumeration<D,S> enumeration, const HagedornParameters<D> &parameters)
        : enumeration_(enumeration)
        , coefficients_(enumeration.size())
        , parameters(parameters)
    { }
    
    /**
     * stores function values of all bases in vector
     */
    void evaluateBases(const Eigen::Matrix<real_t,D,1> &x, std::vector<complex_t> &result) const
    {
        result.resize(enumeration_.size());
        
        LexicalShapeIterator<D,S> it = enumeration_.begin();
        
        //first entry in iterator must be 0
        assert( it.getMultiIndex() == MultiIndex<D>{} );
        
        complex_t phi0 = parameters.evaluateGroundState(x);
        
        while (it.advance()) {
            MultiIndex<D> index = it.getMultiIndex();
            
            std::size_t axis = D;
            
            //find first non-zero axis
            for (std::size_t d = 0; d < D; d++) {
                if (index[d] != 0) {
                    axis = d;
                    break;
                }
            }
            
            assert(axis != D);
            
            //receive precursors
            LexicalShapeIterator<D,S> preit = enumeration_.getBackwardNeighbour(it, axis);
            
            assert(preit.getOrdinal() < it.getOrdinal());
            
            complex_t value = result[preit.getOrdinal()];
            
            Eigen::Matrix<complex_t,D,1> values2;
            for (std::size_t d = 0; d < D; d++) {
                if (preit.getMultiIndex()[d] == 0) {
                    values2(d,0) = complex_t();
                } else {
                    LexicalShapeIterator<D,S> pre2it = enumeration_.getBackwardNeighbour(preit, d);
                    
                    assert(pre2it.getOrdinal() < preit.getOrdinal());
                    
                    values2(d,0) = result[pre2it.getOrdinal()];
                }
            }
            
            result[it.getOrdinal()] = 
                parameters.evaluateBasis(x, axis, preit.getMultiIndex(), value, values2);
        }
    }
};

}

#endif