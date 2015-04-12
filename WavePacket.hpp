#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET_HPP
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET_HPP

#include <Eigen/Core>

// includes matrix inverse
#include <Eigen/LU>

#include <unordered_set>

namespace waveblocks {

template<dim_t D>
struct HagedornParameters
{
    Eigen::Matrix<real_t,D,1> q, p;
    Eigen::Matrix<complex_t,D,D> Q, P;
    
    real_t eps;
    
    HagedornParameters() : q(), p(), 
        Q(Eigen::Matrix<complex_t,D,D>::Identity()), 
        P(Eigen::Matrix<complex_t,D,D>::Identity().conjugate()),
        eps(1.0) {}
    
    HagedornParameters(const HagedornParameters &that) : 
        q(that.q), p(that.p), Q(that.Q), P(that.P), eps(that.eps) {}
};

template<dim_t D>
class HagedornWavepacket
{
public:
    HagedornParameters<D> parameters;
    //E enumeration;
    
    /**
     * coefficients->size() >= enumeration.size()
     */
    std::vector<complex_t> coefficients;
    
    Shape::Enumeration &enumeration_;
    
    HagedornWavepacket() : parameters(), coefficients() {}
    
    HagedornWavepacket(const HagedornWavepacket &that) :
            parameters(that.parameters),
            coefficients(that.coefficients) {}
    
    complex_t evaluateGroundState(const Eigen::Matrix<real_t,D,1> &x) const
    {
        Eigen::Matrix<real_t,D,1> dx = x - parameters.q;
        
        complex_t pr1 = dx.transpose()*parameters.P*parameters.Q.inverse()*dx;
        real_t pr2 = parameters.p.transpose()*dx;
        
        complex_t exponent = complex_t(0.0, 1.0)/(parameters.eps*parameters.eps) * (0.5*pr1 + pr2);
        
        return std::exp(exponent);
    }
    
    complex_t evaluateGroundState(const Eigen::Matrix<real_t,D,1> &x, ContinousSqrt &csqrt) const
    {
        return evaluateGroundState(x) / csqrt(parameters.Q.determinant());
    }
    
    /**
     * evalutes wave function value without storing every bases function value
     */
    complex_t evaluateWaveFunction(const Eigen::Matrix<real_t,D,1> &x) const
    {
        
    }
    
    /**
     * stores function values of all bases in vector
     */
    void evaluateBases(const Eigen::Matrix<real_t,D,1> &x, std::vector<complex_t> &result) const
    {
        Shape<D>::Enumeration::Iterator it = enumeration_.first();
        
        //first entry in iterator must be 0
        assert( it.getMultiIndex() == MultiIndex<D>() );
        
        while (it.hasNext())
        {
            it = it.next();
            
            MultiIndex<D> index = it.getMultiIndex();
            
            dim_t axis = -1;
            
            //find first non-zero axis
            for (dim_t d = 0; d < D; d++) {
                if (index[d] != 0) {
                    axis = d;
                    break;
                }
            }
            
            if (axis == -1)
                throw std::runtime_error("enumerator returns zero entry");
            
            //receive precursors
            Shape<D>::Enumeration::Iterator precursor = it.getBackwardNeighbour(axis);
            
            //check weak ordering fulfilment
            assert(precursor.getOrdinal() < it.getOrdinal());
            
            complex_t value = result[precursor.getOrdinal()];
            
            std::array<complex_t, D> values2;
            for (dim_t d = 0; d < D; d++) {
                if (precursor.getMultiIndex()[d] == 0) {
                    values2[d] = 0;
                } else {
                    Shape<D>::Enumeration::Iterator precursor2 = precursor.getBackwardNeighbour(d);
                    
                    //check weak ordering fulfilment
                    assert(precursor2.getOrdinal() < precursor.getOrdinal());
                    
                    values2[d] = result[precursor2.getOrdinal()];
                }
            }
            
            result[it.getOrdinal()] = evaluateBasis(axis, value, values2);
        }
    }
    
private:
    complex_t evaluateBasis(dim_t axis, complex_t current, 
                              const std::array<complex_t, D> &precursors) const
    {
        
    }
};

}
    
#endif