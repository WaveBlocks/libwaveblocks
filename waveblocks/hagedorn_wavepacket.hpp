#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>

#include "basic_types.hpp"
#include "lexical_shape_enumerator.hpp"

namespace waveblocks {

template<dim_t D, class S>
class HagedornWavepacket
{
private:
    LexicalShapeEnumeration<D,S> enumeration_;
    std::vector<complex_t> coefficients_;
    
public:
    HagedornParameterSet<D> parameters;
    
    HagedornWavepacket(LexicalShapeEnumeration<D,S> enumeration)
        : enumeration_(enumeration)
        , coefficients_(enumeration.size())
        , parameters()
    { }
    
    HagedornWavepacket(LexicalShapeEnumeration<D,S> enumeration, const HagedornParameterSet<D> &parameters)
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
            
            dim_t axis = D;
            
            //find first non-zero axis
            for (dim_t d = 0; d < D; d++) {
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
            for (dim_t d = 0; d < D; d++) {
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