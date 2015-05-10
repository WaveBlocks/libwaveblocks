#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <valarray>
#include <memory>
#include <initializer_list>

#include "basic_types.hpp"
#include "hagedorn_parameter_set.hpp"
#include "sliced_shape_enumeration.hpp"

#include "kahan_sum.hpp"

namespace waveblocks {

/**
 * \tparam D degrees of freedom
 * \tparam C number of components if this wavepacket is vector-valued
 * \parblock
 * Components of vector-valued wavepackets share the same paramaters and enumeration,
 * but use different coefficients.
 * \endparblock
 */
template<dim_t D, dim_t C = 1>
class HagedornWavepacket
{
private:
    real_t eps_;
    std::shared_ptr< HagedornParameterSet<D> > parameters_;
    std::shared_ptr< ShapeEnumeration<D> > enumeration_;
    std::array< std::shared_ptr< std::valarray<complex_t> >, C> coefficients_;
    
    complex_t evaluateGroundState(const Eigen::Matrix<real_t,D,1> &x) const
    {
        const real_t pi = 3.14159265359;
        
        Eigen::Matrix<real_t,D,1> dx = x - parameters_->q;
        
        complex_t pr1 = dx.transpose()*parameters_->P*parameters_->Q.inverse()*dx;
        real_t pr2 = parameters_->p.transpose()*dx;
        
        complex_t exponent = complex_t(0.0, 1.0)/(eps_*eps_) * (0.5*pr1 + pr2);
        
        return 1.0/std::pow(pi*eps_*eps_, D/4.0) * std::exp(exponent);
    }
    
    Eigen::Matrix<complex_t,C,1> evaluateBasis(const Eigen::Matrix<complex_t,D,D> &Qinv,
                                const Eigen::Matrix<complex_t,D,D> &QhQinvt,
                                dim_t axis, 
                                MultiIndex<D> k, 
                                const Eigen::Matrix<complex_t,C,1> &curr_basis, 
                                const Eigen::Matrix<complex_t,C,D> &prev_bases, 
                                const Eigen::Matrix<real_t,D,1> &x) const
    {
        //compute {sqrt(k[i])*phi[k-e[i]]}
        //  e[i]: unit vector aligned to i-th axis
        Eigen::Matrix<complex_t,D,C> prev_bases_scaled = prev_bases.transpose();
        for (dim_t d = 0; d < D; d++)
            prev_bases_scaled.row(d) *= std::sqrt( real_t(k[d]) );
        
        Eigen::Matrix<real_t,D,1> dx = x - parameters_->q;
        
        complex_t temp = Qinv.row(axis)*dx;
        
        Eigen::Matrix<complex_t,C,1> pr1 = curr_basis.transpose() * std::sqrt(2.0)/eps_ * temp;
        Eigen::Matrix<complex_t,C,1> pr2 = QhQinvt.row(axis)*prev_bases_scaled;
        
        return (pr1 - pr2) / std::sqrt( real_t(k[axis])+1.0);
    }
    
public:
    HagedornWavepacket(real_t eps,
                       std::shared_ptr< HagedornParameterSet<D> > parameters, 
                       std::shared_ptr< ShapeEnumeration<D> > enumeration,
                       std::initializer_list< std::shared_ptr< std::valarray<complex_t> > > coefficients)
        : eps_(eps)
        , parameters_(parameters)
        , enumeration_(enumeration)
        , coefficients_()
    { 
        dim_t c = 0;
        for (auto & item : coefficients)
            coefficients_[c++] = item;
        assert (c == C);
    }
    
    HagedornWavepacket(real_t eps,
                       std::shared_ptr< HagedornParameterSet<D> > parameters, 
                       std::shared_ptr< ShapeEnumeration<D> > enumeration,
                       std::array< std::shared_ptr< std::valarray<complex_t> >, C> coefficients)
        : eps_(eps)
        , parameters_(parameters)
        , enumeration_(enumeration)
        , coefficients_(coefficients)
    { }
    
    HagedornWavepacket(real_t eps,
                       std::shared_ptr< HagedornParameterSet<D> > parameters, 
                       std::shared_ptr< ShapeEnumeration<D> > enumeration)
        : eps_(eps)
        , parameters_(parameters)
        , enumeration_(enumeration)
        , coefficients_()
    {
        for (dim_t c = 0; c < C; c++)
            coefficients_[c] = std::make_shared< std::valarray<complex_t> >(enumeration->size());
    }
    
    HagedornWavepacket(const HagedornWavepacket<D,C> &other)
        : eps_(other.eps_)
        , parameters_(other.parameters_)
        , enumeration_(other.enumeration_)
        , coefficients_(other.coefficients_)
    { }
    
    HagedornWavepacket<D,C> &operator=(const HagedornWavepacket<D,C> &other)
    {
        eps_ = other.eps_;
        parameters_ = other.parameters_;
        enumeration_ = other.enumeration_;
        coefficients_ = other.coefficients_;
        
        return *this;
    }
    
    HagedornWavepacket<D,1> component(dim_t component) const
    {
        HagedornWavepacket<D,1> slice(parameters_, enumeration_, coefficients_[component]);
        
        return slice;
    }
    
    real_t eps() const
    {
        return eps_;
    }
    
    std::shared_ptr< HagedornParameterSet<D> > parameters() const
    {
        return parameters_;
    }
    
    std::shared_ptr< ShapeEnumeration<D> > enumeration() const
    {
        return enumeration_;
    }
    
    std::shared_ptr< std::valarray<complex_t> > coefficients(dim_t component) const
    {
        return coefficients_[component];
    }
    
    std::array< std::shared_ptr< std::valarray<complex_t> >, C> coefficients() const
    {
        return coefficients_;
    }
    
    Eigen::Matrix<complex_t,C,1> coefficient(std::size_t ordinal) const
    {
        Eigen::Matrix<complex_t,C,1> coeff;
        for (dim_t c = 0; c < C; c++)
            coeff(c,0) = coefficients_[c]->operator[](ordinal);
                
        return coeff;
    }
    
    class Evaluator
    {
    private:
        real_t eps_;
        std::shared_ptr< HagedornParameterSet<D> > parameters_;
        std::shared_ptr< ShapeEnumeration<D> > enumeration_;
        
        Eigen::Matrix<real_t,D,1> x;
        
        std::vector< Eigen::Matrix<complex_t,C,1> > prev_slice_values;
        std::vector< Eigen::Matrix<complex_t,C,1> > curr_slice_values;
        std::vector< Eigen::Matrix<complex_t,C,1> > next_slice_values;
        
        Eigen::Matrix<complex_t,D,D> Qinv;
        Eigen::Matrix<complex_t,D,D> QhQinvt;
        
        int islice;
        
        complex_t evaluateGroundState() const
        {
            const real_t pi = 3.14159265359;
            
            Eigen::Matrix<real_t,D,1> dx = x - parameters_->q;
            
            complex_t pr1 = dx.transpose()*parameters_->P*parameters_->Q.inverse()*dx;
            real_t pr2 = parameters_->p.transpose()*dx;
            
            complex_t exponent = complex_t(0.0, 1.0)/(eps_*eps_) * (0.5*pr1 + pr2);
            
            return 1.0/std::pow(pi*eps_*eps_, D/4.0) * std::exp(exponent);
        }
        
        Eigen::Matrix<complex_t,C,1> evaluateBasis(dim_t axis, 
                                                   MultiIndex<D> k, 
                                                   const Eigen::Matrix<complex_t,C,1> &curr_basis, 
                                                   const Eigen::Matrix<complex_t,C,D> &prev_bases) const
        {
            //compute {sqrt(k[i])*phi[k-e[i]]}
            //  e[i]: unit vector aligned to i-th axis
            Eigen::Matrix<complex_t,D,C> prev_bases_scaled = prev_bases.transpose();
            for (dim_t d = 0; d < D; d++)
                prev_bases_scaled.row(d) *= std::sqrt( real_t(k[d]) );
            
            Eigen::Matrix<real_t,D,1> dx = x - parameters_->q;
            
            complex_t temp = Qinv.row(axis)*dx;
            
            Eigen::Matrix<complex_t,C,1> pr1 = curr_basis.transpose() * std::sqrt(2.0)/eps_ * temp;
            Eigen::Matrix<complex_t,C,1> pr2 = QhQinvt.row(axis)*prev_bases_scaled;
            
            return (pr1 - pr2) / std::sqrt( real_t(k[axis])+1.0);
        }
        
    public:
        Evaluator(const HagedornWavepacket<D,C> &wavepacket, const Eigen::Matrix<real_t,D,1> &x)
            : eps_(wavepacket.eps_)
            , parameters_(wavepacket.parameters_)
            , enumeration_(wavepacket.enumeration_)
            , x(x)
            , prev_slice_values()
            , curr_slice_values()
            , next_slice_values()
            , Qinv(wavepacket.parameters_->Q.inverse())
            , QhQinvt(wavepacket.parameters_->Q.adjoint()*Qinv.transpose())
            , islice(-1)
        { }
        
        bool next_slice()
        {
            if (++islice >= (int)enumeration_->count_slices())
                return false;
            
            if (islice == 0) {
                next_slice_values.push_back( Eigen::Matrix<complex_t,C,1>::Constant( evaluateGroundState() ) );
                return true;
            } else {
                //exchange slices
                std::swap(prev_slice_values, curr_slice_values);
                std::swap(curr_slice_values, next_slice_values);
                next_slice_values = std::vector< Eigen::Matrix<complex_t,C,1> >( enumeration_->slice(islice).size() );
                
                //loop over all multi-indices within next slice [j = position of multi-index within next slice]
                for (std::size_t j = 0; j < enumeration_->slice(islice).size(); j++) {
                    MultiIndex<D> next_index = enumeration_->slice(islice)[j];
                    //find valid precursor: find first non-zero entry
                    dim_t axis = D;
                    for (dim_t d = 0; d < D; d++) {
                        if (next_index[d] != 0) {
                            axis = d;
                            break;
                        }
                    }
                    
                    assert(axis != D); //assert that multi-index contains some non-zero entries
                    
                    //retrieve the basis value within current slice
                    MultiIndex<D> curr_index = next_index; curr_index[axis] -= 1; //get backward neighbour
                    std::size_t curr_ordinal = enumeration_->slice(islice-1).find(curr_index);
                    
                    assert(curr_ordinal < enumeration_->slice(islice-1).size()); //assert that multi-index has been found within current slice
                    Eigen::Matrix<complex_t,C,1> curr_basis = curr_slice_values[curr_ordinal];
                    
                    //retrieve the basis values within previous slice
                    Eigen::Matrix<complex_t,C,D> prev_bases;
                    for (dim_t d = 0; d < D; d++) {
                        if (curr_index[d] == 0) {
                            //precursor is out of shape therefore set this precursor value to zero
                            prev_bases.col(d) = Eigen::Matrix<complex_t,C,1>::Zero();
                        }
                        else {
                            MultiIndex<D> prev_index = curr_index; prev_index[d] -= 1; //get backward neighbour
                            std::size_t prev_ordinal = enumeration_->slice(islice-2).find(prev_index);
                            assert (prev_ordinal < enumeration_->slice(islice-2).size()); //assert that multi-index has been found within previous slice
                            prev_bases.col(d) = prev_slice_values[prev_ordinal];
                        }
                    }
                    
                    //compute basis value within next slice
                    Eigen::Matrix<complex_t,C,1> next_basis = evaluateBasis(axis, curr_index, curr_basis, prev_bases);
                    next_slice_values[j] = next_basis;
                }
                
                return true;
            }
        }
        
        const std::vector< Eigen::Matrix<complex_t,C,1> > &values() const
        {
            return next_slice_values;
        }
        
        std::size_t slice_offset() const
        {
            assert (islice != -1);
            return enumeration_->slice(islice).offset();
        }
    };
    
    /**
     * 
     */
    Eigen::Matrix<complex_t,C,1> operator()(const Eigen::Matrix<real_t,D,1> &x) const
    {
        //Use Kahan's algorithm to accumulate bases with O(1) numerical error instead of O(Sqrt(N))
        KahanSum< Eigen::Matrix<complex_t,C,1> > psi;
        
        Evaluator evaluator(*this, x);
        
        while (evaluator.next_slice()) {
            const auto &values = evaluator.values();
            
            for (std::size_t j = 0; j < values.size(); j++) {
                psi += coefficient(evaluator.slice_offset() + j).cwiseProduct( values[j] );
            }
        }
        
        return psi();
    }
};

}

#endif