#ifndef WAVEBLOCKS_HAGEDORN_WAVEPACKET
#define WAVEBLOCKS_HAGEDORN_WAVEPACKET

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <array>
#include <memory>
#include <initializer_list>

#include "basic_types.hpp"
#include "math_util.hpp"
#include "hagedorn_parameter_set.hpp"
#include "hagedorn_basis_evaluator.hpp"
#include "shape_enumeration_base.hpp"

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
template<dim_t D>
class HagedornWavepacket
{
public:
    typedef Eigen::Matrix<complex_t,D,1> CMatrixD1;
    typedef Eigen::Matrix<real_t,D,1> RMatrixD1;
    
    /**
     * complex square matrix
     * size: (#dimensions)
     */
    typedef Eigen::Matrix<complex_t,D,D> CMatrixDD;
    
private:
    real_t eps_;
    
    std::shared_ptr< HagedornParameterSet<D> > parameters_;
    std::shared_ptr< ShapeEnumeration<D> > enumeration_;
    
    std::vector<complex_t> coefficients_;
    
public:
    HagedornWavepacket() = default;
    
    /**
     * \brief 
     * 
     * \param[in] paramaters
     * \param[in] enumeration
     * \param coefficients
     * \parblock
     * Vector containing all parameters.<br>
     * Be aware that all content is moved (not copied) into this wavepacket.
     * Therefore the vector is in a indeterminate state afterwards.
     * \endparblock
     */
    HagedornWavepacket(real_t eps,
                       std::shared_ptr< HagedornParameterSet<D> > parameters, 
                       std::shared_ptr< ShapeEnumeration<D> > enumeration,
                       std::vector<complex_t>&& coefficients)
        : eps_(eps)
        , parameters_(parameters)
        , enumeration_(enumeration)
        , coefficients_(std::move(coefficients))
    { }
    
    /**
     * \brief Takes a hagedorn parameter set and shape enumeration. Coefficients are filled with zeros.
     * 
     * \param[in] paramaters 
     * \param[in] enumeration 
     */
    HagedornWavepacket(real_t eps,
                       std::shared_ptr< HagedornParameterSet<D> > parameters, 
                       std::shared_ptr< ShapeEnumeration<D> > enumeration)
        : eps_(eps)
        , parameters_(parameters)
        , enumeration_(enumeration)
        , coefficients_(enumeration.size())
    { }
    
    HagedornWavepacket(const HagedornWavepacket& other) = delete;
    
    /**
     * \brief move constructor
     */
    HagedornWavepacket(HagedornWavepacket&& other)
        : eps_(other.eps_)
        , parameters_(other.parameters_)
        , enumeration_(other.enumeration_)
        , coefficients_(std::move(other.coefficients_))
    { }
    
    HagedornWavepacket &operator=(const HagedornWavepacket& other) = delete;
    
    /**
     * \brief move assignment operator
     * 
     * <blockquote>
     * Move assignement operators typically "steal" the resources held by the argument, 
     * rather than make copies of them, and leave the argument in some valid 
     * but otherwise indeterminate state.
     * </blockquote>
     * 
     */
    HagedornWavepacket &operator=(HagedornWavepacket&& other)
    {
        eps_ = other.eps_;
        parameters_ = other.parameters_;
        enumeration_ = other.enumeration_;
        coefficients_ = std::move(other.coefficients_);
        
        return *this;
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
    
    const std::vector<complex_t>& coefficients() const
    {
        return coefficients_;
    }
    
    std::vector<complex_t>& coefficients()
    {
        return coefficients_;
    }
    
private:
    /**
     * \brief evaluates wavepacket in a memory efficient manner
     * 
     * \tparam N Assigns number of quadrature points. Assign Eigen::Dynamic to do that at runtime.
     */
    template<int N>
    Eigen::Array<complex_t,1,N> evaluate_and_reduce_(const Eigen::Matrix<complex_t,D,N> &x) const
    {
        // determine number of quadrature points since template parameter N could be Eigen::Dynamic
        int npts = x.cols();
        
        // use Kahan's algorithm to accumulate bases with O(1) numerical error instead of O(Sqrt(N))
        KahanSum< Eigen::Array<complex_t,1,N> > psi( Eigen::Matrix<complex_t,1,N>::Zero(1,npts) );
        
        Evaluator<D,N> evaluator(eps_, parameters_, enumeration_, x);
        
        Eigen::Array<complex_t,Eigen::Dynamic,N> prev_basis(0,npts);
        Eigen::Array<complex_t,Eigen::Dynamic,N> curr_basis(0,npts);
        Eigen::Array<complex_t,Eigen::Dynamic,N> next_basis(1,npts);
        
        next_basis = evaluator.seed();
        
        psi += coefficients_[0]*next_basis.row(0).matrix();
        
        for (std::size_t islice = 0; islice < enumeration_->n_slices(); islice++) {
            prev_basis = std::move(curr_basis);
            curr_basis = std::move(next_basis);
            
            next_basis = evaluator.step(islice, prev_basis, curr_basis);
            
            for (long j = 0; j < next_basis.rows(); j++) {
                complex_t cj = coefficients_[enumeration_->slice(islice+1).offset() + j];
                
                //prints: multi-index -> basis -> coefficient
                //std::cout << enumeration_->slice(islice+1)[j] << " -> " << next_basis.row(j).matrix() << " * " << cj << std::endl;
                
                psi += cj*next_basis.row(j).matrix();
            }
        }
        
        return psi();
    }
    
public:
    /**
     * 
     */
    complex_t prefactor() const
    {
        return real_t(1)/parameters_->sqrt_detQ();
    }
    
    /**
     * 
     */
    template<int N>
    Eigen::Matrix<complex_t,1,N> operator()(const Eigen::Matrix<real_t,D,N> &x) const
    {
        Eigen::Matrix<complex_t,D,N> xtmp = x.template cast<complex_t>();
        return evaluate_and_reduce_<N>(xtmp);
    }
    
    /**
     * 
     */
    template<int N>
    Eigen::Matrix<complex_t,1,N> operator()(const Eigen::Matrix<complex_t,D,N> &x) const
    {
        return evaluate_and_reduce_<N>(x);
    }
};

}

#endif