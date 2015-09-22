#pragma once
#include "macros.hpp"
#include "matrixPotentials/bases.hpp"
#include "types.hpp"
#include "utilities/evaluations.hpp"

namespace waveblocks
{
  namespace matrixPotentials
  {
    namespace modules
    {
      namespace hessian
      {
        /**
       * \brief Abstract class for potential evaluation
       * 
       * A matrix potential inheriting an implementation of this module
       * can evaluate its potential, jacobian and hessian in one or multiple points
       * 
       * This makes use of the CRTPattern
       * 
       * \tparam Subtype The type extending this interface (used for static polymorphism)
       * \tparam Basis
       * Which basis (bases::Eigen or bases::Canonical) the potential is given in
       * \tparam N
       * Number of levels (dimension of square matrix when evaluated)
       * \tparam D
       * Dimension of argument space
       */
        template <class Subtype, template <int, int> class Basis, int N, int D>
        struct Abstract {
          using Self = Abstract<Subtype, Basis, N, D>;
          IMPORT_TYPES_FROM( Basis, N, D );
          
          
          hessian_evaluation_type evaluate_hessian_at( const CVector<D> &arg ) const {
            return static_cast<const Subtype*>(this)->evaluate_hessian_at_implementation( arg );
          }
          
          template < template <typename...> class grid_in = std::vector,
                   template <typename...> class grid_out = grid_in >
          grid_out<hessian_evaluation_type> evaluate_hessian(
            const grid_in<CVector<D> > &args ) const {
            return utilities::evaluate_function_in_grid < CVector<D>,
                   hessian_evaluation_type,
                   grid_in,
                   grid_out,
                   function_t > (
                     std::bind( &Self::evaluate_hessian_at, this, std::placeholders::_1 ),
                     args );
          }
          
        };
        
        /**
         * \brief Helper class for easier template specialization
         * 
         * This wraps concrete implementations of the Abstract base class
         * 
         * \tparam N
         * Number of levels (dimension of square matrix when evaluated)
         * \tparam D
         * Dimension of argument space
         */
        template <int N, int D>
        struct BasisSpecific {
            struct Canonical
                : Abstract<BasisSpecific<N, D>::Canonical, bases::Canonical, N, D> {
                
                IMPORT_TYPES_FROM( bases::Canonical, N, D );
                
              protected:
                hessian_type hessian;
                
              public:
                Canonical( 
                           hessian_type hessian )
                  :  hessian( hessian ){}
                  
              public:
                
                
                hessian_evaluation_type evaluate_hessian_at_implementation(
                  const CVector<D> &arg ) const {
                  return utilities::evaluate_function_matrix < N,
                         GMatrix,
                         CVector<D>,
                         hessian_return_type,
                         function_t > ( hessian, arg );
                }
            };
            
            struct Eigen : Abstract<BasisSpecific<N, D>::Eigen, bases::Eigen, N, D> {
            
                IMPORT_TYPES_FROM( bases::Eigen, N, D );
                
              protected:
                hessian_type hessian;
                
              public:
                Eigen( 
                       hessian_type hessian )
                  :  hessian( hessian ) {}
                  
              public:
                
                hessian_evaluation_type evaluate_hessian_at_implementation(
                  const CVector<D> &arg ) const {
                  return utilities::evaluate_function_vector < N,
                         GVector,
                         CVector<D>,
                         hessian_return_type,
                         function_t > ( hessian, arg );
                }
            };
        };
        
        template <int D>
        struct BasisSpecific<1, D> {
            template <template <int, int> class Basis>
            struct General : Abstract<BasisSpecific<1, D>::General<Basis>, Basis, 1, D> {
            
                IMPORT_TYPES_FROM( Basis, 1, D );
                
              protected:
                hessian_type hessian;
                
              public:
                General( 
                         hessian_type hessian )
                  : hessian( hessian ) {}
                  
              public:
                
                hessian_evaluation_type evaluate_hessian_at_implementation(
                  const CVector<D> &arg ) const {
                  return hessian( arg );
                }
            };
            
            using Canonical = General<bases::Canonical>;
            using Eigen = General<bases::Eigen>;
        };
        
        // templated typedef with specialization
        template <template <int, int> class Basis, int N, int D>
        struct _HELPER;
        
        template <int N, int D>
        struct _HELPER<bases::Canonical, N, D> {
          using type = typename BasisSpecific<N, D>::Canonical;
        };
        
        template <int N, int D>
        struct _HELPER<bases::Eigen, N, D> {
          using type = typename BasisSpecific<N, D>::Eigen;
        };
        
        template <template <int, int> class Basis, int N, int D>
        using Standard = typename _HELPER<Basis, N, D>::type;
      }
      
      template <template <int, int> class Basis, int N, int D>
      using Hessian = hessian::Standard<Basis, N, D>;
    }
  }
}
