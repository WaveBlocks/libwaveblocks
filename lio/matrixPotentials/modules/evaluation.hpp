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
      namespace evaluation
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
          
          potential_evaluation_type evaluate_at( const CVector<D> &arg ) const {
            return static_cast<const Subtype*>(this)->evaluate_at_implementation( arg );
          }
          
          template < template <typename...> class grid_in = std::vector,
                   template <typename...> class grid_out = grid_in >
          grid_out<potential_evaluation_type> evaluate(
            const grid_in<CVector<D> > &args ) const {
            return utilities::evaluate_function_in_grid < CVector<D>,
                   potential_evaluation_type,
                   grid_in,
                   grid_out,
                   function_t > (
                     std::bind( &Self::evaluate_at, this, std::placeholders::_1 ), args );
          }
          
          jacobian_evaluation_type evaluate_jacobian_at( const CVector<D> &arg ) const {
            return static_cast<const Subtype*>(this)->evaluate_jacobian_at_implementation( arg );
          }
          
          template < template <typename...> class grid_in = std::vector,
                   template <typename...> class grid_out = grid_in >
          grid_out<jacobian_evaluation_type> evaluate_jacobian(
            const grid_in<CVector<D> > &args ) const {
            return utilities::evaluate_function_in_grid < CVector<D>,
                   jacobian_evaluation_type,
                   grid_in,
                   grid_out,
                   function_t > (
                     std::bind( &Self::evaluate_jacobian_at, this, std::placeholders::_1 ),
                     args );
          }
          
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
          
          template <template <typename...> class Tuple = std::tuple>
          Tuple<potential_evaluation_type, jacobian_evaluation_type, hessian_evaluation_type> taylor_at( const CVector<D> &g ) const {
            return Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>(
                     evaluate_at( g ), evaluate_jacobian_at( g ), evaluate_hessian_at( g ) );
          }
          
          template < template <typename...> class Tuple = std::tuple,
                   template <typename...> class grid_in = std::vector,
                   template <typename...> class grid_out = grid_in >
          grid_out< Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>>taylor( const grid_in<CVector<D> > &args ) const {
            return utilities::evaluate_function_in_grid < CVector<D>,
                   Tuple<>,
                   grid_in,
                   grid_out,
                   function_t > (
                     std::bind( &Self::taylor_at, this, std::placeholders::_1 ), args );
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
                potential_type potential;
                jacobian_type jacobian;
                hessian_type hessian;
                
              public:
                Canonical( potential_type potential,
                           jacobian_type jacobian,
                           hessian_type hessian )
                  : potential( potential ), jacobian( jacobian ), hessian( hessian ){}
                  
              public:
                potential_evaluation_type evaluate_at_implementation(
                  const CVector<D> &arg ) const {
                  return utilities::evaluate_function_matrix < N,
                         GMatrix,
                         CVector<D>,
                         potential_return_type,
                         function_t > ( potential, arg );
                }
                
                jacobian_evaluation_type evaluate_jacobian_at_implementation(
                  const CVector<D> &arg ) const {
                  return utilities::evaluate_function_matrix < N,
                         GMatrix,
                         CVector<D>,
                         jacobian_return_type,
                         function_t > ( jacobian, arg );
                }
                
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
                potential_type potential;
                jacobian_type jacobian;
                hessian_type hessian;
                
              public:
                Eigen( potential_type potential,
                       jacobian_type jacobian,
                       hessian_type hessian )
                  : potential( potential ), jacobian( jacobian ), hessian( hessian ) {}
                  
              public:
                potential_evaluation_type evaluate_at_implementation(
                  const CVector<D> &arg ) const {
                  return utilities::evaluate_function_vector < N,
                         GVector,
                         CVector<D>,
                         potential_return_type,
                         function_t > ( potential, arg );
                }
                
                jacobian_evaluation_type evaluate_jacobian_at_implementation(
                  const CVector<D> &arg ) const {
                  return utilities::evaluate_function_vector < N,
                         GVector,
                         CVector<D>,
                         jacobian_return_type,
                         function_t > ( jacobian, arg );
                }
                
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
                potential_type potential;
                jacobian_type jacobian;
                hessian_type hessian;
                
              public:
                General( potential_type potential,
                         jacobian_type jacobian,
                         hessian_type hessian )
                  : potential( potential ), jacobian( jacobian ), hessian( hessian ) {}
                  
              public:
                potential_evaluation_type evaluate_at_implementation(
                  const CVector<D> &arg ) const {
                  return potential( arg );
                }
                
                jacobian_evaluation_type evaluate_jacobian_at_implementation(
                  const CVector<D> &arg ) const {
                  return jacobian( arg );
                }
                
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
      using Evaluation = evaluation::Standard<Basis, N, D>;
    }
  }
}
