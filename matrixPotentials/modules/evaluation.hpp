#pragma once
#include "../macros.hpp"
#include "../bases.hpp"
#include "../types.hpp"
#include "../utilities/evaluations.hpp"

namespace matrixPotentials {
  namespace modules {
    namespace evaluation {
      template <class Subtype,
               template<int,int> class Basis,
               int N,
               int D >
      struct Abstract
      {
        using Self = Abstract<Subtype,Basis,N,D>;
        IMPORT_TYPES_FROM(Basis,N,D);
        
        potential_evaluation_type evaluate_at( const RVector<D> &arg ) {
          return that.evaluate_at_implementation( arg );
        }
        
        template < template <typename...> class grid_in = std::vector,
                 template <typename...> class grid_out = grid_in >
        grid_out<potential_evaluation_type> evaluate(
          const grid_in<RVector<D> > &args ) {
          return utilities::evaluate_function_in_grid < RVector<D>,
                 potential_evaluation_type,
                 grid_in,
                 grid_out,
                 function_t > (
                   std::bind( &Self::evaluate_at, this, std::placeholders::_1 ), args );
        }
        
        jacobian_evaluation_type evaluate_jacobian_at( const RVector<D> &arg ) {
          return that.evaluate_jacobian_at_implementation( arg );
        }
        
        template < template <typename...> class grid_in = std::vector,
                 template <typename...> class grid_out = grid_in >
        grid_out<jacobian_evaluation_type> evaluate_jacobian(
          const grid_in<RVector<D> > &args ) {
          return utilities::evaluate_function_in_grid < RVector<D>,
                 jacobian_evaluation_type,
                 grid_in,
                 grid_out,
                 function_t > (
                   std::bind( &Self::evaluate_jacobian_at, this, std::placeholders::_1 ),
                   args );
        }
        
        hessian_evaluation_type evaluate_hessian_at( const RVector<D> &arg ) {
          return that.evaluate_hessian_at_implementation( arg );
        }
        
        template < template <typename...> class grid_in = std::vector,
                 template <typename...> class grid_out = grid_in >
        grid_out<hessian_evaluation_type> evaluate_hessian(
          const grid_in<RVector<D> > &args ) {
          return utilities::evaluate_function_in_grid < RVector<D>,
                 hessian_evaluation_type,
                 grid_in,
                 grid_out,
                 function_t > (
                   std::bind( &Self::evaluate_hessian_at, this, std::placeholders::_1 ),
                   args );
        }
        
        template <template <typename...> class Tuple = std::tuple>
        Tuple<> taylor_at(const RVector<D>& g ) {
          return Tuple<>(
                   evaluate_at( g ), evaluate_jacobian_at( g ), evaluate_hessian_at( g ) );
        }
        
        template < template <typename...> class Tuple = std::tuple,
                 template <typename...> class grid_in = std::vector,
                 template <typename...> class grid_out = grid_in >
        grid_out<Tuple<> > taylor( const grid_in<RVector<D> > &args ) {
          return utilities::evaluate_function_in_grid < RVector<D>,
                 Tuple<>,
                 grid_in,
                 grid_out,
                 function_t > (
                   std::bind( &Self::taylor_at, this, std::placeholders::_1 ), args );
        }
      };


      template< int N, int D>
      struct BasisSpecific {
        struct Canonical : Abstract < BasisSpecific<N,D>::Canonical, bases::Canonical ,N,D> {
        
          IMPORT_TYPES_FROM(bases::Canonical,N,D);
          
          
        protected:
          potential_type potential;
          jacobian_type jacobian;
          hessian_type hessian;
          
        public:
          Canonical( potential_type potential,
                                    jacobian_type jacobian,
                                    hessian_type hessian )
            : potential( potential ), jacobian( jacobian ), hessian( hessian ) {}
            
        public:
          potential_evaluation_type evaluate_at_implementation( const RVector<D> &arg ) {
            return utilities::evaluate_function_matrix < N,
                   GMatrix,
                   RVector<D>,
                   potential_return_type,
                   function_t > ( potential, arg );
          }
          
          jacobian_evaluation_type evaluate_jacobian_at_implementation(
            const RVector<D> &arg ) {
            return utilities::evaluate_function_matrix < N,
                   GMatrix,
                   RVector<D>,
                   jacobian_return_type,
                   function_t > ( jacobian, arg );
          }
          
          hessian_evaluation_type evaluate_hessian_at_implementation(
            const RVector<D> &arg ) {
            return utilities::evaluate_function_matrix < N,
                   GMatrix,
                   RVector<D>,
                   hessian_return_type,
                   function_t > ( hessian, arg );
          }
          
        };
        
        struct Eigen : Abstract< BasisSpecific<N,D>::Eigen, bases::Eigen, N, D> {
            
          IMPORT_TYPES_FROM(bases::Eigen,N,D);
          
          
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
          potential_evaluation_type evaluate_at_implementation( const RVector<D> &arg ) {
            return utilities::evaluate_function_vector < N,
                   GVector,
                   RVector<D>,
                   potential_return_type,
                   function_t > ( potential, arg );
          }
          
          jacobian_evaluation_type evaluate_jacobian_at_implementation(
            const RVector<D> &arg ) {
            return utilities::evaluate_function_vector < N,
                   GVector,
                   RVector<D>,
                   jacobian_return_type,
                   function_t > ( jacobian, arg );
          }
          
          hessian_evaluation_type evaluate_hessian_at_implementation(
            const RVector<D> &arg ) {
            return utilities::evaluate_function_vector < N,
                   GVector,
                   RVector<D>,
                   hessian_return_type,
                   function_t > ( hessian, arg );
          }
      };


    };
        
        
      template<int D>
      struct BasisSpecific<1,D> {
        
        template<template<int,int> class Basis>
        struct General : Abstract< BasisSpecific<1,D>::General<Basis>, Basis,1,D>{
        
          IMPORT_TYPES_FROM(Basis,1,D);
          
          
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
          potential_evaluation_type evaluate_at_implementation( const RVector<D> &arg ) {
            return potential(arg);
          }
          
          jacobian_evaluation_type evaluate_jacobian_at_implementation(
            const RVector<D> &arg ) {
            return jacobian(arg);
          }
          
          hessian_evaluation_type evaluate_hessian_at_implementation(
            const RVector<D> &arg ) {
            return hessian(arg);
          }
          
        };
        
        using Canonical = General<bases::Canonical>;
        using Eigen = General<bases::Eigen>;
      };
      
      

      // templated typedef with specialization
      template<template<int,int>class Basis,
        int N,
        int D>
      struct _HELPER;

      template<int N, int D>
      struct _HELPER<bases::Canonical,N,D> {
        using type = typename BasisSpecific<N,D>::Canonical;
      };

      template<int N, int D>
      struct _HELPER<bases::Eigen,N,D> {
        using type = typename BasisSpecific<N,D>::Eigen;
      };

      template<template<int,int>class Basis,
        int N,
        int D>
      using Standard = typename _HELPER<Basis,N,D>::type;

    }
    
    template<template<int,int>class Basis,
        int N,
        int D>
      using Evaluation = evaluation::Standard<Basis,N,D>;
  }
}

