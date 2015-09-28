#pragma once
#include <vector>
#include <functional>
#include "types.hpp"

namespace waveblocks
{
  namespace utilities
  {
    template < int N, int C,
             template <typename, int, int> class M,
             class A,
             class R,
             template <typename...> class F = std::function >
   struct FunctionMatrixEvaluator {
    static M<R, N, C> apply( M<F<R( A )>, N, C> mf, A arg )
    {
      M<R, N, C> m;
      
      for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < C; ++j ) {
          m( i, j ) = mf( i, j )( arg );
        }
      }
      
      return m;
    }
    template<template <typename...> class G_in = std::vector,
             template <typename...> class G_out = G_in>
    static G_out<M<R, N, C> >
     in_grid( M<F<R( A )>, N, C> mf,
        G_in<A> g )
    {
      G_out<M<R, N, C> > result( g.size() );
      auto it = result.begin();
      
      for ( const auto & arg : g ) {
        *it = apply( mf, arg );
        ++it;
      }
      
      return result;
    }
    
   };

  template < template <typename, int, int> class M,
       class A,
       class R,
       template <typename...> class F>
  struct FunctionMatrixEvaluator<1,1,M,A,R,F> {
      static R apply(F<R(A)> f, A arg) {
        return f(arg);
      }
      static  M<R, 1,1> apply( M<F<R( A )>, 1,1 > mf, A arg )
    {
      M<R, 1, 1> m;
          m( 0,0 ) = mf( 0,0 )( arg );
      
      return m;
    }

     template<   template <typename...> class G_in = std::vector,
             template <typename...> class G_out = G_in>
    static G_out<M<R, 1,1> >
     in_grid( M<F<R( A )>, 1,1> mf,
        G_in<A> g )
    {
      G_out<M<R, 1,1> > result( g.size() );
      auto it = result.begin();
      
      for ( const auto & arg : g ) {
        *it = apply( mf, arg );
        ++it;
      }
      
      return result;
    }

      template<  template <typename...> class G_in = std::vector,
             template <typename...> class G_out = G_in>
    static G_out<R >
     in_grid( F<R( A )> mf,
        G_in<A> g )
    {
      G_out<R > result( g.size() );
      auto it = result.begin();
      
      for ( const auto & arg : g ) {
        *it = apply( mf, arg );
        ++it;
      }
      
      return result;
    }
  };

    
    
    /**
     * \brief Evaluate a function in multiple points at once.
     * 
     * Generic template to evaluate a function object in multiple points
     * 
     * \param f 
     * Function object taking input of type A and returning type R
     * \param g
     * grid (template G_in) of A elements 
     * \return
     * grid (template G_out) of R elements
     * 
     * \tparam A
     * Type of the argument to be passed to f
     * \tparam R
     * Return type of f
     * \tparam G_in
     * class template for a container which implements size() and allows for:in loops
     * \tparam G_out
     * class template for a container which implements a constructor with
     *  a single size_t argument and implements begin() to return a forward-iterator 
     * \tparam F
     * class template for a function object
     */
    template < class A,
             class R,
             template <typename...> class G_in = std::vector,
             template <typename...> class G_out = G_in,
             template <typename...> class F = std::function >
    G_out<R> evaluate_function_in_grid( F<R( A )> f, G_in<A> g )
    {
      G_out<R> result( g.size() );
      auto it = result.begin();
      
      for ( const auto & arg : g ) {
        *it = f( arg );
        ++it;
      }
      
      return result;
    }

    /**
     * \brief Evaluate a function square matrix.
     * 
     * Generic template to evaluate allfunction objects in a vector in a single point
     * 
     * \param mf 
     * Matrix (type M<N,N>( of function objects taking input of type A and returning type R
     * \param arg
     * The argument of type A to be passed into all function objects 
     * \return
     * grid (type G_out) of vectors (type V) of N elements of type R which are the results of the function evaluations
     * 
     * \tparam N
     * Dimension of the square matrix to be evaluated
     * \tparam M
     * Class template for a (fixed size) matrix
     * \tparam A
     * Type of the argument to be passed to mf[i]
     * \tparam R
     * Return type of mf[i]
     * \tparam F
     * class template for a function object
     */
    
    /**
     * \brief Evaluate a function square matrix in a grid.
     * 
     * Generic template to evaluate allfunction objects in a vector in multiple points
     * 
     * \param mf 
     * Matrix (type M<N,N>( of function objects taking input of type A and returning type R
     * \param g
     * grid (type G_in) of the arguments of type A to be passed into all function objects 
     * \return
     * grid (type G_out) of vectors (type V) of N elements of type R which are the results of the function evaluations
     * 
     * \tparam N
     * Dimension of the square matrix to be evaluated
     * \tparam M
     * Class template for a (fixed size) matrix
     * \tparam A
     * Type of the argument to be passed to mf[i]
     * \tparam R
     * Return type of mf[i]
     * \tparam G_in
     * class template for a container which implements size() and allows for:in loops
     * \tparam G_out
     * class template for a container which implements a constructor with
     *  a single size_t argument and implements begin() to return a forward-iterator 
     * \tparam F
     * class template for a function object
     */

  }
}
