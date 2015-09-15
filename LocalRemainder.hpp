#pragma once
#include "macros.hpp"
#include "Bases.hpp"
#include "LocalQuadratic.hpp"

template<class Subtype,
         int N,
         int D >
class LocalRemainderAbstract {
  typename EigenBasis<N,D>::potential_evaluation_type evaluate_local_remainder_at(const RVector<D>& arg) {
    return that.evaluate_local_remainder_at_implementation(arg);
  }
};



template <template<template<int,int> class ,int,int> class EvalImpl, int N, int D>
class LocalRemainderHomogenousImplementation : public LocalRemainderAbstract<LocalRemainderHomogenousImplementation<EvalImpl, N,D>,N,D>, EvalImpl<CanonicalBasis,N,D> {
  IMPORT_TYPES_FROM(CanonicalBasis,N,D);

  using leading_level_type = EigenBasis<1,D>;
  
  LocalQuadraticImplementation<EvalImpl<EigenBasis,1,D>,EigenBasis,1,D> quadratic;
public:
  LocalRemainderHomogenousImplementation(
    potential_type pot,
    jacobian_type jac,
    hessian_type hess,
    typename leading_level_type::potential_type lead_pot,
    typename leading_level_type::jacobian_type lead_jac,
    typename leading_level_type::hessian_type lead_hess
  ) : EvalImpl<CanonicalBasis,N,D>(pot,jac,hess), quadratic(lead_pot,lead_jac,lead_hess) {} 
  
    typename EigenBasis<N,D>::potential_evaluation_type evaluate_local_remainder_at_implementation(const RVector<D>& arg) {
      auto u = quadratic.evaluate_local_quadratic_at(arg);
      auto V = EvalImpl<CanonicalBasis,N,D>::potential(arg);
      for (int i = 0; i < N; ++i) {
        V(i,i) = V(i,i) - u;
      }
      return V;
    }  
};

template <template<template<int,int> class,int,int> class EvalImpl, int N, int D>
class LocalRemainderInhomogenousImplementation : public LocalRemainderAbstract<LocalRemainderInhomogenousImplementation<EvalImpl,N,D>,N,D>, EvalImpl<CanonicalBasis,N,D> {
  IMPORT_TYPES_FROM(CanonicalBasis,N,D);
  
  using leading_level_type = EigenBasis<N,D>;

  LocalQuadraticImplementation<EvalImpl<EigenBasis,N,D>,EigenBasis,N,D> quadratic;
public:
  LocalRemainderInhomogenousImplementation(
    typename leading_level_type::potential_type pot,
    typename leading_level_type::jacobian_type jac,
    typename leading_level_type::hessian_type hess
  ) : quadratic(pot,jac,hess) {} 
  
  
    typename EigenBasis<N,D>::potential_evaluation_type evaluate_local_remainder_at_implementation(const RVector<D>& arg) {
      auto u = quadratic.evaluate_local_quadratic_at(arg);
      auto V = EvalImpl<CanonicalBasis,N,D>::potential(arg);
      for (int i = 0; i < N; ++i) {
        V(i,i) = V(i,i) - u(i);
      }
      return V;
    }  
};
