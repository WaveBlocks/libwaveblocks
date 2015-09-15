//~ #pragma once
//~ #include "macros.hpp"
//~ 
//~ template < template <template <int, int> class, int, int> class S,
         //~ template <int, int> class B,
         //~ int N,
         //~ int D >
//~ struct ExponentialImplemenation
//~ {
  //~ GET_TYPES( ExponentialImplemenation, S, B, N, D );
  //~ 
  //~ potential_evaluation_type evaluate_exponential_at( const RVector<D> &arg,
      //~ real_t factor = 1 ) {
    //~ // Compute matrix
    //~ auto values = that.evaluate( arg );
    //~ potential_evaluation_type result;
    //~ 
    //~ // Compute exponential
    //~ Eigen::MatrixExponential<potential_evaluation_type> m_exp( factor * values );
    //~ 
    //~ m_exp.compute( result );
    //~ return result;
  //~ }
  //~ 
  //~ template < template <typename...> class grid_in = std::vector,
           //~ template <typename...> class grid_out = grid_in >
  //~ grid_out<potential_evaluation_type> evaluate_exponential(
    //~ grid_in<RVector<D> > args,
    //~ real_t factor = 1 ) {
    //~ return evaluate_function_in_grid < RVector<D>,
           //~ potential_evaluation_type,
           //~ grid_in,
           //~ grid_out,
           //~ function_t > (
             //~ std::bind( &Self::evaluate_exponential_at,
                        //~ this,
                        //~ std::placeholders::_1,
                        //~ factor ),
             //~ args );
  //~ }
//~ };
//~ 
//~ template<template <template <int, int> class, int, int> class S, template<int, int> class B1, template<int, int> class B2, int N, int D>
//~ struct TransformableAbstract {
//~ 
    //~ using Self = TransformableAbstract<S, B1, B2, N, D>;
    //~ using Subtype = S<B1, N, D>;
    //~ using Basis = B1<N, D>;
    //~ using potential_evaluation_type = typename Basis::potential_evaluation_type;
    //~ using transformation_type = rD_to_cNxN<D, N>;
    //~ using transformation_return_type = CMatrix<N, N>;
    //~ 
  //~ protected:
    //~ friend Subtype;
    //~ 
    //~ 
    //~ typename B2<N, D>::potential_evaluation_type transformed_evaluate_at( const RVector<D> &arg ) {
      //~ return that.transformed_evaluate_at_implementation( arg );
    //~ }
    //~ 
    //~ template < template <typename...> class grid_in = std::vector,
             //~ template <typename...> class grid_out = grid_in >
    //~ grid_out<potential_evaluation_type> transformed_evaluate(
      //~ grid_in<RVector<D> > args ) {
      //~ return evaluate_function_in_grid < RVector<D>,
             //~ potential_evaluation_type,
             //~ grid_in,
             //~ grid_out,
             //~ function_t > (
               //~ std::bind( &Self::transformed_evaluate_at,
                          //~ this,
                          //~ std::placeholders::_1 ),
               //~ args );
    //~ }
    //~ 
    //~ transformation_return_type evaluate_transformation_at( const RVector<D> &arg ) {
      //~ return that.transfomation( arg );
    //~ }
    //~ 
    //~ template < template <typename...> class grid_in = std::vector,
             //~ template <typename...> class grid_out = grid_in >
    //~ grid_out<potential_evaluation_type> evaluate_transformation(
      //~ grid_in<RVector<D> > args ) {
      //~ return evaluate_function_in_grid < RVector<D>,
             //~ potential_evaluation_type,
             //~ grid_in,
             //~ grid_out,
             //~ function_t > (
               //~ std::bind( &Self::evaluate_transformation_at,
                          //~ this,
                          //~ std::placeholders::_1 ),
               //~ args );
    //~ }
//~ };
//~ 
//~ template<template <template <int, int> class, int, int> class S, template<int, int> class B1, template<int, int> class B2, int N, int D>
//~ struct TransformableImplementation;
//~ 
//~ template<template <template <int, int> class, int, int> class S, int N, int D>
//~ struct TransformableImplementation<S, CanonicalBasis, EigenBasis, N, D> : public TransformableAbstract<S, CanonicalBasis, EigenBasis, N, D> {
    //~ using Subtype = S<CanonicalBasis, N, D>;
    //~ using Basis = CanonicalBasis<N, D>;
    //~ using potential_evaluation_type = typename Basis::potential_evaluation_type;
    //~ using potential_return_type = typename Basis::potential_return_type;
    //~ 
    //~ using Supertype = TransformableAbstract<S, CanonicalBasis, EigenBasis, N, D>;
    //~ using transformation_type = typename Supertype::transformation_type;
    //~ using transformation_return_type = typename Supertype::transformation_return_type;
    //~ 
  //~ protected:
    //~ friend Subtype;
    //~ transformation_type transformation;
    //~ 
  //~ public:
    //~ TransformableImplementation( transformation_type transformation ) :
      //~ transformation( transformation ) {}
      //~ 
    //~ typename EigenBasis<N, D>::potential_evaluation_type transformed_evaluate_at_implementation( const RVector<D> &arg ) {
      //~ auto transform = transformation( arg );
      //~ potential_return_type pot = that.evaluate_at( arg );
      //~ return ( transform.inverse() * pot * transform );
    //~ }
//~ };
//~ 
//~ template<template <template <int, int> class, int, int> class S, int N, int D>
//~ struct TransformableImplementation<S, EigenBasis, CanonicalBasis, N, D> : public TransformableAbstract<S, CanonicalBasis, EigenBasis, N, D> {
    //~ using Subtype = S<EigenBasis, N, D>;
    //~ using Basis = EigenBasis<N, D>;
    //~ using potential_evaluation_type = typename Basis::potential_evaluation_type;
    //~ using potential_return_type = typename Basis::potential_return_type;
    //~ using Supertype = TransformableAbstract<S, EigenBasis, CanonicalBasis, N, D>;
    //~ using transformation_type = typename Supertype::transformation_type;
    //~ using transformation_return_type = typename Supertype::transformation_return_type;
    //~ 
  //~ protected:
    //~ friend Subtype;
    //~ transformation_type transformation;
    //~ 
  //~ public:
    //~ TransformableImplementation( transformation_type transformation ) :
      //~ transformation( transformation ) {}
      //~ 
    //~ typename CanonicalBasis<N, D>::potential_evaluation_type transformed_evaluate_at_implementation( const RVector<D> &arg ) {
      //~ auto transform = transformation( arg );
      //~ transformation_return_type eigs;
      //~ 
      //~ for ( int i = 0; i < N; ++i ) {
        //~ eigs( i, i ) = that.potential[i]( arg );
      //~ }
      //~ 
      //~ return ( transform * eigs * transform.inverse() ).real();
    //~ }
//~ };
//~ 
