#pragma once
#include"macros.hpp"

template<template<template<int,int> class ,int,int> class S, template<int,int> class B, int N, int D>
class ExponentialImplemenation {
	GET_TYPES(ExponentialImplemenation,S,B,N,D);
	
	potential_evaluation_type evaluate_exponential_at(RVector<D> arg, real_t factor = 1) {
		// Compute matrix
		auto values = that.evaluate(arg);
		potential_evaluation_type result;
		
		// Compute exponential
		Eigen::MatrixExponential<potential_evaluation_type> m_exp(factor * values);

		m_exp.compute(result);
		return result;
	}
	
	template<template<typename ...> class grid_in = std::vector, template<typename ...> class grid_out = grid_in>
	grid_out<potential_evaluation_type> evaluate_exponential(grid_in<RVector<D>> args, real_t factor = 1) {
		return evaluate_function_in_grid<RVector<D>, potential_evaluation_type, grid_in,grid_out,function_t>(std::bind(&self_type::evaluate_exponential_at,this,std::placeholders::_1,factor),args);
	}
};




//~ template<template<int, int> class B1, template<int,int> class B2, int N, int D>
//~ class TransformableMatrixPotential : public AbstractMatrixPotential<B1,N,D> {
	//~ private:
	//~ typename B1::transformation_type transformation;
	//~ 
	//~ public:
	
	//~ B2::potential_evaluation_type other_evaluate_at(RVector<D> arg) {
		//~ auto transform = transform(arg);
		//~ CMatrix<N,N> eigs;
		//~ for (int d = 0; d < D; ++d) {
			//~ eigs(d,d) = eigen_potential[d](arg);
		//~ }
		//~ return (transform*D*transform.inverse()).real();
	//~ }
	//~ 
	//~ Grid<RMatrix<N,N>> values = canonical_evaluate_at(g);
				//~ Grid<CVector<N>> res(values.size());
				//~ auto it = res.begin();
				//~ for (auto& matrix: values) {
					//~ Eigen::EigenSolver<RMatrix<N,N>> es(matrix,false);
					//~ *it = es.eigenvalues();
					//~ ++it;
				//~ }
				//~ return res;
				//~ 
				 //~ Grid<CMatrix<N,N>> evaluate_eigenvectors_at(Grid<RVector<D>> g) {
			//~ if (CanonicalBasisImplementation) {
				//~ Grid<RMatrix<N,N>> values = canonical_evaluate_at(g);
				//~ Grid<CMatrix<N,N>> res(values.size());
				//~ auto it = res.begin();
				//~ for (auto& matrix: values) {
					//~ Eigen::EigenSolver<RMatrix<N,N>> es(matrix);
					//~ *it = es.eigenvectors();
					//~ ++it;
				//~ }
				//~ return res;
			//~ }
			//~ else {
				//~ return evaluate_function_in_grid<RVector<D>, CMatrix<N,N>, Grid, function_t>(eigen_basis_transform,g);
			//~ }				
        //~ }
//~ };

