#include<vector>
#include<iostream>
#include<array>
#include<list>
#include<functional>
#include<cmath>
#include<Eigen/Core>
#include<Eigen/Eigenvalues>
#include<unsupported/Eigen/MatrixFunctions>
#include"types.hpp"
#include"utilities.hpp"

template<template<int, int> class B, int N, int D>
class MatrixPotential;

//~ // Definition of basis
template<int N, int D>
struct EigenBasis {
	//~ typedef GVector<rD_to_c<D>,N> potential_type;
	//~ typedef GVector<rD_to_cD<D>,N> jacobian_type;
	//~ typedef GVector<rD_to_cDxD<D>,N> hessian_type;
	//~ typedef rD_to_cNxN<D,N> transformation_type;
	//~ typedef rD_to_function_vector<D,N,rD_to_c<D>> local_quadratic_type;
    //~ typedef rD_to_function_vector<D,N,rD_to_c<D>> local_remainder_type;
    //~ 
    //~ typedef CVector<N> potential_evaluation_type;
    //~ typedef GVector<CVector<D>,N> jacobian_evaluation_type;
	//~ typedef GVector<CMatrix<D,D>,N> hessian_evaluation_type;
//~ 
        	//~ 
	//~ static potential_evaluation_type evaluate_at(const MatrixPotential< ::EigenBasis,N,D>& that, RVector<D> arg) {
		//~ return evaluate_function_vector<N,GVector,RVector<D>,complex_t,function_t>(that.potential,arg);
	//~ }
	//~ 
	//~ static jacobian_evaluation_type evaluate_jacobian_at(const MatrixPotential< ::EigenBasis,N,D>& that, RVector<D> arg) {
		//~ return evaluate_function_vector(that.jacobian,arg);
	//~ }
	//~ 
	//~ static hessian_evaluation_type evaluate_hessian_at(const MatrixPotential< ::EigenBasis,N,D>& that, RVector<D> arg) {
		//~ return evaluate_function_matrix(that.hessian, arg);
	//~ }

};

template<int N, int D>
struct CanonicalBasis : MatrixPotential< ::CanonicalBasis, N, D>{
	typedef GMatrix<rD_to_r<D>,N,N> potential_type;
	typedef GMatrix<rD_to_rD<D>,N,N> jacobian_type;
	typedef GMatrix<rD_to_rDxD<D>,N,N> hessian_type;
	typedef rD_to_function_matrix<D,N,rD_to_r<D>> local_quadratic_type;
    typedef rD_to_function_matrix<D,N,rD_to_r<D>> local_remainder_type;
    
	typedef RMatrix<N,N> potential_evaluation_type;
	typedef GMatrix<RVector<D>,N,N> jacobian_evaluation_type;
	typedef GMatrix<RMatrix<D,D>,N,N> hessian_evaluation_type;

	// ::CanonicalBasis is NOT a typo. Clang does not implement the injected-class-name rule yet even though its in the standard.
	// This is a workaround. compare http://stackoverflow.com/a/17687568/3139931
	static potential_evaluation_type evaluate_at(const MatrixPotential< ::CanonicalBasis,N,D>& that, RVector<D> arg) {
		return evaluate_function_matrix<N,GMatrix,RVector<D>,real_t,function_t>(that.potential,arg);
	}
	
	static jacobian_evaluation_type evaluate_jacobian_at(const MatrixPotential< ::CanonicalBasis,N,D>& that, RVector<D> arg) {
		return evaluate_function_matrix(that.jacobian,arg);
	}
	
	static hessian_evaluation_type evaluate_hessian_at(const MatrixPotential< ::CanonicalBasis,N,D>& that, RVector<D> arg) {
		return evaluate_function_matrix(that.hessian, arg);
	}

};
			
//General definition			
template<template<int, int> class B, int N, int D>
class MatrixPotential {
	public:
	using Basis = B<N,D>;
	
	private:
	friend Basis;
	typename Basis::potential_type potential;
	typename Basis::jacobian_type jacobian;
	typename Basis::hessian_type hessian; 
	typename Basis::local_quadratic_type local_quadratic;
	typename Basis::local_remainder_type local_remainder;
	
	public:
	MatrixPotential(typename Basis::potential_type potential,
					typename Basis::jacobian_type jacobian,
					typename Basis::hessian_type hessian
		) : potential(potential),
			jacobian(jacobian),
			hessian(hessian)
		{
			
		}
		
	typename Basis::potential_evaluation_type evaluate_at(RVector<D> arg) {
		return Basis::evaluate_at(*this, arg);	
	}
	
	typename Basis::jacobian_evaluation_type evaluate_jacobian_at(RVector<D> arg) {
		return Basis::evaluate_jacobian_at(*this, arg);
	}
	
	typename Basis::hessian_evaluation_type evaluate_hessian_at(RVector<D> arg) {
		return Basis::evaluate_hessian_at(*this, arg);
	}
};


//~ 
//~ template<template<int, int> class B, int N, int D>
//~ typename B<N,D>::potential_evaluation_type MatrixPotential<B,N,D>::evaluate_at(RVector<D> arg) {
	//~ return evaluate_function_matrix<N,GMatrix,RVector<D>,real_t,function_t>(potential,arg);
//~ }



//~ template<int N, int D>
//~ class MatrixPotential<CanonicalBasis,N,D> {
    //~ private:
		//~ 
        //~ typename CanonicalBasis<N,D>::potential_type potential;
        //~ typename CanonicalBasis<N,D>::jacobian_type jacobian;
        //~ typename CanonicalBasis<N,D>::hessian_type hessian; 
        //~ typename CanonicalBasis<N,D>::transformation_type transformation;
		//~ typename CanonicalBasis<N,D>::local_quadratic_type local_quadratic;
		//~ typename CanonicalBasis<N,D>::local_remainder_type local_remainder;
        //~ 
    //~ public:
		//~ MatrixPotential(typename CanonicalBasis<N,D>::potential_type potential,
						//~ typename CanonicalBasis<N,D>::jacobian_type jacobian,
						//~ typename CanonicalBasis<N,D>::hessian_type hessian,
						//~ typename CanonicalBasis<N,D>::transformation_type transformation
			//~ ) : potential(potential),
				//~ jacobian(jacobian),
				//~ hessian(hessian),
				//~ transformation(transformation)
			//~ {
				//~ calculate_local_quadratic();
				//~ calculate_local_remainder();	
			//~ };
			//~ 
		//~ RMatrix<N,N> evaluate_at(RVector<D>);
		//~ 
		//~ template<template<typename ...> class Grid_in = std::vector, template<typename ...> class Grid_out = Grid_in>
		//~ Grid_out<RMatrix<N,N>> evaluate(Grid_in<RVector<D>> g);
		//~ 
		//~ 
		//~ template<template<typename ...> class Grid_in = std::vector, template<typename ...> class Grid_out = Grid_in>
		//~ Grid_out<CVector<N>> eigen_evaluate(Grid_in<RVector<D>> g);
		//~ 
//~ };
//~ 
//~ template<int N, int D>
//~ RMatrix<N,N> MatrixPotential<CanonicalBasis,N,D>::evaluate_at(RVector<D> arg) {
	//~ return evaluate_function_matrix<N,GMatrix,RVector<D>,real_t,function_t>(potential,arg);
//~ }
//~ 
//~ template<int N, int D>
//~ template<template<typename ...> class Grid_in, template<typename ...> class Grid_out>
//~ Grid_out<RMatrix<N,N>> MatrixPotential<CanonicalBasis, N,D>::evaluate(Grid_in<RVector<D>> g) {
	//~ // std:bind magic 
	//~ return evaluate_function_in_grid<RVector<D>, RMatrix<N,N>, Grid_in, Grid_out, function_t>(std::bind(&MatrixPotential<CanonicalBasis,N,D>::evaluate_at, this, std::placeholders::_1),g);	
//~ }
//~ 
//~ template<int N, int D>
//~ template<template<typename ...> class Grid_in, template<typename ...> class Grid_out>
	//~ Grid_out<CVector<N>> MatrixPotential<CanonicalBasis, N, D>::eigen_evaluate(Grid_in<RVector<D>> g) {
	//~ Grid_out<RMatrix<N,N>> values = evaluate(g);
	//~ Grid_out<CVector<N>> res(values.size());
	//~ auto it = res.begin();
	//~ for (auto& matrix: values) {
		//~ Eigen::EigenSolver<RMatrix<N,N>> es(matrix,false);
		//~ *it = es.eigenvalues();
		//~ ++it;
	//~ }
	//~ return res;
//~ }

//~ template<int N, int D>
//~ class MatrixPotential<EigenBasis,N,D> {
    //~ private:
		//~ 
        //~ typename EigenBasis<N,D>::potential_type potential;
        //~ typename EigenBasis<N,D>::jacobian_type jacobian;
        //~ typename EigenBasis<N,D>::hessian_type hessian; 
        //~ typename EigenBasis<N,D>::transformation_type transformation;
		//~ typename EigenBasis<N,D>::local_quadratic_type local_quadratic;
		//~ typename EigenBasis<N,D>::local_remainder_type local_remainder;
        //~ 
    //~ public:
		//~ MatrixPotential(typename EigenBasis<N,D>::potential_type potential,
						//~ typename EigenBasis<N,D>::jacobian_type jacobian,
						//~ typename EigenBasis<N,D>::hessian_type hessian,
						//~ typename EigenBasis<N,D>::transformation_type transformation
			//~ ) : potential(potential),
				//~ jacobian(jacobian),
				//~ hessian(hessian),
				//~ transformation(transformation)
			//~ {
				//~ calculate_local_quadratic();
				//~ calculate_local_remainder();	
			//~ };
			//~ 
//~ 
		//~ RMatrix<N,N> canonical_evaluate_at(RVector<D> arg);
		//~ 
		//~ template<template<typename ...> class Grid_in = std::vector, template<typename ...> class Grid_out = Grid_in>
		//~ Grid_out<CVector<N>> evaluate(Grid_in<RVector<D>> g);
//~ };


// Evaluation	
	
	//~ template<int N, int D>
	//~ RMatrix<N,N> MatrixPotential<EigenBasis,N,D>::canonical_evaluate_at(RVector<D> arg) {
		//~ auto transform = transformation(arg);
		//~ CMatrix<N,N> eigs;
		//~ for (int d = 0; d < D; ++d) {
			//~ eigs(d,d) = potential[d](arg);
		//~ }
		//~ return (transform*D*transform.inverse()).real();			
	//~ }
//~ 
	//~ 
        //~ 
	//~ template<int N, int D>
	//~ template<template<typename ...> class Grid>
	//~ Grid<CVector<N>> MatrixPotential<EigenBasis,N, D>::evaluate(Grid<RVector<D>> g) {
		//~ return evaluate_function_vector_in_grid<N, GVector, RVector<D>, complex_t, Grid, function_t>(potential,g);
	//~ }
        //~ 
        
        
        
        
               //~ 
	//~ // Jacobian
        //~ // Evaluates the jacobian of the canonical matrix in multiple points
        //~ template<template<typename ...> class Grid = std::vector>
        //~ Grid<GMatrix<RVector<D>,N,N>>canonical_evaluate_jacobian_at(Grid<RVector<D>> g) {
			//~ if (canonicalBasis) {
				//~ return evaluate_function_matrix_in_grid(canonical_jacobian,g);
			//~ }
        //~ }
        //~ 
        //~ // Evaluates the jacobians of the eigen functions in multiple points
        //~ template<template<typename ...> class Grid = std::vector>
        //~ Grid<GVector<RVector<D>,N>> eigen_evaluate_jacobian_at(Grid<RVector<D>> g) {
			//~ if (!canonicalBasis) {
				//~ return evaluate_function_vector_in_grid(eigen_jacobian,g);
			//~ }
        //~ }
//~ 
        //~ 
	//~ // Hessian
        //~ // Evaluates the hessian of the canonical matrix in multiple points
        //~ template<template<typename ...> class Grid = std::vector>
        //~ Grid<GMatrix<RMatrix<D,D>,N,N>> canonical_evaluate_hessian_at(Grid<RVector<D>> g) {
			//~ if (canonicalBasis) {
				//~ return evaluate_function_matrix_in_grid(canonical_hessian, g);
			//~ }
        //~ }
        //~ 
        //~ // Evaluates the hessians of the eigen functions in multiple points
        //~ template<template<typename ...> class Grid = std::vector>
        //~ Grid<GVector<RMatrix<D,D>,N>> eigen_evaluate_hessian_at(Grid<RVector<D>> g) {
			//~ if (!canonicalBasis) {
				//~ return evaluate_function_vector_in_grid(eigen_hessian, g);
			//~ }
        //~ }
        //~ 
	//~ // Eigenvectors
		//~ // Evaluates the eigen basis transformation in a grid of points
        //~ template<template<typename ...> class Grid = std::vector>
        //~ Grid<CMatrix<N,N>> evaluate_eigenvectors_at(Grid<RVector<D>> g) {
			//~ if (canonicalBasis) {
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
//~ 
	//~ // Exponential
		//~ // Computes the exponential of the canonical matrix
        //~ RMatrix<N,N> canonical_evaluate_exponential(RVector<D> arg, real_t factor = 1) {
			//~ // Compute matrix
			//~ RMatrix<N,N> values = canonical_evaluate(arg);
			//~ RMatrix<N,N> result;
			//~ 
			//~ // Compute exponential
			//~ Eigen::MatrixExponential<RMatrix<N,N>> m_exp(factor * values);
//~ 
			//~ m_exp.compute(result);
			//~ return result;
        //~ }
        //~ 
        //~ // Computes the exponential of the canonical matrix in a grid
		//~ template<template<typename ...> class Grid = std::vector>
		//~ Grid<RMatrix<N,N>> canonical_evaluate_exponential_at(Grid<RVector<D>> g, real_t factor = 1) {
			//~ return evaluate_function_in_grid<RVector<D>, RMatrix<N,N>, Grid, function_t>(std::bind(&MatrixPotential<Basis,N,D>::canonical_evaluate_exponential,this,std::placeholders::_1),g);	
		//~ }
        //~ 
	//~ // Taylor
        //~ // Returns the potential, jacobian and hessian
        //~ template<template<typename ...> class Tuple = std::tuple, template<typename ...> class Grid = std::vector>
        //~ Tuple<> canonical_taylor(Grid<RVector<D>> g){
            //~ return Tuple<>(canonical_evaluate_at(g),canonical_evaluate_jacobian_at(g),canonical_evaluate_hessian_at(g));
        //~ }
//~ 
	//~ // Local quadratic
		//~ // Computes a function which returns the quadratic approximation function matrix for a given point q.
        //~ void canonical_calculate_local_quadratic() {
			//~ if (canonicalBasis) {
				//~ canonical_local_quadratic = [=] (RVector<D> q) {
					//~ GMatrix<rD_to_r<D>,N,N> result_matrix;
							//~ 
					//~ for (int l = 0; l < N; ++l) {
						//~ for (int m = 0; m < N; ++m) {
							//~ result_matrix(l,m) = [=](RVector<D> x) {
								//~ real_t V = canonical_potential(l,m)(q);
								//~ RVector<D> J =  canonical_jacobian(l,m)(q);
								//~ RMatrix<D,D> H = canonical_hessian(l,m)(q);
//~ 
								//~ real_t result = V;
								//~ for(int i = 0; i < D; ++i) {
									//~ real_t xmqi = x[i] - q[i];
									//~ result += J[i]*(xmqi);
									//~ for (int j = 0; j < D; ++j) {
										//~ result += 0.5*xmqi*H(i,j) * (x[j]-q[j]);
									//~ }
								//~ }
								//~ return result;
							//~ };
						//~ }
					//~ }
					//~ return result_matrix;
				//~ };
			//~ }
		//~ }
		//~ 
		//~ // Computes the local quadratic approximations of the eigenfunctions
		//~ void eigen_calculate_local_quadratic() {
			//~ if (!canonicalBasis) {
				//~ eigen_local_quadratic = [=] (RVector<D> q) {
					//~ GVector<rD_to_c<D>,N> result_vector;
							//~ 
					//~ for (int l = 0; l < N; ++l) {
						//~ result_vector[l] = [=](RVector<D> x) {
							//~ complex_t V = eigen_potential(l)(q);
							//~ CVector<D> J =  eigen_jacobian(l)(q);
							//~ CMatrix<D,D> H = eigen_hessian(l)(q);
//~ 
							//~ complex_t result = V;
							//~ for(int i = 0; i < D; ++i) {
								//~ real_t xmqi = x[i] - q[i];
								//~ result += J[i]*(xmqi);
								//~ for (int j = 0; j < D; ++j) {
									//~ result += 0.5*xmqi*H(i,j) * (x[j]-q[j]);
								//~ }
							//~ }
							//~ return result;
						//~ };
					//~ }
					//~ return result_vector;
				//~ };
			//~ }
		//~ }
	//~ 
	//~ // Local remainder	
		//~ // Computes a function which returns the local remainder around a point q.
        //~ void canonical_calculate_local_remainder() {
            //~ canonical_local_remainder = [=] (RVector<D> q) {
				//~ GMatrix<rD_to_r<D>,N,N> result;
			    //~ GMatrix<rD_to_r<D>,N,N> local_quadratic = canonical_local_quadratic(q);
			    //~ 
				//~ for (int i = 0; i < N; ++i) {
					//~ for (int j = 0; j < N; ++j) {
						//~ result(i,j) = [=](RVector<D> x) {
							//~ return canonical_potential(i,j)(x) - local_quadratic(i,j)(x);
						//~ };
					//~ }
				//~ }
				//~ return result;
			//~ };
		//~ }
//~ 
        //~ // Evaluates the local remainder matrix in a grid.
        //~ template<template<typename ...> class Grid = std::vector>
		//~ Grid<RMatrix<N,N>> canonical_evaluate_local_remainder_at(Grid<RVector<D>> g, RVector<D> position) {
            //~ return evaluate_function_matrix_in_grid<N,GMatrix,RVector<D>,real_t,Grid,function_t>(canonical_local_remainder(position),g);
        //~ }
//~ 
		//~ 
		//~ // Computes a function which returns the local remainder around a point q.
        //~ void eigen_calculate_local_remainder() {
            //~ eigen_local_remainder = [=] (RVector<D> q) {
				//~ GVector<rD_to_c<D>,N> result;
			    //~ GVector<rD_to_c<D>,N> local_quadratic = eigen_local_quadratic(q);
			    //~ 
				//~ for (int i = 0; i < N; ++i) {
					//~ result[i] = [=](RVector<D> x) {
						//~ return eigen_potential(i)(x) - local_quadratic(i)(x);
					//~ };
				//~ }
				//~ return result;
			//~ };
		//~ }
		//~ 
		//~ // Evaluates the local remainder vector in a grid.
        //~ template<template<typename ...> class Grid = std::vector>
		//~ Grid<RVector<N>> eigen_evaluate_local_remainder_at(Grid<RVector<D>> g, RVector<D> position) {
            //~ return evaluate_function_vector_in_grid<N,GVector,RVector<D>,real_t,Grid,function_t>(eigen_local_remainder(position),g);
        //~ }
				//~ 
//~ 
//~ };


//~ // Scalar specialization
//~ template<class Basis>
//~ class MatrixPotential<Basis> {
    //~ private:
        //~ rD_to_r<D> canonical_potential;
        //~ rD_to_rD<D> canonical_jacobian;
        //~ rD_to_rDxD<D> canonical_hessian;
        //~ 
        //~ rD_to_function<D,rD_to_r<D>> canonical_local_quadratic;
        //~ rD_to_function<D,rD_to_r<D>> canonical_local_remainder;
        //~ 
    //~ public:
        //~ MatrixPotential(rD_to_r<D> potential,
						//~ rD_to_rD<D> jacobian,
						//~ rD_to_rDxD<D>  hessian
			//~ ) : canonical_potential(potential),
				//~ canonical_jacobian(jacobian),
				//~ canonical_hessian(hessian)
			//~ {
				//~ canonical_calculate_local_quadratic();
				//~ canonical_calculate_local_remainder();	
			//~ };
		//~ 
	//~ // Evaluation
		//~ // Evaluates the canonical matrix in one point
        //~ real_t evaluate(RVector<D> arg) {
			//~ return canonical_potential(arg);
        //~ }
        //~ 
		//~ // Evaluates the canonical matrix in multiple points
        //~ template<template<typename ...> class Grid = std::vector>
        //~ Grid<real_t> evaluate_at(Grid<RVector<D>> g) {
			//~ return evaluate_function_in_grid<RVector<D>, real_t, Grid, function_t>(canonical_potential,g);	
        //~ }
        //~ 
    //~ // Jacobian    
        //~ // Evaluates the jacobian of the canonical matrix in multiple points
        //~ template<template<typename ...> class Grid = std::vector>
        //~ Grid<RVector<D>> evaluate_jacobian_at(Grid<RVector<D>> g) {
			//~ return evaluate_function_in_grid(canonical_jacobian,g);
        //~ }
    //~ // Hessian
        //~ // Evaluates the hessian of the canonical matrix in multiple points
        //~ template<template<typename ...> class Grid = std::vector>
        //~ Grid<RMatrix<D,D>> evaluate_hessian_at(Grid<RVector<D>> g) {
			//~ return evaluate_function_in_grid(canonical_hessian, g);
        //~ }
        //~ 
	//~ // Exponential
		//~ // Computes the exponential of the canonical matrix
        //~ real_t evaluate_exponential(RVector<D> arg, real_t factor = 1) {
			//~ return exp(canonical_potential(arg));
        //~ }
        //~ 
        //~ // Computes the exponential of the canonical matrix in a grid
		//~ template<template<typename ...> class Grid = std::vector>
		//~ Grid<real_t> evaluate_exponential_at(Grid<RVector<D>> g, real_t factor = 1) {
			//~ return evaluate_function_in_grid<RVector<D>, real_t, Grid, function_t>(std::bind(&MatrixPotential<1,D>::canonical_evaluate_exponential,this,std::placeholders::_1),g);	
		//~ }
        //~ 
	//~ // Taylor
        //~ // Returns the potential, jacobian and hessian
        //~ template<template<typename ...> class Tuple = std::tuple, template<typename ...> class Grid = std::vector>
        //~ Tuple<> taylor(Grid<RVector<D>> g){
            //~ return Tuple<>(evaluate_at(g),evaluate_jacobian_at(g),evaluate_hessian_at(g));
        //~ }
        //~ 
	//~ // Local quadratic
		//~ // Computes a function which returns the quadratic approximation function matrix for a given point q.
        //~ void canonical_calculate_local_quadratic() {
			//~ canonical_local_quadratic = [=] (RVector<D> q) {
				//~ return [=](RVector<D> x) {
					//~ real_t V = canonical_potential(q);
					//~ RVector<D> J =  canonical_jacobian(q);
					//~ RMatrix<D,D> H = canonical_hessian(q);
//~ 
					//~ real_t result = V;
					//~ for(int i = 0; i < D; ++i) {
						//~ real_t xmqi = x[i] - q[i];
						//~ result += J[i]*(xmqi);
						//~ for (int j = 0; j < D; ++j) {
							//~ result += 0.5*xmqi*H(i,j) * (x[j]-q[j]);
						//~ }
					//~ }
					//~ return result;
				//~ };
			//~ };
		//~ }
		//~ 
	//~ // Local remainder	
		//~ // Computes a function which returns the local remainder around a point q.
        //~ void canonical_calculate_local_remainder() {
            //~ canonical_local_remainder = [=] (RVector<D> q) {
				//~ rD_to_r<D> result;
			    //~ rD_to_r<D> local_quadratic = canonical_local_quadratic(q);
//~ 
				//~ return [=](RVector<D> x) {
					//~ return canonical_potential(x) - local_quadratic(x);
				//~ };
			//~ };
		//~ }
//~ 
        //~ // Evaluates the local remainder matrix in a grid.
        //~ template<template<typename ...> class Grid = std::vector>
		//~ Grid<real_t> evaluate_local_remainder_at(Grid<RVector<D>> g, RVector<D> position) {
            //~ return evaluate_function_in_grid<RVector<D>,real_t,Grid,function_t>(canonical_local_remainder(position),g);
        //~ }
//~ };

int main() {
	// Testing
	const int D = 3;
	const int N = 2;
	GMatrix<rD_to_r<D>,N,N> canonical_potential;
	GMatrix<rD_to_rD<D>,N,N> canonical_jacobian;
	GMatrix<rD_to_rDxD<D>,N,N> canonical_hessian;
	
	canonical_potential(0,0) = [=](RVector<D> x){ return x[0]*x[1] - x[2] + 1;};
	canonical_jacobian(0,0) = [=](RVector<D> x) { return RVector<D>({x[1],x[0],-1});};
	canonical_hessian(0,0) = [=](RVector<D> x) { 
		RMatrix<D,D> result;
		result(0,0) = 0;
		result(0,1) = 1;
		result(0,2) = 0;
		result(1,0) = 1;
		result(1,0) = 0;
		result(1,1) = 0;
		result(1,2) = 0;
		result(2,0) = 0;
		result(2,1) = 0;
		result(2,2) = 0;
		return result;
	};
		
	canonical_potential(0,1) = [=](RVector<D> x){ return 2*x[0];};
	canonical_jacobian(0,1) = [=](RVector<D> x) { return RVector<D>({2,0,0});};
	canonical_hessian(0,1) = [=](RVector<D> x) { 
		RMatrix<D,D> result;
		result(0,0) = 0;
		result(0,1) = 0;
		result(0,2) = 0;
		result(1,0) = 0;
		result(1,0) = 0;
		result(1,1) = 0;
		result(1,2) = 0;
		result(2,0) = 0;
		result(2,1) = 0;
		result(2,2) = 0;
		return result;
	};
	
	
	canonical_potential(1,0) = [=](RVector<D> x){ return x[0]*x[0];};
	canonical_jacobian(1,0) = [=](RVector<D> x) { return RVector<D>({2*x[0],0,0});};
	canonical_hessian(1,0) = [=](RVector<D> x) { 
		RMatrix<D,D> result;
		result(0,0) = 2;
		result(0,1) = 0;
		result(0,2) = 0;
		result(1,0) = 0;
		result(1,0) = 0;
		result(1,1) = 0;
		result(1,2) = 0;
		result(2,0) = 0;
		result(2,1) = 0;
		result(2,2) = 0;
		return result;
	};
	canonical_potential(1,1) = [=](RVector<D> x){ return x[1]*x[1];};
	canonical_jacobian(1,1) = [=](RVector<D> x) { return RVector<D>({0,2*x[1],0});};
	canonical_hessian(1,1) = [=](RVector<D> x) { 
		RMatrix<D,D> result;
		result(0,0) = 0;
		result(0,1) = 0;
		result(0,2) = 0;
		result(1,0) = 0;
		result(1,0) = 2;
		result(1,1) = 0;
		result(1,2) = 0;
		result(2,0) = 0;
		result(2,1) = 0;
		result(2,2) = 0;
		return result;
	};	
	
	GVector<rD_to_c<D>,N> eigen_potential;
	GVector<rD_to_cD<D>,N> eigen_jacobian;
	GVector<rD_to_cDxD<D>,N> eigen_hessian;
	eigen_potential(1) = [=](RVector<D> x){ return x[1]*x[1];};
	eigen_jacobian(1) = [=](RVector<D> x) { return CVector<D>({0,2*x[1],0});};
	eigen_hessian(1) = [=](RVector<D> x) { 
		CMatrix<D,D> result;
		result(0,0) = 0;
		result(0,1) = 0;
		result(0,2) = 0;
		result(1,0) = 0;
		result(1,0) = 2;
		result(1,1) = 0;
		result(1,2) = 0;
		result(2,0) = 0;
		result(2,1) = 0;
		result(2,2) = 0;
		return result;
	};	
	
	eigen_potential(0) = [=](RVector<D> x){ return x[0]*x[0];};
	eigen_jacobian(0) = [=](RVector<D> x) { return CVector<D>({2*x[0],0,0});};
	eigen_hessian(0) = [=](RVector<D> x) { 
		CMatrix<D,D> result;
		result(0,0) = 2;
		result(0,1) = 0;
		result(0,2) = 0;
		result(1,0) = 0;
		result(1,0) = 0;
		result(1,1) = 0;
		result(1,2) = 0;
		result(2,0) = 0;
		result(2,1) = 0;
		result(2,2) = 0;
		return result;
	};
	
	rD_to_cNxN<D,N> basis_transform;
	

	
	//~ MatrixPotential<CanonicalBasis,2,3> eigen_test(eigen_potential,  basis_transform, eigen_jacobian, eigen_hessian);
	
	MatrixPotential<CanonicalBasis,2,3> test(canonical_potential, canonical_jacobian, canonical_hessian);
	//~ MatrixPotential<EigenBasis,2,3> testEigs(eigen_potential, eigen_jacobian, eigen_hessian);
	
	//~ test.canonical_evaluate_local_remainder_at<std::vector>(std::vector<RVector<3>>(1,RVector<3>({1,2,3})),RVector<3>{0,1,2});
	
	std::vector<RVector<3>> evalPoints(3);
	evalPoints[0] = RVector<3>{1,1,3};
	evalPoints[1] = RVector<3>{0,2,0};
	evalPoints[2] = RVector<3>{2,1,3};
	test.evaluate_at(evalPoints[0]);
	//~ testEigs.evaluate_at(evalPoints[0]);
	//~ test.test();
	//~ test.evaluate_at<std::vector>(evalPoints);
	//~ test.eigen_evaluate_at(evalPoints);
	//~ test.evaluate_eigenvectors_at(evalPoints);
	//~ test.canonical_evaluate_exponential(evalPoints[0]);
	
	//~ test.evaluate(evalPoints);
	//~ RMatrix<3,2> mat;
	//~ mat.block<3,1>(0,0) = evalPoints[0];
	//~ mat.block<3,1>(0,1) = evalPoints[1];
	//~ 
	//~ //test.evaluate<Helper<2>::type,std::vector>(MatrixToGrid<RMatrix<3,2>>(mat));
	//~ 
	//~ std::cout<<mat<<std::endl;
		
}

