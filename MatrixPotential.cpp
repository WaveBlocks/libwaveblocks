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
#include"evaluate.hpp"
#include"local_remainder.hpp"
#include"macros.hpp"


template<template<int,int> class B, int N, int D>
struct MatrixPotential : public EvaluationImplementation<MatrixPotential,B,N,D>, public LocalRemainderImplementation<MatrixPotential,B,N,D>{

	using Basis = B<N,D>;
	
	/* ::MatrixPotential is not a typo. Refering to a class template
	 * inside of that is interpreted as refering to the instantiation by
	 * some compilers such as clang
	 * 
	 * this workaround for compilers that haven't fully implemented 
	 * the C++11 injected-class-name rule can be found at
	 * http://stackoverflow.com/a/17687568/3139931
	 */
	friend EvaluationImplementation< ::MatrixPotential,B,N,D>;
	friend LocalRemainderImplementation< ::MatrixPotential,B,N,D>;

	MatrixPotential(typename Basis::potential_type potential,
		typename Basis::jacobian_type jacobian,
		typename Basis::hessian_type hessian) 
		:
		EvaluationImplementation< ::MatrixPotential,B,N,D>(potential,jacobian, hessian),
		LocalRemainderImplementation< ::MatrixPotential,B,N,D>()
	{
	}
	
};



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
	

	
	//~ MatrixPotential<CanonicalBasisImplementation,2,3> eigen_test(eigen_potential,  basis_transform, eigen_jacobian, eigen_hessian);
	
	MatrixPotential<CanonicalBasis,2,3> test(canonical_potential, canonical_jacobian, canonical_hessian);
	//~ MatrixPotential<EigenBasis,2,3> testEigs(eigen_potential, eigen_jacobian, eigen_hessian);
	
	//~ test.canonical_evaluate_local_remainder_at<std::vector>(std::vector<RVector<3>>(1,RVector<3>({1,2,3})),RVector<3>{0,1,2});
	
	std::vector<RVector<3>> evalPoints(3);
	evalPoints[0] = RVector<3>{1,1,3};
	evalPoints[1] = RVector<3>{0,2,0};
	evalPoints[2] = RVector<3>{2,1,3};
	test.evaluate_at(evalPoints[0]);
	test.evaluate(evalPoints);
	//~ CanonicalBasisImplementation<N,D>::evaluate_at(test,evalPoints[0]);
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

