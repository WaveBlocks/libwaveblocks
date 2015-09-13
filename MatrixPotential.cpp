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

//General definition			
template<template<int, int> class B, int N, int D>
class MatrixPotential {
	
	private:
	using Basis = B<N,D>;
	using potential_type = typename Basis::potential_type;
	using jacobian_type = typename Basis::jacobian_type;
	using hessian_type = typename Basis::hessian_type;
	using local_quadratic_type = typename Basis::local_quadratic_type;
	using local_remainder_type = typename Basis::local_remainder_type;
	
	using potential_evaluation_type = typename Basis::potential_evaluation_type;
	using jacobian_evaluation_type = typename Basis::jacobian_evaluation_type;
	using hessian_evaluation_type = typename Basis::hessian_evaluation_type;
	using self_type = MatrixPotential<B,N,D>;

	friend Basis;

	
	
	
	
	potential_type potential;
	jacobian_type jacobian;
	hessian_type hessian; 
	local_quadratic_type local_quadratic;
	local_remainder_type local_remainder;
	
	public:
	MatrixPotential(potential_type potential,
					jacobian_type jacobian,
					hessian_type hessian
		) : potential(potential),
			jacobian(jacobian),
			hessian(hessian)
		{
			calculate_local_quadratic();
			calculate_local_remainder();
		}
		
	potential_evaluation_type evaluate_at(RVector<D> arg) {
		return Basis::evaluate_at(*this, arg);	
	}
	
	template<template<typename ...> class grid_in = std::vector, template<typename ...> class grid_out = grid_in>
	grid_out<potential_evaluation_type> evaluate(grid_in<RVector<D>> args) {
		return evaluate_function_in_grid<RVector<D>, potential_evaluation_type, grid_in,grid_out,function_t>(std::bind(&MatrixPotential<B,N,D>::evaluate_at,this,std::placeholders::_1),args);
	}
	
	
	
	jacobian_evaluation_type evaluate_jacobian_at(RVector<D> arg) {
		return Basis::evaluate_jacobian_at(*this, arg);
	}
	
	template<template<typename ...> class grid_in = std::vector, template<typename ...> class grid_out = grid_in>
	grid_out<jacobian_evaluation_type> evaluate_jacobian(grid_in<RVector<D>> args) {
		return evaluate_function_in_grid<RVector<D>, jacobian_evaluation_type, grid_in,grid_out,function_t>(std::bind(&MatrixPotential<B,N,D>::evaluate_jacobian_at,this,std::placeholders::_1),args);
	}
	
	
	
	hessian_evaluation_type evaluate_hessian_at(RVector<D> arg) {
		return Basis::evaluate_hessian_at(*this, arg);
	}
	
	template<template<typename ...> class grid_in = std::vector, template<typename ...> class grid_out = grid_in>
	grid_out<hessian_evaluation_type> evaluate_hessian(grid_in<RVector<D>> args) {
		return evaluate_function_in_grid<RVector<D>, hessian_evaluation_type, grid_in,grid_out,function_t>(std::bind(&MatrixPotential<B,N,D>::evaluate_hessian_at,this,std::placeholders::_1),args);
	}
	
	
	
	template<template<typename ...> class Tuple = std::tuple>
	Tuple<> taylor_at(RVector<D> g){
		return Tuple<>(evaluate_at(g), evaluate_jacobian_at(g), evaluate_hessian_at(g));
	}
	
	template<template<typename ...> class Tuple = std::tuple, template<typename ...> class grid_in = std::vector, template<typename ...> class grid_out = grid_in>
	grid_out<Tuple<>> taylor(grid_in<RVector<D>> args) {
		return evaluate_function_in_grid<RVector<D>, Tuple<>, grid_in,grid_out,function_t>(std::bind(&MatrixPotential<B,N,D>::taylor_at,this,std::placeholders::_1),args);
	}
	
	
	void calculate_local_quadratic() {
		local_quadratic = [=] (RVector<D> q) {
			potential_type result_matrix;
					
			for (int l = 0; l < N; ++l) {
				for (int m = 0; m < N; ++m) {
					result_matrix(l,m) = [=](RVector<D> x) {
						auto V = potential(l,m)(q);
						auto J =  jacobian(l,m)(q);
						auto H = hessian(l,m)(q);

						auto result = V;
						for(int i = 0; i < D; ++i) {
							auto xmqi = x[i] - q[i];
							result += J[i]*(xmqi);
							for (int j = 0; j < D; ++j) {
								result += 0.5*xmqi*H(i,j) * (x[j]-q[j]);
							}
						}
						return result;
					};
				}
			}
			return result_matrix;
		};
	}
	
	void calculate_local_remainder() {
		local_remainder = [=] (RVector<D> q) {
			potential_type result;
			auto local_quadratic_q = local_quadratic(q);
			
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					result(i,j) = [=](RVector<D> x) {
						return potential(i,j)(x) - local_quadratic_q(i,j)(x);
					};
				}
			}
			return result;
		};
	}
	
	potential_evaluation_type evaluate_local_remainder_at(RVector<D> g, RVector<D> position) {
		return Basis::evaluate_local_remainder_at(*this,g,position);
	}
	
	template<template<typename ...> class grid_in = std::vector, template<typename ...> class grid_out = grid_in>
	grid_out<potential_evaluation_type> evaluate_local_remainder(grid_in<RVector<D>> args) {
		return evaluate_function_in_grid<RVector<D>, potential_evaluation_type, grid_in,grid_out,function_t>(std::bind(&MatrixPotential<B,N,D>::evaluate_local_remainder_at,this,std::placeholders::_1),args);
	}
	
	potential_evaluation_type evaluate_exponential_at(RVector<D> arg, real_t factor = 1) {
		// Compute matrix
		auto values = evaluate(arg);
		potential_evaluation_type result;
		
		// Compute exponential
		Eigen::MatrixExponential<potential_evaluation_type> m_exp(factor * values);

		m_exp.compute(result);
		return result;
	}
	
	template<template<typename ...> class grid_in = std::vector, template<typename ...> class grid_out = grid_in>
	grid_out<potential_evaluation_type> evaluate_exponential(grid_in<RVector<D>> args, real_t factor = 1) {
		return evaluate_function_in_grid<RVector<D>, potential_evaluation_type, grid_in,grid_out,function_t>(std::bind(&MatrixPotential<B,N,D>::evaluate_exponential_at,this,std::placeholders::_1,factor),args);
	}
};

template<template<int, int> class B1, template<int,int> class B2, int N, int D>
class TransformableMatrixPotential : public MatrixPotential<B1,N,D> {
	private:
	typename B1::transformation_type transformation;
	
	public:
	
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
};



template<int N, int D>
struct CanonicalBasis {
	using potential_type = GMatrix<rD_to_r<D>,N,N>;
	using jacobian_type = GMatrix<rD_to_rD<D>,N,N>;
	using hessian_type = GMatrix<rD_to_rDxD<D>,N,N>;
	using transformation_type = rD_to_cNxN<D,N>;
	using local_quadratic_type = rD_to_function_matrix<D,N,rD_to_r<D>>;
    using local_remainder_type = rD_to_function_matrix<D,N,rD_to_r<D>>;
    
	using potential_evaluation_type = RMatrix<N,N>;
	using jacobian_evaluation_type = GMatrix<RVector<D>,N,N>;
	using hessian_evaluation_type = GMatrix<RMatrix<D,D>,N,N>;

	private:
	friend MatrixPotential< ::CanonicalBasis,N,D>;
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
	
	static real_t evaluate_local_remainder_at(const MatrixPotential< ::CanonicalBasis,N,D>& that, RVector<D> g, RVector<D> position) {
		return evaluate_function_matrix<N,GMatrix,RVector<D>,real_t,function_t>(that.local_remainder(position),g);
	}

};

template<int N, int D>
struct EigenBasis {
	using potential_type = GVector<rD_to_c<D>,N>;
	using jacobian_type = GVector<rD_to_cD<D>,N>;
	using hessian_type = GVector<rD_to_cDxD<D>,N>;
	using transformation_type = rD_to_cNxN<D,N>;
	using local_quadratic_type = rD_to_function_vector<D,N,rD_to_c<D>>;
    using local_remainder_type = rD_to_function_vector<D,N,rD_to_c<D>>;
    using transformation_type = rD_to_cNxN<D,N>;
    
    using potential_evaluation_type = CVector<N>;
    using jacobian_evaluation_type = GVector<CVector<D>,N>;
	using hessian_evaluation_type = GVector<CMatrix<D,D>,N>;

	private:
	friend MatrixPotential< ::CanonicalBasis,N,D>;

	static potential_evaluation_type evaluate_at(const MatrixPotential< ::EigenBasis,N,D>& that, RVector<D> arg) {
		return evaluate_function_vector<N,GVector,RVector<D>,complex_t,function_t>(that.potential,arg);
	}
	
	static jacobian_evaluation_type evaluate_jacobian_at(const MatrixPotential< ::EigenBasis,N,D>& that, RVector<D> arg) {
		return evaluate_function_vector(that.jacobian,arg);
	}
	
	static hessian_evaluation_type evaluate_hessian_at(const MatrixPotential< ::EigenBasis,N,D>& that, RVector<D> arg) {
		return evaluate_function_matrix(that.hessian, arg);
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
	

	
	//~ MatrixPotential<CanonicalBasis,2,3> eigen_test(eigen_potential,  basis_transform, eigen_jacobian, eigen_hessian);
	
	MatrixPotential<CanonicalBasis,2,3> test(canonical_potential, canonical_jacobian, canonical_hessian);
	//~ MatrixPotential<EigenBasis,2,3> testEigs(eigen_potential, eigen_jacobian, eigen_hessian);
	
	//~ test.canonical_evaluate_local_remainder_at<std::vector>(std::vector<RVector<3>>(1,RVector<3>({1,2,3})),RVector<3>{0,1,2});
	
	std::vector<RVector<3>> evalPoints(3);
	evalPoints[0] = RVector<3>{1,1,3};
	evalPoints[1] = RVector<3>{0,2,0};
	evalPoints[2] = RVector<3>{2,1,3};
	test.evaluate_at(evalPoints[0]);
	test.evaluate(evalPoints);
	//~ CanonicalBasis<N,D>::evaluate_at(test,evalPoints[0]);
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

