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

template<int N, int D>
class MatrixPotential {
    private:
        GMatrix<rD_to_r<D>,N,N> canonical_potential;
        GMatrix<rD_to_rD<D>,N,N> canonical_jacobian;
        GMatrix<rD_to_rDxD<D>,N,N> canonical_hessian;
        
        GVector<rD_to_c<D>,N> eigen_potential;
        rD_to_cNxN<D,N> eigen_basis_transform;
        GVector<rD_to_cD<D>,N> eigen_jacobian;
        GVector<rD_to_cDxD<D>,N> eigen_hessian;
        bool canonicalBasis;
        
        rD_to_function_matrix<D,N,rD_to_r<D>> canonical_local_quadratic;
        rD_to_function_matrix<D,N,rD_to_r<D>> canonical_local_remainder;
        
		rD_to_function_vector<D,N,rD_to_c<D>> eigen_local_quadratic;
		rD_to_function_vector<D,N,rD_to_c<D>> eigen_local_remainder;

        
    public:

        MatrixPotential(GMatrix<rD_to_r<D>,N,N> potential,
						GMatrix<rD_to_rD<D>,N,N> jacobian,
						GMatrix<rD_to_rDxD<D>,N,N> hessian
			) : canonical_potential(potential),
				canonical_jacobian(jacobian),
				canonical_hessian(hessian),
				canonicalBasis(true)
			{
				canonical_calculate_local_quadratic();
				canonical_calculate_local_remainder();	
			};
		MatrixPotential(GVector<rD_to_c<D>,N> potential,
						rD_to_cNxN<D,N> basis_transform,
						GVector<rD_to_cD<D>,N> jacobian,
						GVector<rD_to_cDxD<D>,N> hessian
			) : eigen_potential(potential),
				eigen_jacobian(jacobian),
				eigen_hessian(hessian),
				eigen_basis_transform(basis_transform),
				canonicalBasis(false)
			{
				eigen_calculate_local_quadratic();
				eigen_calculate_local_remainder();	
			};

		// Evaluates the canonical matrix in multiple points
        template<template<typename ...> class Grid = std::vector>
        Grid<RMatrix<N,N>> canonical_evaluate_at(Grid<RVector<D>> g) {
			// std:bind magic 
			return evaluate_function_in_grid<RVector<D>, RMatrix<N,N>, Grid, function_t>(std::bind(&MatrixPotential<N,D>::canonical_evaluate,this,std::placeholders::_1),g);	
        }
        
        // Evaluates the canonical matrix in one point
        RMatrix<N,N> canonical_evaluate(RVector<D> arg) {
			if (canonicalBasis) {
				return evaluate_function_matrix<N,GMatrix,RVector<D>,real_t,function_t>(canonical_potential,arg);
			}
			else {
				auto transform = eigen_basis_transform(arg);
				CMatrix<N,N> eigs;
				for (int d = 0; d < D; ++d) {
					eigs(d,d) = eigen_potential[d](arg);
				}
				return (transform*D*transform.inverse()).real();			
			}
        }
        
        // Evaluates the jacobian of the canonical matrix in multiple points
        template<template<typename ...> class Grid = std::vector>
        Grid<GMatrix<RVector<D>,N,N>>canonical_evaluate_jacobian_at(Grid<RVector<D>> g) {
			if (canonicalBasis) {
				return evaluate_function_matrix_in_grid(canonical_jacobian,g);
			}
        }

        // Evaluates the hessian of the canonical matrix in multiple points
        template<template<typename ...> class Grid = std::vector>
        Grid<GMatrix<RMatrix<D,D>,N,N>> canonical_evaluate_hessian_at(Grid<RVector<D>> g) {
			if (canonicalBasis) {
				return evaluate_function_matrix_in_grid(canonical_hessian, g);
			}
        }
        
		// Evaluates the jacobians of the eigen functions in multiple points
        template<template<typename ...> class Grid = std::vector>
        Grid<GVector<RVector<D>,N>> eigen_evaluate_jacobian_at(Grid<RVector<D>> g) {
			if (!canonicalBasis) {
				return evaluate_function_vector_in_grid(eigen_jacobian,g);
			}
        }

        // Evaluates the hessians of the eigen functions in multiple points
        template<template<typename ...> class Grid = std::vector>
        Grid<GVector<RMatrix<D,D>,N>> eigen_evaluate_hessian_at(Grid<RVector<D>> g) {
			if (!canonicalBasis) {
				return evaluate_function_vector_in_grid(eigen_hessian, g);
			}
        }
        
        // Evaluates the eigen functions in multiple points
		template<template<typename ...> class Grid = std::vector>
        Grid<CVector<N>> eigenvalue_evaluate_at(Grid<RVector<D>> g) {
			if (canonicalBasis) {
				Grid<RMatrix<N,N>> values = canonical_evaluate_at(g);
				Grid<CVector<N>> res(values.size());
				auto it = res.begin();
				for (auto& matrix: values) {
					Eigen::EigenSolver<RMatrix<N,N>> es(matrix,false);
					*it = es.eigenvalues();
					++it;
				}
				return res;
			}
			else {
				return evaluate_function_vector_in_grid<N, GVector, RVector<D>, complex_t, Grid, function_t>(eigen_potential,g);
			}
        }

		// Evaluates the eigen basis transformation in a grid of points
        template<template<typename ...> class Grid = std::vector>
        Grid<CMatrix<N,N>> evaluate_eigenvectors_at(Grid<RVector<D>> g) {
			if (canonicalBasis) {
				Grid<RMatrix<N,N>> values = canonical_evaluate_at(g);
				Grid<CMatrix<N,N>> res(values.size());
				auto it = res.begin();
				for (auto& matrix: values) {
					Eigen::EigenSolver<RMatrix<N,N>> es(matrix);
					*it = es.eigenvectors();
					++it;
				}
				return res;
			}
			else {
				return evaluate_function_in_grid<RVector<D>, CMatrix<N,N>, Grid, function_t>(eigen_basis_transform,g);
			}				
        }

		// TODO
        template<template<typename ...> class Grid = std::vector>
        Grid<RMatrix<N,N>> canonical_evaluate_exponential_at(Grid<RVector<D>> g, real_t factor = 1) {
			function_t<RMatrix<N,N>(RVector<D>)> canonical_expPotential = [=](RVector<D> args){
                // Compute matrix
                RMatrix<N,N> values = evaluate_at(args);
                RMatrix<N,N> result;
                // Compute exponential
                Eigen::MatrixExponential<RMatrix<N,N>> m_exp(factor * values);

                return m_exp.compute(result);
            };
            return evaluate_function_matrix_in_grid<N,GMatrix,RVector<D>,real_t,Grid,function_t>(canonical_expPotential,g);
        }

        // TODO
        template<template<typename ...> class Tuple = std::tuple, template<typename ...> class Grid = std::vector>
        Tuple<> canonical_evaluate_local_quadratic(Grid<RVector<D>> g){
            return Tuple<>(evaluate_eigenvalues_at(g),evaluate_jacobian_at(g),evaluate_hessian_at(g));
        }

		
		// Computes a function which returns the local remainder around a point q.
        void canonical_calculate_local_remainder() {
            canonical_local_remainder = [=] (RVector<D> q) {
				GMatrix<rD_to_r<D>,N,N> result;
			    GMatrix<rD_to_r<D>,N,N> local_quadratic = canonical_local_quadratic(q);
			    
				for (int i = 0; i < N; ++i) {
					for (int j = 0; j < N; ++j) {
						result(i,j) = [=](RVector<D> x) {
							return canonical_potential(i,j)(x) - local_quadratic(i,j)(x);
						};
					}
				}
				return result;
			};
		}

        // Evaluates the local remainder matrix in a grid.
        template<template<typename ...> class Grid = std::vector>
		Grid<RMatrix<N,N>> canonical_evaluate_local_remainder_at(Grid<RVector<D>> g, RVector<D> position) {
            return evaluate_function_matrix_in_grid<N,GMatrix,RVector<D>,real_t,Grid,function_t>(canonical_local_remainder(position),g);
        }

		// Computes a function which returns the quadratic approximation function matrix for a given point q.
        void canonical_calculate_local_quadratic() {
			if (canonicalBasis) {
				canonical_local_quadratic = [=] (RVector<D> q) {
					GMatrix<rD_to_r<D>,N,N> result_matrix;
							
					for (int l = 0; l < N; ++l) {
						for (int m = 0; m < N; ++m) {
							result_matrix(l,m) = [=](RVector<D> x) {
								real_t V = canonical_potential(l,m)(q);
								RVector<D> J =  canonical_jacobian(l,m)(q);
								RMatrix<D,D> H = canonical_hessian(l,m)(q);

								real_t result = V;
								for(int i = 0; i < D; ++i) {
									real_t xmqi = x[i] - q[i];
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
		}
		
		// Computes the local quadratic approximations of the eigenfunctions
		void eigen_calculate_local_quadratic() {
			if (!canonicalBasis) {
				eigen_local_quadratic = [=] (RVector<D> q) {
					GVector<rD_to_c<D>,N> result_vector;
							
					for (int l = 0; l < N; ++l) {
						result_vector[l] = [=](RVector<D> x) {
							complex_t V = eigen_potential(l)(q);
							CVector<D> J =  eigen_jacobian(l)(q);
							CMatrix<D,D> H = eigen_hessian(l)(q);

							complex_t result = V;
							for(int i = 0; i < D; ++i) {
								real_t xmqi = x[i] - q[i];
								result += J[i]*(xmqi);
								for (int j = 0; j < D; ++j) {
									result += 0.5*xmqi*H(i,j) * (x[j]-q[j]);
								}
							}
							return result;
						};
					}
					return result_vector;
				};
			}
		}
		
		// Computes a function which returns the local remainder around a point q.
        void eigen_calculate_local_remainder() {
            eigen_local_remainder = [=] (RVector<D> q) {
				GVector<rD_to_c<D>,N> result;
			    GVector<rD_to_c<D>,N> local_quadratic = eigen_local_quadratic(q);
			    
				for (int i = 0; i < N; ++i) {
					result[i] = [=](RVector<D> x) {
						return eigen_potential(i)(x) - local_quadratic(i)(x);
					};
				}
				return result;
			};
		}
				

};

// Scalar specialization
template<int D>
class MatrixPotential<1,D> {
    private:
        rD_to_r<D> canonical_potential;
        rD_to_rD<D> canonical_jacobian;
        rD_to_rDxD<D> canonical_hessian;
        
        rD_to_function<D,rD_to_r<D>> canonical_local_quadratic;
        rD_to_function<D,rD_to_r<D>> canonical_local_remainder;
        
    public:
        MatrixPotential(rD_to_r<D> potential,
						rD_to_rD<D> jacobian,
						rD_to_rDxD<D>  hessian
			) : canonical_potential(potential),
				canonical_jacobian(jacobian),
				canonical_hessian(hessian)
			{
				canonical_calculate_local_quadratic();
				canonical_calculate_local_remainder();	
			};
		
		// Evaluates the canonical matrix in multiple points
        template<template<typename ...> class Grid = std::vector>
        Grid<real_t> evaluate_at(Grid<RVector<D>> g) {
			return evaluate_function_in_grid<RVector<D>, real_t, Grid, function_t>(canonical_potential,g);	
        }
        
        // Evaluates the canonical matrix in one point
        real_t evaluate(RVector<D> arg) {
			return canonical_potential(arg);
        }
        // Evaluates the jacobian of the canonical matrix in multiple points
        template<template<typename ...> class Grid = std::vector>
        Grid<RVector<D>> evaluate_jacobian_at(Grid<RVector<D>> g) {
			return evaluate_function_in_grid(canonical_jacobian,g);
        }

        // Evaluates the hessian of the canonical matrix in multiple points
        template<template<typename ...> class Grid = std::vector>
        Grid<RMatrix<D,D>> evaluate_hessian_at(Grid<RVector<D>> g) {
			return evaluate_function_in_grid(canonical_hessian, g);
        }
		
		// Computes a function which returns the local remainder around a point q.
        void canonical_calculate_local_remainder() {
            canonical_local_remainder = [=] (RVector<D> q) {
				rD_to_r<D> result;
			    rD_to_r<D> local_quadratic = canonical_local_quadratic(q);

				return [=](RVector<D> x) {
					return canonical_potential(x) - local_quadratic(x);
				};
			};
		}

        // Evaluates the local remainder matrix in a grid.
        template<template<typename ...> class Grid = std::vector>
		Grid<real_t> evaluate_local_remainder_at(Grid<RVector<D>> g, RVector<D> position) {
            return evaluate_function_in_grid<RVector<D>,real_t,Grid,function_t>(canonical_local_remainder(position),g);
        }

		// Computes a function which returns the quadratic approximation function matrix for a given point q.
        void canonical_calculate_local_quadratic() {
			canonical_local_quadratic = [=] (RVector<D> q) {
				return [=](RVector<D> x) {
					real_t V = canonical_potential(q);
					RVector<D> J =  canonical_jacobian(q);
					RMatrix<D,D> H = canonical_hessian(q);

					real_t result = V;
					for(int i = 0; i < D; ++i) {
						real_t xmqi = x[i] - q[i];
						result += J[i]*(xmqi);
						for (int j = 0; j < D; ++j) {
							result += 0.5*xmqi*H(i,j) * (x[j]-q[j]);
						}
					}
					return result;
				};
			};
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
	

	
	MatrixPotential<2,3> eigen_test(eigen_potential,  basis_transform, eigen_jacobian, eigen_hessian);
	
	MatrixPotential<2,3> test(canonical_potential, canonical_jacobian, canonical_hessian);
	
	test.canonical_evaluate_local_remainder_at<std::vector>(std::vector<RVector<3>>(1,RVector<3>({1,2,3})),RVector<3>{0,1,2});
	
	std::vector<RVector<3>> evalPoints(3);
	evalPoints[0] = RVector<3>{1,1,1};
	evalPoints[1] = RVector<3>{0,0,0};
	evalPoints[2] = RVector<3>{2,1,3};
	test.canonical_evaluate_at<std::vector>(evalPoints);
	test.eigenvalue_evaluate_at(evalPoints);
	test.evaluate_eigenvectors_at(evalPoints);
		
}

