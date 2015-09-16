#include <vector>
#include <iostream>
#include "types.hpp"
#include "bases.hpp"
#include "potentials.hpp"
   
using namespace matrixPotentials;  
int main()
{
  // Testing
  const int D = 3;
  const int N = 2;
  GMatrix<rD_to_r<D>, N, N> canonical_potential;
  GMatrix<rD_to_rD<D>, N, N> canonical_jacobian;
  GMatrix<rD_to_rDxD<D>, N, N> canonical_hessian;
  
  canonical_potential( 0,
  0 ) = [ = ]( RVector<D> x ) {
    return x[0] * x[1] - x[2] + 1;
  };
  canonical_jacobian(
  0, 0 ) = [ = ]( RVector<D> x ) {
    return RVector<D>( {x[1], x[0], -1} );
  };
  canonical_hessian( 0, 0 ) = [ = ]( RVector<D> x ) {
    RMatrix<D, D> result;
    result( 0, 0 ) = 0;
    result( 0, 1 ) = 1;
    result( 0, 2 ) = 0;
    result( 1, 0 ) = 1;
    result( 1, 0 ) = 0;
    result( 1, 1 ) = 0;
    result( 1, 2 ) = 0;
    result( 2, 0 ) = 0;
    result( 2, 1 ) = 0;
    result( 2, 2 ) = 0;
    return result;
  };
  
  canonical_potential( 0, 1 ) = [ = ]( RVector<D> x ) {
    return 2 * x[0];
  };
  canonical_jacobian( 0,
  1 ) = [ = ]( RVector<D> x ) {
    return RVector<D>( {2, 0, 0} );
  };
  canonical_hessian( 0, 1 ) = [ = ]( RVector<D> x ) {
    RMatrix<D, D> result;
    result( 0, 0 ) = 0;
    result( 0, 1 ) = 0;
    result( 0, 2 ) = 0;
    result( 1, 0 ) = 0;
    result( 1, 0 ) = 0;
    result( 1, 1 ) = 0;
    result( 1, 2 ) = 0;
    result( 2, 0 ) = 0;
    result( 2, 1 ) = 0;
    result( 2, 2 ) = 0;
    return result;
  };
  
  canonical_potential( 1, 0 ) = [ = ]( RVector<D> x ) {
    return x[0] * x[0];
  };
  canonical_jacobian(
  1, 0 ) = [ = ]( RVector<D> x ) {
    return RVector<D>( {2 * x[0], 0, 0} );
  };
  canonical_hessian( 1, 0 ) = [ = ]( RVector<D> x ) {
    RMatrix<D, D> result;
    result( 0, 0 ) = 2;
    result( 0, 1 ) = 0;
    result( 0, 2 ) = 0;
    result( 1, 0 ) = 0;
    result( 1, 0 ) = 0;
    result( 1, 1 ) = 0;
    result( 1, 2 ) = 0;
    result( 2, 0 ) = 0;
    result( 2, 1 ) = 0;
    result( 2, 2 ) = 0;
    return result;
  };
  canonical_potential( 1, 1 ) = [ = ]( RVector<D> x ) {
    return x[1] * x[1];
  };
  canonical_jacobian(
  1, 1 ) = [ = ]( RVector<D> x ) {
    return RVector<D>( {0, 2 * x[1], 0} );
  };
  canonical_hessian( 1, 1 ) = [ = ]( RVector<D> x ) {
    RMatrix<D, D> result;
    result( 0, 0 ) = 0;
    result( 0, 1 ) = 0;
    result( 0, 2 ) = 0;
    result( 1, 0 ) = 0;
    result( 1, 0 ) = 2;
    result( 1, 1 ) = 0;
    result( 1, 2 ) = 0;
    result( 2, 0 ) = 0;
    result( 2, 1 ) = 0;
    result( 2, 2 ) = 0;
    return result;
  };
  
  GVector<rD_to_c<D>, N> eigen_potential;
  GVector<rD_to_cD<D>, N> eigen_jacobian;
  GVector<rD_to_cDxD<D>, N> eigen_hessian;
  eigen_potential( 1 ) = [ = ]( RVector<D> x ) {
    return x[1] * x[1];
  };
  eigen_jacobian( 1 ) = [ = ]( RVector<D> x ) {
    return CVector<D>( {0, 2 * x[1], 0} );
  };
  eigen_hessian( 1 ) = [ = ]( RVector<D> x ) {
    CMatrix<D, D> result;
    result( 0, 0 ) = 0;
    result( 0, 1 ) = 0;
    result( 0, 2 ) = 0;
    result( 1, 0 ) = 0;
    result( 1, 0 ) = 2;
    result( 1, 1 ) = 0;
    result( 1, 2 ) = 0;
    result( 2, 0 ) = 0;
    result( 2, 1 ) = 0;
    result( 2, 2 ) = 0;
    return result;
  };
  
  eigen_potential( 0 ) = [ = ]( RVector<D> x ) {
    return x[0] * x[0];
  };
  eigen_jacobian( 0 ) = [ = ]( RVector<D> x ) {
    return CVector<D>( {2 * x[0], 0, 0} );
  };
  eigen_hessian( 0 ) = [ = ]( RVector<D> x ) {
    CMatrix<D, D> result;
    result( 0, 0 ) = 2;
    result( 0, 1 ) = 0;
    result( 0, 2 ) = 0;
    result( 1, 0 ) = 0;
    result( 1, 0 ) = 0;
    result( 1, 1 ) = 0;
    result( 1, 2 ) = 0;
    result( 2, 0 ) = 0;
    result( 2, 1 ) = 0;
    result( 2, 2 ) = 0;
    return result;
  };
  
  rD_to_c<D> leading_level = [=] (RVector<D> x) {
    return x[0]*x[1]+x[2];
  };
  
  rD_to_cD<D> leading_jac = [=] ( RVector<D> x ) {
      return CVector<D>( {2 * x[0], 0, 0} );
  };
  rD_to_cDxD<D> leading_hess = [=] (RVector<D> x) {
        CMatrix<D, D> result;
    result( 0, 0 ) = 2;
    result( 0, 1 ) = 0;
    result( 0, 2 ) = 0;
    result( 1, 0 ) = 0;
    result( 1, 0 ) = 0;
    result( 1, 1 ) = 0;
    result( 1, 2 ) = 0;
    result( 2, 0 ) = 0;
    result( 2, 1 ) = 0;
    result( 2, 2 ) = 0;
    return result;
  };
  
  MatrixPotential<N,D> s(canonical_potential,canonical_jacobian,canonical_hessian,leading_level,leading_jac,leading_hess);
  std::cout<<s.evaluate_local_remainder_at(RVector<D>({1,2,3}),RVector<D>({1,2,3}))<<std::endl;
}
