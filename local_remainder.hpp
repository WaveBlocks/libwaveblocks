#pragma once
#include "macros.hpp"

template <template <template <int, int> class, int, int> class S,
          template <int, int> class B,
          int N,
          int D>
struct LocalRemainderAbstract {

  GET_TYPES(LocalRemainderAbstract, S, B, N, D);

  void calculate_local_quadratic() {
    that.calculate_local_quadratic_implementation();
  }

  void calculate_local_remainder() {
    that.calculate_local_quadratic_implementation();
  }
};

template <template <template <int, int> class, int, int> class S,
          template <int, int> class B,
          int N,
          int D>
class LocalRemainderImplementation : LocalRemainderAbstract<S, B, N, D> {
 public:
  GET_TYPES(LocalRemainderImplementation, S, B, N, D);

 protected:
  local_remainder_type local_remainder;
  local_quadratic_type local_quadratic;

 public:
  LocalRemainderImplementation() {
    calculate_local_quadratic_implementation();
    calculate_local_remainder_implementation();
  }

  void calculate_local_quadratic_implementation() {
    this->local_quadratic = [=](RVector<D> q) {
      potential_type result_matrix;

      for (int l = 0; l < N; ++l) {
        for (int m = 0; m < N; ++m) {
          result_matrix(l, m) = [=](RVector<D> x) {
            auto V = that.potential(l, m)(q);
            auto J = that.jacobian(l, m)(q);
            auto H = that.hessian(l, m)(q);

            auto result = V;

            for (int i = 0; i < D; ++i) {
              auto xmqi = x[i] - q[i];
              result += J[i] * (xmqi);

              for (int j = 0; j < D; ++j) {
                result += 0.5 * xmqi * H(i, j) * (x[j] - q[j]);
              }
            }

            return result;
          };
        }
      }

      return result_matrix;
    };
  }

  void calculate_local_remainder_implementation() {
    this->local_remainder = [=](RVector<D> q) {
      potential_type result;
      auto local_quadratic_q = this->local_quadratic(q);

      for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
          result(i, j) = [=](RVector<D> x) {
            return that.potential(i, j)(x) - local_quadratic_q(i, j)(x);
          };
        }
      }

      return result;
    };
  }
};
