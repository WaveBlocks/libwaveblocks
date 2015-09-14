#pragma once
#include <vector>
#include <functional>
#include "types.hpp"

// Function evaluations
template <class A,
          class R,
          template <typename...> class G_in = std::vector,
          template <typename...> class G_out = G_in,
          template <typename...> class F = std::function>
G_out<R> evaluate_function_in_grid(F<R(A)> f, G_in<A> g) {
  G_out<R> result(g.size());
  auto it = result.begin();

  for (const auto& arg : g) {
    *it = f(arg);
    ++it;
  }

  return result;
}

template <int N,
          template <typename, int> class V,
          class A,
          class R,
          template <typename...> class F = std::function>
V<R, N> evaluate_function_vector(V<F<R(A)>, N> mf, A arg) {
  V<R, N> m;

  for (int i = 0; i < N; ++i) {
    m[i] = mf[i](arg);
  }

  return m;
}

template <int N,
          template <typename, int> class V,
          class A,
          class R,
          template <typename...> class G_in = std::vector,
          template <typename...> class G_out = G_in,
          template <typename...> class F = std::function>
G_out<V<R, N> > evaluate_function_vector_in_grid(V<F<R(A)>, N> mf, G_in<A> g) {
  G_out<V<R, N> > result(g.size());
  auto it = result.begin();

  for (const auto& arg : g) {
    *it = evaluate_function_vector<N, V, A, R, F>(mf, arg);
    ++it;
  }

  return result;
}

template <int N,
          template <typename, int, int> class M,
          class A,
          class R,
          template <typename...> class F = std::function>
M<R, N, N> evaluate_function_matrix(M<F<R(A)>, N, N> mf, A arg) {
  M<R, N, N> m;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      m(i, j) = mf(i, j)(arg);
    }
  }

  return m;
}

template <int N,
          template <typename, int, int> class M,
          class A,
          class R,
          template <typename...> class G_in = std::vector,
          template <typename...> class G_out = G_in,
          template <typename...> class F = std::function>
G_out<M<R, N, N> > evaluate_function_matrix_in_grid(M<F<R(A)>, N, N> mf,
                                                    G_in<A> g) {
  G_out<M<R, N, N> > result(g.size());
  auto it = result.begin();

  for (const auto& arg : g) {
    *it = evaluate_function_matrix<N, M, A, R, F>(mf, arg);
    ++it;
  }

  return result;
}
