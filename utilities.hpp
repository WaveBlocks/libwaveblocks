#pragma once
#include <vector>
#include <functional>
#include "types.hpp"

template <class Matrix> class MatrixToGrid;

template <class Matrix> class MatrixToGridIterator {
private:
  MatrixToGrid<Matrix> &g;

  static const int N = Matrix::RowsAtCompileTime;
  static const int M = Matrix::ColsAtCompileTime;
  using grid_element_type = GVector<typename Matrix::Scalar, N>;
  int i;

public:
  MatrixToGridIterator(MatrixToGridIterator &adaptor, int i)
      : g(adaptor), i(i) {}
  bool operator!=(MatrixToGridIterator<Matrix> other) { return other.i != i; }

  void operator++() { ++i; }
  const grid_element_type &operator*() { return g[i]; }
};

template <class Matrix> class MatrixToGrid {
private:
  static const int N = Matrix::RowsAtCompileTime;
  static const int M = Matrix::ColsAtCompileTime;
  using grid_element_type = GVector<typename Matrix::Scalar, N>;

  const Matrix &matrix;

public:
  size_t size() { return M; }

  MatrixToGrid(size_t size) {}
  MatrixToGrid(const Matrix &matrix) : matrix(matrix) {}

  const grid_element_type &operator[](int i) {
    const grid_element_type &x = matrix.template block<N, 1>(0, i);
    return x;
  }

  MatrixToGridIterator<Matrix> begin() {
    return MatrixToGridIterator<Matrix>(*this, 0);
  }
  MatrixToGridIterator<Matrix> end() {
    return MatrixToGridIterator<Matrix>(*this, M);
  }
};

template <int M> struct Helper {
  template <class I>
  using type =
      MatrixToGrid<GMatrix<typename I::Scalar, I::RowsAtCompileTime, M> >;
};

template <class Matrix>
using grid_element_type =
    GVector<typename Matrix::Scalar, Matrix::RowsAtCompileTime>;

// Copies for the time being. Would like to do something cleverer
template <class Matrix, template <typename...> class Grid = std::vector>
Grid<grid_element_type<Matrix> > matrixToGrid(const Matrix &m) {
  const int N = Matrix::RowsAtCompileTime;
  const int M = Matrix::ColsAtCompileTime;

  Grid<grid_element_type<Matrix> > result(M);
  auto it = result.begin();
  for (int i = 0; i < M; ++i) {
    *it = m.template block<N, 1>(0, i);
    ++it;
  }
  return result;
}

template <class Matrix, template <typename...> class Grid = std::vector>
Matrix gridToMatrix(Grid<grid_element_type<Matrix> > g) {
  const int N = Matrix::RowsAtCompileTime;
  const int M = Matrix::ColsAtCompileTime;

  Matrix m;
  auto it = g.begin();
  for (int i = 0; i < M; ++i) {
    m.block<N, 1>(0, i) = *it;
    ++it;
  }

  return m;
}

// Function evaluations
template <class A, class R, template <typename...> class G_in = std::vector,
          template <typename...> class G_out = G_in,
          template <typename...> class F = std::function>
G_out<R> evaluate_function_in_grid(F<R(A)> f, G_in<A> g) {
  G_out<R> result(g.size());
  auto it = result.begin();
  for (auto &arg : g) {
    *it = f(arg);
    ++it;
  }
  return result;
}

template <int N, template <typename, int> class V, class A, class R,
          template <typename...> class F = std::function>
V<R, N> evaluate_function_vector(V<F<R(A)>, N> mf, A arg) {
  V<R, N> m;
  for (int i = 0; i < N; ++i) {
    m[i] = mf[i](arg);
  }

  return m;
}

template <int N, template <typename, int> class V, class A, class R,
          template <typename...> class G_in = std::vector,
          template <typename...> class G_out = G_in,
          template <typename...> class F = std::function>
G_out<V<R, N> > evaluate_function_vector_in_grid(V<F<R(A)>, N> mf, G_in<A> g) {
  G_out<V<R, N> > result(g.size());
  auto it = result.begin();
  for (auto &arg : g) {
    *it = evaluate_function_vector<N, V, A, R, F>(mf, arg);
    ++it;
  }
  return result;
}

template <int N, template <typename, int, int> class M, class A, class R,
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

template <int N, template <typename, int, int> class M, class A, class R,
          template <typename...> class G_in = std::vector,
          template <typename...> class G_out = G_in,
          template <typename...> class F = std::function>
G_out<M<R, N, N> > evaluate_function_matrix_in_grid(M<F<R(A)>, N, N> mf,
                                                    G_in<A> g) {
  G_out<M<R, N, N> > result(g.size());
  auto it = result.begin();
  for (auto &arg : g) {
    *it = evaluate_function_matrix<N, M, A, R, F>(mf, arg);
    ++it;
  }
  return result;
}
