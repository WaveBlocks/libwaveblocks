#pragma once
#include "types.hpp"

namespace lio { 
  namespace utilities {
    template <class Matrix>
    using grid_element_type =
      GVector<typename Matrix::Scalar, Matrix::RowsAtCompileTime>;

    // Adaptor class
    template <class Matrix>
    class MatrixToGrid;

    template <class Matrix>
    class MatrixToGridIterator
    {
      private:
        MatrixToGrid<Matrix> &g;
        
        static const int N = Matrix::RowsAtCompileTime;
        static const int M = Matrix::ColsAtCompileTime;
        using grid_element_type = GVector<typename Matrix::Scalar, N>;
        int i;
        
      public:
        MatrixToGridIterator( MatrixToGrid<Matrix> &adaptor, int i )
          : g( adaptor ), i( i ) {}
        bool operator!=( MatrixToGridIterator<Matrix> other ) {
          return other.i != i;
        }
        
        void operator++() {
          ++i;
        }
        grid_element_type operator*() {
          return g[i];
        }
    };

    template <class Matrix>
    class MatrixToGrid
    {
      private:
        static const int N = Matrix::RowsAtCompileTime;
        static const int M = Matrix::ColsAtCompileTime;
        using grid_element_type = GVector<typename Matrix::Scalar, N>;
        
        const Matrix &matrix;
        
      public:
        size_t size() {
          return M;
        }
        
        MatrixToGrid( const Matrix &matrix ) : matrix( matrix ) {}
        
        grid_element_type operator[]( int i ) {
          return matrix.template block<N, 1>( 0, i );
        }
        
        MatrixToGridIterator<Matrix> begin() {
          return MatrixToGridIterator<Matrix>( *this, 0 );
        }
        MatrixToGridIterator<Matrix> end() {
          return MatrixToGridIterator<Matrix>( *this, M );
        }
        
        template <class I>
        using type = MatrixToGrid<Matrix>;
    };




    // Copiers

    // Copies matrix into grid
    template <class Matrix, template <typename...> class Grid = std::vector>
    Grid<grid_element_type<Matrix> > matrix_to_grid( const Matrix &m )
    {
      const int N = Matrix::RowsAtCompileTime;
      const int M = Matrix::ColsAtCompileTime;
      
      Grid<grid_element_type<Matrix> > result( M );
      auto it = result.begin();
      
      for ( int i = 0; i < M; ++i ) {
        *it = m.template block<N, 1>( 0, i );
        ++it;
      }
      
      return result;
    }

    // copies grid into matrix
    template <class Matrix, template <typename...> class Grid = std::vector>
    Matrix grid_to_matrix( Grid<grid_element_type<Matrix> > g )
    {
      const int N = Matrix::RowsAtCompileTime;
      const int M = Matrix::ColsAtCompileTime;
      
      Matrix m;
      auto it = g.begin();
      
      for ( int i = 0; i < M; ++i ) {
        m.block<N, 1>( 0, i ) = *it;
        ++it;
      }
      
      return m;
    }
  }
}
