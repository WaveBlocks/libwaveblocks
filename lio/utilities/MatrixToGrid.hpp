#pragma once
#include "types.hpp"

namespace lio
{
  namespace utilities
  {
    template <class Matrix>
    using grid_element_type =
      GVector<typename Matrix::Scalar, Matrix::RowsAtCompileTime>;
      
    template <class Matrix>
    class MatrixToGrid;
    
     /**
     * \brief Forward iterator for the MatrixToGrid class
     * \tparam Matrix The class that the MatrixToGrid is adapting
     */
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
    
    /**
     * \brief Adaptor which accepts an Eigen::Matrix and emulates some behavior of a std::vector
     * 
     * In particular allows to use for:in loops to iterate over each column of the matrix 
     * as well as random accesses to each column.
     * 
     * \tparam Matrix Eigen::Matrix to which to adapt
     */
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
        
        // "grid class template"
        template <class I>
        using type = MatrixToGrid<Matrix>;
    };
    
    
      /**
     * \brief Copies a Eigen::Matrix into a grid
     * 
     * Builds a grid of row vectors from a matrix
     * 
     * \param m The matrix to copy into a grid
     * 
     * \return
     * Grid of row vectors
     * 
     * \tparam Matrix The type of matrix to adapt
     * \tparam Grid Class template to use for return
     */
    template <class Matrix, template <typename...> class Grid = std::vector>
    Grid<grid_element_type<Matrix> > matrix_to_grid( const Matrix &m )
    {
      const int N = m.rows();
      const int M = m.cols();
      
      Grid<grid_element_type<Matrix> > result( M );
      auto it = result.begin();
      
      for ( int i = 0; i < M; ++i ) {
        *it = m.template block<N, 1>( 0, i );
        ++it;
      }
      
      return result;
    }
    
      /**
     * \brief Copies a grid into a Eigen::Matrix
     * 
     * Builds a matrix from a grid of row vectors 
     * 
     * \param g Grid of row vectors
     * \return
     * Matrix with these row vectors
     * 
     * \tparam Matrix The type of matrix to use for return
     * \tparam Grid Class template of the grid to adapt
     */    template <class Matrix, template <typename...> class Grid = std::vector>
    Matrix grid_to_matrix( Grid<grid_element_type<Matrix> > g )
    {
      const int N = Matrix::RowsAtCompileTime;
      
      Matrix m;
      auto it = g.begin();
      
      for ( int i = 0; it != g.end(); ++i, ++it ) {
        m.block<N, 1>( 0, i ) = *it;
      }
      
      return m;
    }
  }
}
