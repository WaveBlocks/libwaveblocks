#ifndef YAML_MATRIX_HPP
#define YAML_MATRIX_HPP

#include <yaml-cpp/yaml.h>


namespace YAML {
    template<class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    struct convert<Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>> {

        static Node encode(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & M) {
            YAML::Node matrixnode;
            for (int r = 0; r < M.rows(); r++) {
                YAML::Node rownode;
                for (int c = 0; c < M.cols(); c++) {
                    rownode[c] = M(r,c);
                }
                matrixnode[r] = rownode;
            }
            return matrixnode;
        }

        static bool decode(const Node& matrixnode, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & M) {
            if (matrixnode.size() != (std::size_t)M.rows()) throw std::runtime_error("Matrix row size mismatch.");
            for (int r = 0; r < M.rows(); r++) {
                YAML::Node const& rownode = matrixnode[r];
                if (rownode.size() != (std::size_t)M.cols()) throw std::runtime_error("Matrix col size mismatch.");
                for (int c = 0; c < M.cols(); c++) {
                    M(r,c) = rownode[c].as<_Scalar>();
                }
            }
            return true;
        }
    };
}

#endif
