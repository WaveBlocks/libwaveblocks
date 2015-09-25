#ifndef YAML_MATRIX_HPP
#define YAML_MATRIX_HPP

#include <yaml-cpp/yaml.h>


namespace YAML {
    template<class S, int R, int C>
    struct convert<Eigen::Matrix<S,R,C>> {

        static Node encode(const Eigen::Matrix<S,R,C> & M) {
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

        static bool decode(const Node& matrixnode, Eigen::Matrix<S,R,C> & M) {
            if (matrixnode.size() != M.rows()) throw std::runtime_error("Matrix row size mismatch.");
            for (int r = 0; r < M.rows(); r++) {
                YAML::Node const& rownode = matrixnode[r];
                if (rownode.size() != M.cols()) throw std::runtime_error("Matrix col size mismatch.");
                for (int c = 0; c < M.cols(); c++) {
                    M(r,c) = rownode[c].as<S>();
                }
            }
        }
    };
}

#endif
