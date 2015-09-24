#ifndef YAML_MATRIX_HPP
#define YAML_MATRIX_HPP

#include <yaml-cpp/yaml.h>

#include "../basic_types.hpp"


namespace YAML {
    template<int D>
    struct convert<Eigen::Matrix<waveblocks::complex_t,D,D> > {

        static Node encode(const Eigen::Matrix<waveblocks::complex_t,D,D> & rhs) {
            YAML::Node node;
            for (int i = 0; i < D; i++) {
                YAML::Node row;
                for (int j = 0; j < D; j++) {
                    row[j] = rhs(i,j);
                }
                node[i] = row;
            }
            return node;
        }

        static bool decode(const Node& node, Eigen::Matrix<waveblocks::complex_t,D,D> & rhs) {
            if (node.size() != D)
                throw std::runtime_error("encountered matrix column dimension != parameter set dimensionality");

            for (waveblocks::dim_t i = 0; i < D; i++) {
                YAML::Node const& row = node[i];

                if (row.size() != D)
                    throw std::runtime_error("encountered matrix row dimension != parameter set dimensionality");

                for (waveblocks::dim_t j = 0; j < D; j++) {
                    rhs(i,j) = row[j].as<waveblocks::complex_t>();
                }
            }
        }
    };
}

#endif
