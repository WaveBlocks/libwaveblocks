#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>

#include <yaml-cpp/yaml.h>

#include "../wavepackets/shapes/shape_commons.hpp"


namespace waveblocks {
    namespace yaml {
        template<dim_t D>
        struct ShapeDecoder
        {
            AbstractShape<D>* operator()(YAML::Node const& config)
            {
                if (config["dimensionality"].as<int>() != D) {
                    throw std::runtime_error("encountered shape dimensionality != requested shape dimensionality");
                }

                std::string type = config["type"].as<std::string>();

                if (type == "hyperbolic") {
                    int sparsity = int(config["sparsity"].as<double>()+0.5);

                    if (config["limits"]) {
                        return new LimitedHyperbolicCutShape<D>(sparsity, decode_limits(config["limits"]));
                    }
                    else {
                        return new HyperbolicCutShape<D>(sparsity);
                    }
                }
                else if (type == "hypercubic") {
                    return new HyperCubicShape<D>(decode_limits(config["limits"]));
                }

                throw std::runtime_error("unknown shape type");
            }

        private:
            std::array<int,D> decode_limits(YAML::Node const& node)
            {
                std::array<int,D> limits;

                if (node.size() != D) {
                    throw std::runtime_error("size of limits != shape dimensionality");
                }

                for (std::size_t i = 0; i < node.size(); i++) {
                    limits[i] = node[i].as<int>();
                }

                return limits;
            }
        };
    }
}
