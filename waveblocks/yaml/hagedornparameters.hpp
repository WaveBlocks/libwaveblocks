#pragma once

#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "../types.hpp"
#include "../wavepackets/hawp_paramset.hpp"


namespace YAML {
    template<int D>
    struct convert<waveblocks::wavepackets::HaWpParamSet<D>> {

        static Node encode(const waveblocks::wavepackets::HaWpParamSet<D> & PI) {
            YAML::Node pinode;

            pinode["q"] = PI.q();
            pinode["p"] = PI.p();
            pinode["Q"] = PI.Q();
            pinode["P"] = PI.P();
            pinode["S"] = PI.S();

            return pinode;
        }

        static bool decode(const Node& pinode, waveblocks::wavepackets::HaWpParamSet<D> & PI) {
            if (pinode["q"]) {
                PI.q(pinode["q"].as<waveblocks::RMatrix<D,1>>());
            }
            if (pinode["p"]) {
                PI.p(pinode["p"].as<waveblocks::RMatrix<D,1>>());
            }
            if (pinode["Q"]) {
                PI.Q(pinode["Q"].as<waveblocks::CMatrix<D,D>>());
            }
            if (pinode["P"]) {
                PI.P(pinode["P"].as<waveblocks::CMatrix<D,D>>());
            }
            if (pinode["S"]) {
                PI.S(pinode["S"].as<waveblocks::complex_t>());
            }
            return true;
        }
    };
}
