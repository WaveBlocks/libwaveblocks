#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

#include <Eigen/Core>

#include <waveblocks/basic_types.hpp>
#include <waveblocks/wavepackets/hawp_paramset.hpp>

#include <waveblocks/yaml/complex.hpp>
#include <waveblocks/yaml/matrix.hpp>
#include <waveblocks/yaml/hagedornparameters.hpp>


using namespace waveblocks;

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;


    // Create a new yaml document
    std::cout << "Creating a yaml file:" << std::endl;

    const dim_t D = 2;
    wavepackets::HaWpParamSet<D> PI = wavepackets::HaWpParamSet<D>();

    YAML::Node root;
    root["PI"] = PI;

    YAML::Emitter out;
    out.SetIndent(4);
    out.SetSeqFormat(YAML::Flow);
    out << root;

    std::cout << out.c_str() << std::endl;

    std::ofstream fout("pi.yaml");
    fout << out.c_str() << std::endl;


    // Read the file
    std::cout << "Loading a yaml file:" << std::endl;

    YAML::Node config = YAML::LoadFile("pi.yaml");
    wavepackets::HaWpParamSet<D> PIL = config["PI"].as<wavepackets::HaWpParamSet<D>>();

    std::cout << "Loaded values are:\n" << PIL << std::endl;
}
