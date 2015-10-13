#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

#include <Eigen/Core>

#include <waveblocks/basic_types.hpp>
#include <waveblocks/yaml/complex.hpp>

#include <waveblocks/yaml/matrix.hpp>

using namespace waveblocks;

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;

    // Load from an existing yaml document
    std::cout << "Loading a yaml file:" << std::endl;

    YAML::Node config = YAML::LoadFile("config.yaml");
    const complex_t c1 = config["c1"].as<complex_t>();

    std::cout << "Complex number c1: " << c1 << std::endl;


    // Create a new yaml document
    std::cout << "Creating a yaml file:" << std::endl;

    YAML::Node anode;

    const complex_t cn = complex_t(2, 3);
    anode["cN"] = cn;

    CMatrix<3,3> q;
    q <<  1,2,3,
          4,5,6,
          7,8,9;
    anode["q"] = q;

    YAML::Node root;
    root["part0"] = config;
    root["part1"] = anode;

    YAML::Emitter out;
    out.SetIndent(4);
    //out.SetMapStyle(YAML::Flow);
    out << root;

    std::cout << out.c_str() << std::endl;

    std::ofstream fout("test.yaml");
    fout << out.c_str() << std::endl;
}
