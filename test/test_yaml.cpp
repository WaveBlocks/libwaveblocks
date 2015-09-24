#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

#include <waveblocks/basic_types.hpp>
#include <waveblocks/yaml/complex.hpp>

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

    const complex_t cn = complex_t(2, 3);

    YAML::Node anode;
    anode["cN"] = cn;

    YAML::Emitter out;
    out << config;
    out << anode;

    std::cout << out.c_str() << std::endl;

    std::ofstream fout("test.yaml");
    fout << out.c_str() << std::endl;
}
