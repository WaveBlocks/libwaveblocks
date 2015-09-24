#include <iostream>
#include <yaml-cpp/yaml.h>

#include <waveblocks/basic_types.hpp>
#include <waveblocks/yaml/complex.hpp>

using namespace waveblocks;

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;

    YAML::Node config = YAML::LoadFile("config.yaml");

    const complex_t c = config["co"].as<complex_t>();
    std::cout << "Complex number c: " << c << std::endl;
}
