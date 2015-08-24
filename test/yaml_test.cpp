#include "waveblocks/yaml/shape.hpp"
#include "waveblocks/yaml/hawp_paramset.hpp"

#include "yaml-cpp/yaml.h"

#include <complex>
#include <sstream>

#include <cstdlib>

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    
    using namespace waveblocks;
    
    // Load shape
    {
        YAML::Node config = YAML::LoadFile("../../test/shape.yaml");
        
        for (std::size_t i = 0; i < config.size(); i++) {
            yaml::ShapeDecoder<3> shape_decoder;
            
            AbstractShape<3>* shape = shape_decoder(config[i]);
            
            std::cout << *shape << std::endl;
            
            delete shape;
        }
    }
    
    // Load parameter set
    {
        YAML::Node config = YAML::LoadFile("../../test/paramset.yaml");
        
        yaml::HaWpParamSetDecoder<3> paramset_decoder;
        
        std::cout << paramset_decoder(config) << std::endl;
    }
    
    return 0;
}