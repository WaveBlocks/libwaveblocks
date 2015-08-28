#include "waveblocks/yaml/shape_decoder.hpp"
#include "waveblocks/yaml/hawp_paramset_decoder.hpp"
#include "waveblocks/yaml/hawp_scalar_decoder.hpp"

#include "yaml-cpp/yaml.h"

#include "waveblocks/tiny_multi_index.hpp"

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
    
    // Load scalar wavepacket
    {
        const dim_t D = 3;
        typedef TinyMultiIndex<std::size_t, D> MultiIndex;
        
        YAML::Node config = YAML::LoadFile("../../test/wavepacket.yaml");
        
        yaml::ScalarHaWpDecoder<D,MultiIndex> wp_decoder;
        wp_decoder.logging_verbosity = 1;
        
        ScalarHaWp<D,MultiIndex> wp = wp_decoder(config, "../../test/");
        
        std::cout << wp.shape()->n_entries() << std::endl;
        std::cout << wp.parameters() << std::endl;
    }
    
    return 0;
}