#ifndef WAVEBLOCKS_YAML_HAGEDORN_PARAMETER_SET_HPP
#define WAVEBLOCKS_YAML_HAGEDORN_PARAMETER_SET_HPP

#include <stdexcept>

#include <yaml-cpp/yaml.h>
#include <Eigen/Core>

#include "../hawp_paramset.hpp"

#include "complex.hpp"

namespace waveblocks
{
namespace yaml
{

template<dim_t D>
struct HaWpParamSetDecoder
{
public:
    HaWpParamSet<D> operator()(YAML::Node const& config)
    {
        HaWpParamSet<D> params;

        params.q(decode_rvector(config["q"]));
        params.p(decode_rvector(config["p"]));
        params.Q(decode_cmatrix(config["Q"]));
        params.P(decode_cmatrix(config["P"]));

//         if (config["sqrt_detQ"]) {
//             params.sqrt_detQ = config["sqrt_detQ"].as<std::complex<double> >();
//         }
//         
        return params;
    }
    
private:
    Eigen::Matrix<real_t,D,1> decode_rvector(YAML::Node const& config)
    {
        Eigen::Matrix<real_t,D,1> result(D,1);
        
        if (config.size() != D)
            throw std::runtime_error("encountered vector dimension != parameter set dimensionality");
        
        for (dim_t i = 0; i < D; i++) {
            result(i,0) = config[i].as<double>();
        }
        
        return result;
    }
    
    Eigen::Matrix<complex_t,D,D> decode_cmatrix(YAML::Node const& config)
    {
        Eigen::Matrix<complex_t,D,D> result(D,D);
        
        if (config.size() != D)
            throw std::runtime_error("encountered matrix column dimension != parameter set dimensionality");
        
        for (dim_t i = 0; i < D; i++) {
            YAML::Node const& in_row = config[i];
            
            if (in_row.size() != D)
                throw std::runtime_error("encountered matrix row dimension != parameter set dimensionality");
            
            for (dim_t j = 0; j < D; j++) {
                result(i,j) = in_row[j].as< std::complex<double> >();
            }
        }
        
        return result;
    }
};

}
}

#endif