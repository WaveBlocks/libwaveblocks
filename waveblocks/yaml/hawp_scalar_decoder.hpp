#pragma once

#include <stdexcept>
#include <memory>

#include "hawp_paramset_decoder.hpp"
#include "shape_decoder.hpp"

#include "../hawp_commons.hpp"
#include "../shape_enumerator.hpp"
#include "../csv/hawp_coeffs_loader.hpp"


namespace waveblocks
{
namespace yaml
{

template<dim_t D, class MultiIndex>
struct ScalarHaWpDecoder
{
public:
    unsigned logging_verbosity = 0; // logging disabled by default

    ScalarHaWp<D,MultiIndex> operator()(YAML::Node const& config)
    {
        if (config["dimensionality"] && D != config["dimensionality"].as<int>())
            throw std::runtime_error("incompatible wavepacket. reason: dimensionality does not match");

        if (config["n_components"] && 1 != config["n_components"].as<int>())
            throw std::runtime_error("incompatible wavepacket. reason: number of components > 1");

        if (1 != config["parameters"].size())
            throw std::runtime_error("incompatible wavepacket. reason: number of parameter sets > 1");

        if (1 != config["shapes"].size())
            throw std::runtime_error("incompatible wavepacket. reason: number of shapes > 1");

//         if (1 != config["coefficients"].size())
//             throw std::runtime_error("incompatible wavepacket. reason: number of coefficients vectors > 1");
//

        // (1) decode eps
        ScalarHaWp<D,MultiIndex> wp;
        wp.eps() = config["eps"].as<double>();

        // (2) decode parameters
        HaWpParamSetDecoder<D> paramset_decoder;
        wp.parameters() = paramset_decoder(config["parameters"][0]);

        // (3) decode and enumerate shape
        ShapeDecoder<D> shape_decoder;
        ShapeEnumerator<D,MultiIndex> shape_enumerator;
        AbstractShape<D>* shape = shape_decoder(config["shapes"][0]);
        wp.shape() = shape_enumerator.enumerate(shape);
        delete shape;

        return wp;
    }
};

} // namespace yaml

} // namespace waveblocks
