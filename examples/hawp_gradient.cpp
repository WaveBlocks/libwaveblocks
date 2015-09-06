
#include <iostream>

#include <Eigen/Core>

#include <waveblocks/hawp_commons.hpp>
#include <waveblocks/hawp_gradient_operator.hpp>
#include <waveblocks/tiny_multi_index.hpp>
#include <waveblocks/yaml/hawp_scalar_decoder.hpp>

#include <waveblocks/util/timer.hpp>

#include <yaml-cpp/yaml.h>

#include <boost/format.hpp>

using namespace waveblocks;

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    
    const dim_t D = 3;
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;
    
    if (argc != 3) {
        std::cout << "Program requires 2 arguments " << std::endl;
        std::cout << "  1. argument: wavepacket file (YAML)" << std::endl;
        std::cout << "  2. argument: coefficients file (CSV)" << std::endl;
        return -1;
    }
    
    std::string wp_file = argv[1];
    std::string coeffs_file = argv[2];
    
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "     LOAD WAVEPACKET" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    
    ScalarHaWp<D,MultiIndex> wp;
    try {
        YAML::Node config = YAML::LoadFile(wp_file);
        
        yaml::ScalarHaWpDecoder<D,MultiIndex> wp_decoder;
        wp = wp_decoder(config);
        
        csv::HaWpCoefficientsLoader<D,MultiIndex> coeffs_loader;
        wp.coefficients() = coeffs_loader(wp.shape(), coeffs_file);
        
    } catch (YAML::Exception & error) {
        std::cout << "YAML::Exception: " << error.what() << std::endl;
        return -1;
    } catch (std::runtime_error & error) {
        std::cout << "std::runtime_error: " << error.what() << std::endl;
        return -1;
    }
    
    std::cout << "eps: " << wp.eps() << std::endl;
    std::cout << "shape: " << wp.shape()->n_entries() << " entries, " << wp.shape()->n_slices() << " slices" << std::endl;
    std::cout << wp.parameters() << std::endl;
    
    std::cout << "coefficients {" << std::endl;
    for (int i = 0; i < wp.shape()->n_slices(); i++) {
        auto const& slice = wp.shape()->slice(i);
        
        for (std::size_t j = 0; j < slice.size(); j++) {
            std::cout << "   " << slice[j] << ": " << wp.coefficients()[slice.offset() + j] << std::endl;
        }
    }
    std::cout << "}" << std::endl;
    
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "     EVALUATE WAVEPACKET & GRADIENT" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n   ", "[", "]");
    
    const int numQ = 1;
    RMatrix<D,numQ> grid(D,1);
    std::cout << boost::format("Define %i quadrature points (columns)") % grid.cols() << std::endl;
    for (int i = 0; i < D; i++) {
        grid(i,0) = -1.0 + 2.0*(i)/(D-1);
    }
    
    std::cout << "   " << grid.format(CleanFmt) << std::endl;
    std::cout << std::endl;
    
    std::cout << boost::format("Evaluate basis functions on %i quadrature points") % grid.cols() << std::endl;
    {
        CMatrix<Eigen::Dynamic, numQ> result; // size = ( number of basis functions X number of quadrature points)
        
        result = wp.evaluate_basis(grid);
        
        for (int i = 0; i < wp.shape()->n_slices(); i++) {
            auto const& slice = wp.shape()->slice(i);
            
            for (std::size_t j = 0; j < slice.size(); j++) {
                std::cout << "   " << slice[j] << ": " << result.row(slice.offset() + j) << std::endl;
            }
        }
    }
    std::cout << std::endl;
    
    std::cout << boost::format("Evaluate wavepacket on %i quadrature points") % grid.cols() << std::endl;
    Timer timer;
    {
        CMatrix<1,numQ> result;
        
        timer.start();
        result = wp.evaluate(grid);
        timer.stop();
        std::cout << "   " << result.format(CleanFmt) << std::endl;
        std::cout << "   time: " << timer.millis() << " [ms] " << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << boost::format("Evaluate gradient (%i components) on %i quadrature points") % D % grid.cols() << std::endl;
    {
        HaWpGradientOperator<D,MultiIndex> nabla;
        HaWpGradient<D,MultiIndex> gradwp = nabla(wp);
        
        std::cout << "   Evaluate each component one by one ... " << std::endl;
        double cummtime = 0.0;
        for (std::size_t c = 0; c < gradwp.n_components(); c++) {
            timer.start();
            CMatrix<1,numQ> result = gradwp[c].evaluate(grid);
            timer.stop();
            cummtime += timer.millis();
            
            std::cout << "   " << result.format(CleanFmt) << std::endl;
        }
        std::cout << "   time: " << cummtime << " [ms] " << std::endl;
        std::cout << std::endl;
        
        std::cout << "   Evaluate all components at once ... " << std::endl;
        timer.start();
        CMatrix<Eigen::Dynamic,numQ> result = gradwp.evaluate(grid);
        timer.stop();
        std::cout << "   " << result.format(CleanFmt) << std::endl;
        std::cout << "   time: " << timer.millis() << " [ms] " << std::endl;
    }
    std::cout << std::endl;
}