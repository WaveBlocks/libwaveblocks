#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

#include <Eigen/Core>

#include <waveblocks/basic_types.hpp>
#include <waveblocks/hawp_paramset.hpp>

#include <waveblocks/yaml/complex.hpp>
#include <waveblocks/yaml/matrix.hpp>
#include <waveblocks/yaml/hagedornparameters.hpp>

#include <waveblocks/hdf5/eigen3-hdf5.hpp>


using namespace waveblocks;


int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;

    const dim_t D = 2;

    // Read the yaml file
    std::cout << "Loading a yaml file:" << std::endl;

    YAML::Node config = YAML::LoadFile("pi.yaml");
    HaWpParamSet<D> PI = config["PI"].as<HaWpParamSet<D>>();

    std::cout << "Loaded values are:\n" << PI << std::endl;

    // Save values into a hdf5 file
    H5::H5File file("pi.hdf5", H5F_ACC_TRUNC);

    // TODO: Better API for saving grouped stuff
    auto group = H5::Group(file.createGroup("/PI"));

    EigenHDF5::save(group, "q", PI.q);
    EigenHDF5::save(group, "p", PI.p);
    EigenHDF5::save(group, "Q", PI.Q);
    EigenHDF5::save(group, "P", PI.P);

    file.close();
}
