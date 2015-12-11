#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

#include <Eigen/Core>

#include <waveblocks/types.hpp>

#include <waveblocks/yaml/complex.hpp>
#include <waveblocks/yaml/matrix.hpp>

using namespace waveblocks;

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;


    // Create a new yaml document
    std::cout << "Creating a yaml file:" << std::endl;

    RMatrix<1,1> s;
    s << 42;

    RMatrix<3,1> col;
    col << 1,2,3;

    RMatrix<1,3> row;
    row << 1,2,3;

    RMatrix<2,2> R;
    R << 1, 2,
         3, 4;

    CMatrix<3,3> C;
    C << 1.1,            1.2,                1.3,
         complex_t(1,1), complex_t(1.2,3.4), 0,
         3.1,            3.2,                3.3;

    YAML::Node root;

    root["scalar"] = s;
    root["row"] = row;
    root["col"] = col;
    root["R"] = R;
    root["C"] = C;

    YAML::Emitter out;
    out.SetIndent(4);
    //out.SetMapFormat(YAML::Flow);
    out.SetSeqFormat(YAML::Flow);
    out << root;

    std::cout << out.c_str() << std::endl;

    std::ofstream fout("matrices.yaml");
    fout << out.c_str() << std::endl;
}
