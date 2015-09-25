#include <iostream>
#include <Eigen/Core>
#include <waveblocks/hdf5/eigen3-hdf5.hpp>


void save_matrix(const Eigen::Matrix3d & mat)
{
    H5::H5File file("test_eigen.hdf5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "M3x3", mat);
    file.close();
}


void load_matrix(Eigen::Matrix3d & mat)
{
    H5::H5File file("test_eigen.hdf5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "M3x3", mat);
    file.close();
}


int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;

    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n          ", "[", "]");

    Eigen::Matrix3d M;
    M << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    std::cout << "Matrix M: " << M.format(CleanFmt) << '\n';

    std::cout << "Save matrix M to file: test_eigen.hdf" << std::endl;
    save_matrix(M);

    std::cout << "Load matrix M from file: test_eigen.hdf" << std::endl;
    load_matrix(M);

    std::cout << "Matrix M: " << M.format(CleanFmt) << '\n';
}
