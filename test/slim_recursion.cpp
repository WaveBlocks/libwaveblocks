#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include "util/time.hpp"

#include "waveblocks.hpp"
#include "sample_wavefunc.hpp"

using namespace waveblocks;

template<dim_t D, class S>
void testSlicedShapeEnumeration(const SlicedShapeEnumeration<D,S> &enumeration)
{
    std::cout << "check enumeration {" << std::endl;

    {
        std::size_t ordinal = 0;
        for (auto index : enumeration) {
            auto ifound = enumeration.find(index);
            if (ordinal != ifound) {
                std::cout << "   [FAILURE] find("<<index<<") = "<<ifound<<" != "<<ordinal << std::endl;
            }

            if (index != enumeration[ordinal]) {
                std::cout << "   [FAILURE] at("<<ordinal<<") != "<<index << std::endl;
            }
            ordinal++;
        }

        if (ordinal != enumeration.size()) {
            std::cout << "   [FAILURE] size() != "<<ordinal << std::endl;
        }
    }

    {
        std::size_t ordinal = 0;
        for (std::size_t islice = 0; islice < enumeration.slices().count(); islice++) {
            auto slice = enumeration.slice(islice);

            if (ordinal != slice.offset()) {
                std::cout << "   [FAILURE] slice_"<<islice<<".offset() != "<< ordinal << std::endl;
            }

            std::size_t ientry = 0;
            for (auto index : slice) {
                auto ifound = slice.find(index);
                if (ientry != ifound) {
                    std::cout << "   [FAILURE] slice_"<<islice<<".find("<<index<< ") = "<<ifound<<" != "<<ientry << std::endl;
                }

                if (index != slice[ientry]) {
                    std::cout << "   [FAILURE] slice_"<<islice<<".at("<<ientry<< ") != "<<index << std::endl;
                }

                if (index != enumeration[ordinal]) {
                    std::cout << "   [FAILURE] order of slice_"<<islice<<" iteration != order of full iteration" << std::endl;
                }

                ++ientry;
                ++ordinal;
            }

            if (ientry != slice.size()) {
                std::cout << "   [FAILURE] slice_"<<islice<<".size() != "<<ientry << std::endl;
            }
        }
    }

    std::cout << "}" << std::endl;
}

template<dim_t D>
MultiIndex<D> createFilledMultiIndex(int entry)
{
    MultiIndex<D> index;
    for (dim_t i = 0; i < D; i++)
        index[i] = entry;
    return index;
}

int main(int argc, char *argv[])
{
    const dim_t D = 3;
    typedef HyperCubicShape<D> S;
    
    S shape(createFilledMultiIndex<D>(15));
    
    auto parameters = createSampleParameters<D>();
    auto enumeration = std::make_shared< SlicedShapeEnumeration<D,S> >(shape);
    auto coefficients = createSampleCoefficients<D,S>(enumeration);
    
    HagedornWavepacket<D,S> wavepacket(parameters, coefficients);

    // evaluate wavepacket at a chosen location
    {
        std::cout << "chosen evaluation {" << std::endl;

        double start = getRealTime();
        Eigen::Matrix<real_t,D,1> x;
        for (dim_t d = 0; d < D; d++)
            x(d,0) = (d+1)/real_t(2*D);
        double stop = getRealTime();

        complex_t psi = wavepacket[x];

        std::cout << "   psi: " << psi << '\n';
        std::cout << "   time: " << (stop - start) << '\n';
        std::cout << "}" << std::endl;
    }

    /**
     * read space delimited csv files with columns x_1 x_2 .. x_d real(psi[x]) imag(psi[x])
     */
    std::cout << "compare results to reference files {" << std::endl;
    for (int iarg = 1; iarg < argc; iarg++) {
        std::cout << "   [FILE] " << argv[iarg] << std::endl;
        std::ifstream in(argv[iarg]);

        std::size_t lines = 0;
        while (in.good()) {
            ++lines;

            //read position
            Eigen::Matrix<real_t,D,1> x;
            for (dim_t d = 0; d < D; d++)
                in >> x(d,0);

            //read reference value
            real_t ref_real, ref_imag;
            in >> ref_real;
            in >> ref_imag;
            complex_t ref(ref_real, ref_imag);

            //compute wavepacket value
            complex_t psi = wavepacket[x];

            auto error = std::norm(psi - ref)/std::norm(ref);

            if (error > 10e-10) {
                std::cout << "      [FAILURE] mismatch at line " << lines << ". error = " << error << std::endl;
            }
        }

        std::cout << "      [INFO] processed " << lines << " lines" << std::endl;
    }
    std::cout << "}" << std::endl;

    return 0;
}
