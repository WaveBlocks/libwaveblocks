#ifndef WAVEBLOCKS_TEST_SHAPE_ENUMERATION
#define WAVEBLOCKS_TEST_SHAPE_ENUMERATION

#include <iostream>
#include <string>

#include "waveblocks/shape_enumeration_base.hpp"

namespace waveblocks {

template<dim_t D>
void checkShapeEnumeration(const ShapeEnumeration<D> &enumeration, std::string tag)
{
    std::cout << "check " << tag << " {" << std::endl;

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
        for (std::size_t islice = 0; islice < enumeration.count_slices(); islice++) {
            auto & slice = enumeration.slice(islice);

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

}

#endif
