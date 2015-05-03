#include "waveblocks.hpp"

using namespace waveblocks;

template<dim_t D, class S>
void testShapeExtension(bool strict,
                        const SlicedShapeEnumeration<D,S> &enumeration,
                        const ShapeExtensionEnumeration<D,S> &extension)
{
    std::cout << "check shape extension implementation {" << std::endl;
    for (auto entry : enumeration) {
        //check whether forward neighbours exists
        for (dim_t d = 0; d < D; d++) {
            MultiIndex<D> neighbour = entry;
            neighbour[d] += 1;

            std::size_t extord = extension.find(neighbour);
            if (extord < extension.size()) {
                //extension contains neighbour node
                // => check that shape does NOT
                if (enumeration.find(neighbour) < enumeration.size()) {
                    std::cout << "   [FAILURE] extension contains a node "<<neighbour<<" at "<<extord<<" that is already inside the shape" << std::endl;
                }
            } else {
                //extension does not contain neighbour node
                // => check that shape does contain it
                if (enumeration.find(neighbour) >= enumeration.size()) {
                    std::cout << "   [FAILURE] extension does not contain " << neighbour << std::endl;
                }
            }
        }
    }

    std::size_t nUnneeded = 0;

    //check that for every node inside extension:
    // -> at least one neighbour is contained in shape
    for (auto entry : extension) {
        dim_t count = 0;
        for (dim_t d = 0; d < D; d++) {
            if (entry[d] == 0)
                continue;

            MultiIndex<D> neighbour = entry;
            neighbour[d] -= 1;
            if (enumeration.find(neighbour) < enumeration.size())
                ++count;
        }

        if (count == 0) {
            ++nUnneeded;
            if (strict)
                std::cout << "   [WARNING] extension contains unneeded entry: "<<entry<<std::endl;
        }
    }

    std::cout << "   [INFO] unneeded nodes ratio : "<<nUnneeded<<"/"<<extension.size()<<" = "<<
              double(nUnneeded)/double(extension.size())*100.<<"%"<<std::endl;

    std::cout << "}" << std::endl;
}

int main()
{
    const dim_t D = 30;

    HyperbolicCutShape<D> shape(15.0);

    SlicedShapeEnumeration<D,HyperbolicCutShape<D>> enumeration(shape);

    ShapeExtensionEnumeration<D,HyperbolicCutShape<D>> extension(shape);

    testShapeExtension(false, enumeration, extension);

    return 0;
}
