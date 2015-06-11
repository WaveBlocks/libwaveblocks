
#include <array>

#include "waveblocks/shape_superset.hpp"
#include "waveblocks/shape_hyperbolic.hpp"
#include "waveblocks/shape_hypercubic.hpp"

#include "waveblocks/tiny_multi_index.hpp"
#include "waveblocks/shape_enum.hpp"
#include "waveblocks/shape_enumerator.hpp"
#include "waveblocks/shape_enum_union.hpp"

using namespace waveblocks;

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;
    
    const dim_t D = 5;
    
    typedef HyperbolicCutShape<D> S1;
    typedef HyperCubicShape<D> S2;
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;
    
    S1 shape1(7.0);
    S2 shape2(3);
    
    typedef SupersetShape<D,S1,S2> SS;
    
    SS superset(shape1, shape2);
    
    ShapeEnumerator<D,MultiIndex> enumerator;
    
    ShapeEnum<D,MultiIndex> enum1 = enumerator.generate(shape1);
    ShapeEnum<D,MultiIndex> enum2 = enumerator.generate(shape2);
    
    ShapeEnum<D,MultiIndex> enum_static_union = enumerator.generate(superset);
    
    ShapeEnum<D,MultiIndex> enum_dynamic_union = shape_enum::strict_union({&enum1, &enum2});
    
    if (enum_static_union.n_entries() != enum_dynamic_union.n_entries()) {
        std::cout << "[FAIL] size of static union != dynamic union" << std::endl;
        return 1;
    }
    
    if (enum_static_union.limits() != enum_dynamic_union.limits()) {
        std::cout << "[FAIL] limits of static union != dynamic union" << std::endl;
        return 1;
    }
    
    for (int islice = 0; islice < std::max(enum_static_union.n_slices(),enum_dynamic_union.n_slices()); islice++) {
        if (enum_static_union.slice(islice).offset() != enum_dynamic_union.slice(islice).offset()) {
            std::cout << "[FAIL] offsets of static union != dynamic union" << std::endl;
            return 1;
        }
    }
    
    if (enum_static_union != enum_dynamic_union) {
        std::cout << "[FAIL] content of static union != dynamic union" << std::endl;
        return 1;
    }
    
    std::cout << "[SUCCESS] test passed" << std::endl;
    
    return 0;
}