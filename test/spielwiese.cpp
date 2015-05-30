#include <Eigen/Core>
#include <iostream>

#include <waveblocks/shape_hyperbolic.hpp>
#include <waveblocks/shape_extended.hpp>
#include <waveblocks/tiny_multi_index.hpp>
#include <waveblocks/lexical_shape_enumerator.hpp>

using namespace waveblocks;

int main()
{
    const dim_t D = 4;
    
    typedef TinyMultiIndex<std::size_t,D> MultiIndex;
    typedef HyperbolicCutShape<D> S;
    typedef ExtendedShape<D,S> ES;
    
    S shape(7.0);
    
    ES extended(shape);
    
    {
        LexicalIndexGenerator<D,MultiIndex,S> gen(shape);
        do {
            std::cout << gen.index() << std::endl;
        } while (gen.forward());
    }
    
    {
        LexicalIndexGenerator<D,MultiIndex,ES> gen(extended);
        do {
            std::cout << gen.index() << std::endl;
        } while (gen.forward());
    }
    
    return 0;
}