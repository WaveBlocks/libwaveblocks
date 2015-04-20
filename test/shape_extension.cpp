#include "waveblocks/shape_extension.hpp"

#include "waveblocks/hyperbolic_shape.hpp"

using namespace waveblocks;

int main()
{
    const dim_t D = 2;
    
    HyperbolicCutShape<D> shape(10);
    
    ShapeExtensionEnumeration<D,HyperbolicCutShape<D>> extension(shape);
    
    
    
    return 0;
}