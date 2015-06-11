
#include <memory>
#include <array>

#include <Eigen/Core>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/shape_superset.hpp"
#include "waveblocks/shape_hyperbolic.hpp"
#include "waveblocks/shape_hypercubic.hpp"


#include "waveblocks/tiny_multi_index.hpp"
#include "waveblocks/shape_enumeration_subset.hpp"
#include "waveblocks/shape_enumeration_operations.hpp"
#include "waveblocks/shape_enum.hpp"
#include "waveblocks/shape_enumerator.hpp"

#include "waveblocks/hawp_evaluator.hpp"

#include "sample_wavepacket.hpp"

using namespace waveblocks;

int main(int argc, char* argv[])
{
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
    
    ShapeEnum<D,MultiIndex> enum_union = enumerator.generate(superset);
    
    for (auto & slice : enum_union.slices()) {
        for (auto index : slice) {
            std::cout << index << std::endl;
        }
    }
    
    double eps = 0.9;
    
    Eigen::Matrix<complex_t,D,1> x;
    for (dim_t d = 0; d < D; d++) {
        x(d,0) = complex_t( (d+1)/real_t(2*D), (D-d)/real_t(2*D) );
    }
    
    const int N = 1;
    
    HaWpEvaluator<D,MultiIndex,N> evaluator{eps, createSampleParameters<D>(), enum_union, x};
    
    HaWpEvaluator<D,MultiIndex,N>::CArrayXN prev_slice(0,N);
    HaWpEvaluator<D,MultiIndex,N>::CArrayXN curr_slice(0,N);
    HaWpEvaluator<D,MultiIndex,N>::CArrayXN next_slice = evaluator.seed();
    
    for (int islice = 0; islice < enum_union.n_slices()-1; islice++) {
        prev_slice = std::move(curr_slice);
        curr_slice = std::move(next_slice);
        
        next_slice = evaluator.step(islice, prev_slice, curr_slice);
        
        HaWpEvaluator<D,MultiIndex,N>::CArrayXN subset = copy_subset(next_slice, enum_union.slice(islice+1), enum1.slice(islice+1));
    }
    
    return 0;
}