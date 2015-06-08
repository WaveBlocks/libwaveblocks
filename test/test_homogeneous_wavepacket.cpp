
#include <memory>
#include <array>

#include <Eigen/Core>

#include "waveblocks/basic_types.hpp"
#include "waveblocks/shape_superset.hpp"
#include "waveblocks/shape_hyperbolic.hpp"
#include "waveblocks/shape_hypercubic.hpp"


#include "waveblocks/tiny_multi_index.hpp"
#include "waveblocks/shape_enumeration_base.hpp"
#include "waveblocks/shape_enumeration_default.hpp"
#include "waveblocks/shape_enumeration_subset.hpp"

#include "waveblocks/hagedorn_wavepacket.hpp"
#include "waveblocks/hagedorn_basis_evaluator.hpp"

#include "sample_wavepacket.hpp"

using namespace waveblocks;

int main(int argc, char* argv[])
{
    const dim_t D = 5;
    
    typedef HyperbolicCutShape<D> S1;
    typedef HyperCubicShape<D> S2;
    
    S1 shape1(7.0);
    S2 shape2(3);
    
    typedef SupersetShape<D,S1,S2> SS;
    
    SS superset(shape1, shape2);
    
    std::shared_ptr< ShapeEnumeration<D> > superset_enum(new DefaultShapeEnumeration<D, TinyMultiIndex<std::size_t, D>, SS>(superset));
    
    std::array< std::shared_ptr< ShapeEnumeration<D> >, 2> subset_enum_list{
        std::shared_ptr< ShapeEnumeration<D> >{new DefaultShapeEnumeration<D, TinyMultiIndex<std::size_t, D>, S1>(shape1)},
        std::shared_ptr< ShapeEnumeration<D> >{new DefaultShapeEnumeration<D, TinyMultiIndex<std::size_t, D>, S2>(shape2)}
    };
    
    for (auto & slice : superset_enum->slices()) {
        for (auto index : slice) {
            std::cout << index << std::endl;
        }
    }
    
    
    HagedornWavepacket<D> wavepacket{0.9, createSampleParameters<D>(), superset_enum, createSampleCoefficients<D>(superset_enum)};
    
    Eigen::Matrix<complex_t,D,1> x;
    for (dim_t d = 0; d < D; d++) {
        x(d,0) = complex_t( (d+1)/real_t(2*D), (D-d)/real_t(2*D) );
    }
    
    const int N = 1;
    
    Evaluator<D,N> evaluator{wavepacket.eps(), wavepacket.parameters(), wavepacket.enumeration(), x};
    
    Evaluator<D,N>::CArrayXN prev_slice(0,N);
    Evaluator<D,N>::CArrayXN curr_slice(0,N);
    Evaluator<D,N>::CArrayXN next_slice = evaluator.seed();
    
    for (std::size_t islice = 0; islice < superset_enum->n_slices()-1; islice++) {
        prev_slice = std::move(curr_slice);
        curr_slice = std::move(next_slice);
        
        next_slice = evaluator.step(islice, prev_slice, curr_slice);
        
        for (std::size_t c = 0; c < subset_enum_list.size(); c++) {
            Evaluator<D,N>::CArrayXN subset = copy_subset<D,complex_t,N>(next_slice, superset_enum->slice(islice+1), subset_enum_list[c]->slice(islice+1));
            //std::cout << subset << std::endl;
        }
    }
    
    return 0;
}