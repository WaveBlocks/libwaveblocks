#include <Eigen/Core>

#include <iostream>

#include "waveblocks/shape_enum.hpp"
#include "waveblocks/shape_enum_union.hpp"
#include "waveblocks/shape_enum_subset.hpp"
#include "waveblocks/shape_enumerator.hpp"

#include "waveblocks/tiny_multi_index.hpp"

#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/shape_hyperbolic.hpp"
#include "waveblocks/shape_superset.hpp"

#include "waveblocks/hawp.hpp"
#include "waveblocks/hawp_paramset.hpp"
#include "waveblocks/hawp_evaluator.hpp"

using namespace waveblocks;

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    
    bool STATIC_SHAPE_UNION = false;
    bool SIMPLIFIED_EVALUATION = false;
    
    const dim_t D = 5;
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;
    typedef HyperCubicShape<D> S1;
    typedef HyperCubicShape<D> S2;
    typedef LimitedHyperbolicCutShape<D> S3;
    
    S1 shape1({4,4,4,2,2});
    S2 shape2({2,2,4,4,4});
    S3 shape3(9.0, {4,4,4,4,4});
    
    ShapeEnumerator<D,MultiIndex> shape_enumerator;
    
    ShapeEnum<D,MultiIndex> shape_enum1 = shape_enumerator.generate(shape1);
    ShapeEnum<D,MultiIndex> shape_enum2 = shape_enumerator.generate(shape2);
    ShapeEnum<D,MultiIndex> shape_enum3 = shape_enumerator.generate(shape3);
    
    std::array< ShapeEnum<D,MultiIndex>*, 3 > shape_enum_list;
    shape_enum_list[0] = &shape_enum1;
    shape_enum_list[1] = &shape_enum2;
    shape_enum_list[2] = &shape_enum3;
    
    ShapeEnum<D,MultiIndex> enum_union;
    
    if (STATIC_SHAPE_UNION) {
        // This approach generates the enumeration using a description
        // of the union shape. You need to know the exact shape type of
        // all wavepacket components at compile-type.
        SupersetShape<D, S1, S2, S3> shape_union(shape1, shape2, shape3);
        enum_union = shape_enumerator.generate(shape_union);
    } else {
        // This approach merges all sub-enumerations to create the union.
        // The advantage of this approach is, that you dont need to know the exact shape type
        // of all wavepacket components neither.
        enum_union = shape_enum::strict_union(shape_enum_list);
    }
    
    double eps = 0.9;
    HaWpParamSet<D> paramset;
    // paramset = ...
    
    // define quadrature points
    int npts = 5; // number of quadrature points
    const int N = Eigen::Dynamic; // number of quadrature points (template parameter)
    
    Eigen::Matrix<complex_t, D, N> x(D,npts);
    // x = ...
    
    // set up evaluator
    HaWpEvaluator<D,MultiIndex,N> evaluator = hawp::basis(eps, &paramset, &enum_union).at(x);
    
    if (SIMPLIFIED_EVALUATION) {
        // This approach computes all bases at once.
        // Be aware that this approach consumes a lot of memory!
        
        
        // evaluate complete basis for all wavepacket components (using union of all individual enumerations)
        HaWpBasisVector<N> union_basis = evaluator.all();
        
        for (std::size_t c = 0; c < shape_enum_list.size(); c++) {
            ShapeEnum<D,MultiIndex>* enum_comp = shape_enum_list[c];
            
            // gather complete basis for a specific wavepacket component
            HaWpBasisVector<N> comp_basis = shape_enum::copy_subset(union_basis, &enum_union, enum_comp);
            
            // do something with basis ...
        }
    } else {
        // This approach does not evaluate all bases at once.
        // It evaluates slice by slice. This saves memory.
        
        HaWpBasisVector<N> prev_basis;
        HaWpBasisVector<N> curr_basis = evaluator.seed();
        
        for (int islice = 0; islice < enum_union.n_slices(); islice++) {
            const ShapeSlice<D,MultiIndex>& slice_union = enum_union.slice(islice);
            
            for (std::size_t c = 0; c < shape_enum_list.size(); c++) {
                const ShapeSlice<D,MultiIndex>& slice_comp = shape_enum_list[c]->slice(islice);
                
                HaWpBasisVector<N> comp_basis = shape_enum::copy_subset(curr_basis, slice_union, slice_comp);
                
                // do something with basis ...
            }
            
            HaWpBasisVector<N> next_basis = evaluator.step(islice, prev_basis, curr_basis);
            
            prev_basis = std::move(curr_basis);
            curr_basis = std::move(next_basis);
        }
    }
    
    return 0;
}