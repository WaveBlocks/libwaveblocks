#include <array>
#include <string>
#include <iostream>

#include <waveblocks/shape_hypercubic.hpp>
#include <waveblocks/shape_hyperbolic.hpp>
#include <waveblocks/shape_extended.hpp>

#include <waveblocks/shape_enum.hpp>
#include <waveblocks/shape_enum_extended.hpp>
#include <waveblocks/shape_enumerator.hpp>
#include <waveblocks/tiny_multi_index.hpp>

#include <waveblocks/stdarray2stream.hpp>

using namespace waveblocks;

template<dim_t D, class MultiIndex>
bool check_extension(ShapeEnum<D,MultiIndex> const& base_enum, ShapeEnum<D,MultiIndex> const& ext_enum, bool strict, bool verbose)
{
    bool passed = true;
    
    //ensure that for every node inside base shape, the extended shape contains all forward neighbours
    for (int islice = 0; islice < base_enum.n_slices(); islice++) {
        auto & base_slice = base_enum.slice(islice);
        auto & ext_slice = ext_enum.slice(islice+1);
        
        for (auto node : base_slice) {
            for (dim_t d = 0; d < D; d++) {
                MultiIndex fnode = node; fnode[d] += 1;
                
                std::size_t ordinal;
                if (!ext_slice.try_find(fnode, ordinal)) {
                    std::cout << "   [FAILURE] extended shape does not contain " << fnode << std::endl;
                    passed = false;
                }
            }
        }
    }
    
    //count unneeded nodes contained in extended shape
    std::size_t total = 1;
    std::size_t nUnneeded = 0;
    
    for (int islice = 1; islice < ext_enum.n_slices(); islice++) {
        auto & base_slice = base_enum.slice(islice-1);
        auto & ext_slice = ext_enum.slice(islice);
        
        for (auto node : ext_slice) {
            ++total;
            
            bool miss = true;
            for (dim_t d = 0; d < D; d++) {
                if (node[d] > 0) {
                    MultiIndex bnode = node; bnode[d] -= 1;
                    
                    std::size_t ordinal;
                    if (base_slice.try_find(bnode, ordinal)) {
                        miss = false;
                        break;
                    }
                }
            }
            
            if (miss) {
                ++nUnneeded;
                if (strict)
                    passed = false;
                else if (verbose)
                    std::cout << "   [WARNING] extended shape contains unneeded entry: "<<node<<std::endl;
            }
        }
    }
    
    if (!strict) {
        std::cout << "   [INFO] overhead : "<<nUnneeded<<"/"<<total<<" = "<<
        double(nUnneeded)/double(total)*100.<<"%"<<std::endl;
    }
    
    return passed;
}

/**
 * \tparam D dimensions
 * \tparam S type of test shape
 * \param[in] base test shape
 * \param[in] verbose if enabled: print all unneeded nodes of extended shape
 */
template<dim_t D, class S>
bool test(const S &base, bool verbose)
{
    bool passed = true;
    
    
    std::cout << "{" << std::endl;
    
    typedef TinyMultiIndex<std::size_t,D> MultiIndex;
    
    // (1) CREATE EXTENSION AT COMPILE-TIME
    {
        ExtendedShape<D,S> extended(base);
        
        ShapeEnumerator<D,MultiIndex> enumerator;
        
        ShapeEnum<D,MultiIndex> ext_enum = enumerator.generate(extended);
        ShapeEnum<D,MultiIndex> base_enum = enumerator.generate(base);
        
        passed &= check_extension<D,MultiIndex>(base_enum, ext_enum, false, verbose);
    }
    
    // (2) CREATE EXTENSION AT RUNTIME
    {
        ShapeEnumerator<D,MultiIndex> enumerator;
        ShapeEnum<D,MultiIndex> base_enum = enumerator.generate(base);
        ShapeEnum<D,MultiIndex> ext_enum = shape_enum::extend(&base_enum);
        
        passed &= check_extension<D,MultiIndex>(base_enum, ext_enum, true, verbose);
    }
    
    std::cout << "}" << std::endl;
    
    return passed;
}

int main()
{
    bool passed = true;
    
    // check 1D-case
    {
        const dim_t D = 1;
        
        typedef HyperCubicShape<D> S;
        S shape(5); 
        passed &= test<D,S>(shape, false);
    }
    
    // check multidimensional case
    {
        const dim_t D = 4;
        
        {
            typedef HyperbolicCutShape<D> S;
            S shape(16);
            passed &= test<D,S>(shape, false);
        }
        
        {
            typedef HyperCubicShape<D> S;
            S shape(5);
            passed &= test<D,S>(shape, false);
        }
        
        {
            typedef LimitedHyperbolicCutShape<D> S;
            S shape(15, 8);
            passed &= test<D,S>(shape, false);
        }
        
        
    }
    
    std::cout << "Result: " << (passed ? "TEST PASSED" : "TEST FAILED") << std::endl;
    
    return 0;
}