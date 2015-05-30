#include <array>
#include <string>
#include <iostream>

#include <waveblocks/shape_hypercubic.hpp>
#include <waveblocks/shape_hyperbolic.hpp>
#include <waveblocks/shape_extended.hpp>
#include <waveblocks/lexical_shape_enumerator.hpp>
#include <waveblocks/stdarray2stream.hpp>

using namespace waveblocks;

/**
 * \tparam D dimensions
 * \tparam S type of test shape
 * \param[in] base test shape
 * \param[in] strict if enabled: print all unneeded nodes of extended shape
 */
template<dim_t D, class S>
void check_extended_shape(const S &base, bool strict)
{
    ExtendedShape<D,S> extended(base);
    
    {
        std::cout << "check extended shape implementation (" << extended.description() << ") of " << base.description() << " {" << std::endl;
        
        LexicalIndexGenerator<D,std::array<int,D>,S> base_gen(base);
        
        //ensure that for every node inside base shape, the extended shape contains all forward neighbours
        do {
            const std::array<int,D> node = base_gen.index();
            for (dim_t d = 0; d < D; d++) {
                std::array<int,D> fnode = node; fnode[d] += 1;
                
                if (fnode[d] > extended.limit(fnode,d)) {
                    std::cout << "   [FAILURE] extended shape does not contain " << fnode << std::endl;
                }
            }
        } while (base_gen.forward());
    }
    
    //count unneeded nodes contained in extended shape
    {
        std::size_t nUnneeded = 0, total = 0;
        
        LexicalIndexGenerator<D,std::array<int,D>,ExtendedShape<D,S>> ext_gen(extended);
        
        do {
            const std::array<int,D> node = ext_gen.index();
            int count = 0;
            bool enable = false;
            for (dim_t d = 0; d < D; d++) {
                if (node[d] != 0) {
                    enable = true;
                    std::array<int,D> bnode = node; bnode[d] -= 1;
                    
                    if (bnode[d] <= base.limit(bnode,d)) {
                        ++count;
                    }
                }
            }
            
            ++total;
            if (enable && count == 0) {
                ++nUnneeded;
                if (strict)
                    std::cout << "   [WARNING] extended shape contains unneeded entry: "<<node<<std::endl;
            }
        } while(ext_gen.forward());
        
        std::cout << "   [INFO] overhead : "<<nUnneeded<<"/"<<total<<" = "<<
                double(nUnneeded)/double(total)*100.<<"%"<<std::endl;
    }
    
    std::cout << "}" << std::endl;
}

int main()
{
    const dim_t D = 4;
    
    {
        typedef HyperCubicShape<D> S;
        S shape(5);
        check_extended_shape<D,S>(shape, false);
    }
    
    {
        typedef LimitedHyperbolicCutShape<D> S;
        S shape(15.0, 8);
        check_extended_shape<D,S>(shape, false);
    }
    
    {
        typedef HyperbolicCutShape<D> S;
        S shape(10.0);
        check_extended_shape<D,S>(shape, false);
    }
    
    return 0;
}