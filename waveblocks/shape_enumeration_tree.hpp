#ifndef WAVEBLOCKS_SHAPE_ENUMERATION_TREE
#define WAVEBLOCKS_SHAPE_ENUMERATION_TREE

#include <array>
#include <vector>
#include <stdexcept>

#include "basic_types.hpp"

namespace waveblocks {

template<dim_t D, class S>
class SparseManhattanTree
{
private:
    typedef unsigned size_type;
    
    struct node_
    {
        /**
         * \brief size of all subtrees
         */
        size_type size;
        
        /**
         * \brief number of elements in all subtrees
         */
        size_type entries;
    };
    
    int sum_;
    std::vector<node_> tree_;
    
    /**
     * \brief stores for each subtree type its capacity
     * 
     * formula (recursive):
     *  cap(0,k) = 1
     *  cap(d,0) = 1
     *  cap(d,k) = cap(d-1,k) + cap(d,k-1)
     * 
     * Access: cap(d,k) = tree_cap[k*D + d]
     */
    std::vector<size_type> tree_cap_;
    
    inline size_type &capacity(int d, int k)
    {
        return tree_cap_[k*D + d];
    }
    
public:
    /**
     * \param sum 
     */
    SparseManhattanTree(int sum)
        : tree_cap_(D*(sum+1))
    {
        // compute lookup table of subtree capacities
        for (int k = 0; k <= sum; k++) {
            capacity(0,k) = 1;
        }
        
        for (int d = 0; d < D; d++) {
            capacity(d,0) = 1;
        }
        
        for (int k = 1; k <= sum; k++) {
            for (int d = 1; d < D; d++) {
                capacity(d,k) = capacity(d-1,k) + capacity(d,k-1);
            }
        }
        
        // fill tree
        
    }
    
    size_type size() const
    {
        return tree_[0].entries;
    }
    
    std::array<int,D> operator[](size_type ordinal) const
    {
        if (ordinal >= tree_[0].entries())
            throw std::out_of_range("ordinal is out of range");
        
        auto it = tree_.cbegin() + 1;
        
        std::array<int,D> index;
        
        int kmax = sum_;
        for (int depth = 0; depth < D-1; depth++) {
            for (int k = kmax; k >= 0; k--) {
                if (ordinal < it->entries) {
                    // enter subtree
                    kmax = k;
                    index[depth] = kmax - k;
                    ++it;
                    break;
                } else {
                    // skip subtree
                    ordinal -= it->entries;
                    it += it->size;
                }
            }
        }
        
        index[D-1] = kmax;
        
        return index;
    }
    
    size_type find(const std::array<int,D> &index) const
    {
        // check input
        {
            int sum = 0;
            for (int i = 0; i < D; i++)
                sum += index[i];
            
            if (sum != sum_)
                throw std::domain_error("slice does not contain multi-index");
        }
        
        auto it = tree_.cbegin() + 1;
        
        int kmax = sum_;
        
        size_type ordinal = 0;
        
        for (int depth = 0; depth < D-1; depth++) {
            // go to right subtree
            for (int i = 0; i < index[depth]; i++) {
                ordinal += it->entries;
                it += it->size;
            }
            
            kmax -= index[depth];
            
            ++it; // descend to first child
            
            // requested subtree is empty => slice does not contain multi-index
            // check multi-index with contains(index)!
            if (it->entries == 0)
                throw std::domain_error("slice does not contain multi-index");
        }
    }
};

}

#endif