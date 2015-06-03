#ifndef WAVEBLOCKS_DEFAULT_SHAPE_ENUMERATION_HPP
#define WAVEBLOCKS_DEFAULT_SHAPE_ENUMERATION_HPP

#include "shape_enumeration_base.hpp"

#include <stdexcept>
#include <algorithm>

namespace waveblocks {

template<dim_t D, class MultiIndex, class S>
class DefaultShapeEnumeration;

/**
 * \brief Default implementation of a shape enumeration. 
 * 
 * This class takes a shape description object and builds a lookup-table to
 * perform queries. It lets the user freely choose an appropriate type to
 * represent multi-indices internally.
 * 
 * <b>Implementation details</b>
 * 
 * This class uses a vector to do <i>ordinal -> multi-index</i> queries.
 * For <i>multi-index -> ordinal</i> it does binary search by default. 
 * The constructor provides an option to use an hashmap instead. 
 * However it turned out that binary search is slightly faster than 
 * a hashmap regardless of the hash-function.
 * 
 * \tparam D number of multi-index dimensions
 * \tparam MultiIndex 
 * \parblock
 * Type to internally represent a multi-index. <br>
 * A custom type should provide same interface as std::array<int,D>. <br>
 * Furthermore a custom type must specialize std::less, std::hash, std::equal_to.
 * \endparblock
 * \tparam S shape description class
 * 
 * \see TinyMultiIndex A compressed multi-index type that represents all multi-indices using a single integer.
 */
template<dim_t D, class MultiIndex>
class DefaultShapeSlice : public ShapeSlice<D>
{
    template<dim_t D_, class M_, class S_>
    friend class DefaultShapeEnumeration;
private:
    std::size_t islice;
    
    std::size_t offset_;
    
    /**
     * flag whether a dictionary is used.
     */
    bool use_dict_;
    
    std::vector< MultiIndex > table_;
    
    std::unordered_map< MultiIndex, std::size_t > dict_;
    
    inline MultiIndex forward_(MultiIndex index, dim_t axis) const
    {
        index[axis] += 1;
        return index;
    }
    
    inline MultiIndex backward_(MultiIndex index, dim_t axis) const
    {
        index[axis] -= 1;
        return index;
    }
    
public:
    DefaultShapeSlice() = default;
    
    std::size_t offset() const override
    {
        return offset_;
    }
    
    std::size_t size() const override
    {
        return table_.size();
    }
    
    std::size_t slice_index() const override
    {
        return islice;
    }
    
    std::array<int,D> operator[](std::size_t ordinal) const override
    {
        assert(ordinal < size());
        
        return static_cast< std::array<int,D> >( table_[ordinal] );
    }
    
    std::size_t find(const std::array<int,D> &_index) const override
    {
        
        
        MultiIndex index(_index);
        
        if (use_dict_) {
            
            auto it = dict_.find(index);
            if (it == dict_.end())
                throw std::invalid_argument("slice does not contain multi-index");
            else
                return it->second;
            
        } else {
            std::less< MultiIndex > comp;
            
            auto it = std::lower_bound(table_.begin(), table_.end(), index, comp);
            
            if (*it == index)
                return it - table_.begin();
            else
                throw std::invalid_argument("slice does not contain multi-index");
        }
    }
    
    virtual std::array<std::size_t,D> findBackwardNeighbours(const std::array<int,D> &_index) const override
    {
        std::array<std::size_t,D> ordinals{}; //zero initialize
        
        MultiIndex index(_index);
        
        std::less< MultiIndex > comp;
        
        // find last non-zero entry
        dim_t dlast = D-1;
        while (dlast >= 0 && index[dlast] == 0) {
            --dlast;
        }
        
        if (dlast >= 0) {
            auto lower = table_.begin();
            
            auto upper = std::lower_bound(lower, table_.end(), backward_(index, dlast), comp);
            ordinals[dlast] = upper - table_.begin();
            
            for (dim_t i = 0; i < dlast; i++) {
                if (index[i] != 0) {
                    lower = std::lower_bound(lower, upper, backward_(index, i), comp);
                    ordinals[i] = lower - table_.begin();
                }
            }
        }
        
        return ordinals;
    }
};

template<dim_t D, class MultiIndex, class S>
class DefaultShapeEnumeration : public ShapeEnumeration<D>
{
private:
    std::size_t size_;
    std::size_t n_slices_;
    S shape_;
    std::vector< DefaultShapeSlice<D, MultiIndex> > slices_;
    
public:
    DefaultShapeEnumeration(const S& shape, bool use_dict_ = false)
        : shape_(shape)
        , slices_()
    {
        // check multi-index type for compatibility
        {
            for (dim_t d = 0; d < D; d++) {
                if (shape_.bbox(d) > MultiIndex::limit(d))
                    throw std::runtime_error("multi-index type is not suitable. reason: overflow");
            }
        }
        
        // initialize slice vector
        {
            std::size_t sum = 0;
            for (dim_t d = 0; d < D; d++) {
                sum += shape_.bbox(d)+1;
            }
            
            slices_.resize(sum);
            for (std::size_t i = 0; i < sum; i++) {
                slices_[i].islice = i;
                slices_[i].use_dict_ = use_dict_;
            }
        }
        
        // enumerate shape & store all multi-indices
        {
            MultiIndex index{}; //zero initialize
            std::size_t islice = 0;
            
            while (true) {
                // iterate over last axis
                for (dim_t i = 0; i <= shape_.template limit<MultiIndex>(index,D-1); i++) {
                    index[D-1] = i;
                    
                    if (use_dict_)
                        slices_[islice+i].dict_[index] = slices_[islice+i].table_.size();
                    
                    slices_[islice+i].table_.push_back(index);
                }
                index[D-1] = 0;
                
                // iterate over other axes
                if (D > 1) {
                    dim_t j = D-2;
                    while ((int)index[j] == shape_.template limit<MultiIndex>(index,j)) {
                        islice -= index[j];
                        index[j] = 0;
                        if (j == 0)
                            goto enumeration_complete;
                        else
                            j = j-1;
                    }
                    islice += 1;
                    index[j] += 1;
                }
            }
enumeration_complete:
            (void)0;
        }
        
        // determine number of slices & slice offsets
        {
            std::size_t offset = 0;
            n_slices_ = 0;
            while (n_slices_ < slices_.size() && slices_[n_slices_].table_.size() != 0) {
                slices_[n_slices_].offset_ = offset;
                
                offset += slices_[n_slices_].table_.size();
                
                ++n_slices_;
            }
            
            slices_.resize(n_slices_);
            
            size_ = offset;
        }
        
        // end of constructor
    }
    
    std::size_t n_slices() const override
    {
        return n_slices_;
    }
    
    std::size_t size() const override
    {
        return size_;
    }
    
    const ShapeSlice<D>& slice(std::size_t islice) const override
    {
        const static DefaultShapeSlice<D,MultiIndex> empty_slice{};
        
        if (islice >= slices_.size()) {
            return empty_slice; // empty slice
        } else {
            return slices_[islice];
        }
    }
    
    bool contains(const std::array<int,D> &index) const override
    {
        return index[0] <= shape_.template limit< std::array<int,D> >(index, 0);
    }
    
    /**
     * 
     */
    int bbox(dim_t axis) const override
    {
        return shape_.bbox(axis);
    }
};

}

#endif