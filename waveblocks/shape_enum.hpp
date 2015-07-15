#ifndef WAVEBLOCKS_SHAPE_ENUM_HPP
#define WAVEBLOCKS_SHAPE_ENUM_HPP

#include <vector>
#include <stdexcept>

#include "basic_types.hpp"

namespace waveblocks {

/**
 * \ingroup ShapeEnum
 * 
 * The \f$ s \f$-th slice of a shape enumeration contains all multi-indices 
 * \f$ \boldsymbol{k} \in \mathfrak{K} \f$ 
 * that satisfy \f$ \displaystyle\sum_{d=1}^{D} k_d = s \f$.
 * 
 */
template<dim_t D, class MultiIndex>
class ShapeSlice
{
private:
    std::size_t offset_;
    
    std::vector< MultiIndex > table_;
    
//     std::unordered_map< MultiIndex, std::size_t > dict_;
    
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
    typedef typename std::vector<MultiIndex>::const_iterator const_iterator;
    
    ShapeSlice() = default;
    
    ShapeSlice(const ShapeSlice& that) = default;
    
    ShapeSlice(ShapeSlice&& that)
        : offset_(that.offset_)
        , table_(std::move(that.table_))
    { }
    
    ShapeSlice &operator=(const ShapeSlice & that) = default;
    
    ShapeSlice &operator=(ShapeSlice&& that)
    {
        offset_ = that.offset_;
        table_ = std::move(that.table_);
        return *this;
    }
    
    /**
     * \brief constructs an empty slice
     */
    ShapeSlice(std::size_t offset)
        : offset_(offset)
        , table_()
    {}
    
    ShapeSlice(std::vector<MultiIndex>&& table, std::size_t offset)
        : offset_(offset)
        , table_(std::move(table))
    { }
    
    std::vector< MultiIndex > & _table()
    {
        return table_;
    }
    
    std::vector< MultiIndex > const& _table() const
    {
        return table_;
    }
    
    /**
     * \return number of nodes in all previous slices
     */
    std::size_t offset() const
    {
        return offset_;
    }
    
    /**
     * \return number of nodes in this slice
     */
    std::size_t size() const
    {
        return table_.size();
    }
    
    
    const_iterator begin() const
    {
        return table_.begin();
    }
    
    const_iterator end() const
    {
        return table_.end();
    }
    
    const_iterator cbegin() const
    {
        return table_.cbegin();
    }
    
    const_iterator cend() const
    {
        return table_.cend();
    }
    
    /**
     * \brief Returns the multi-index of the node at position <i>ordinal</i>.
     * 
     * Notice that the first node in the slice has ordinal 0 (not 1 or offset()).
     * 
     * Portable programs should never call this function with an argument that is <i>out-of-range</i>,
     * since this causes \e undefined \e behaviour.
     * 
     * <b>complexity: </b>logarithmic in the number of slice-nodes
     * 
     * \param[in] ordinal position of a node in this slice
     * \return multi-index of the specified node
     */
    MultiIndex operator[](std::size_t ordinal) const
    {
        assert(ordinal < size());
        
        return table_[ordinal];
    }
    
    bool try_find(const MultiIndex& index, std::size_t& ordinal) const
    {
        std::less< MultiIndex > comp;
        
        auto it = std::lower_bound(table_.begin(), table_.end(), index, comp);
        
        if (it != table_.end() && *it == index) {
            ordinal = it - table_.begin();
            return true;
        }
        else {
            return false;
        }
    }
    
    /**
     * \brief Returns the position of the node with multi-index \p index.
     * 
     * Notice that the first node in the slice has position 0 (not 1 or offset()).
     * 
     * Portable programs should never call this function with an node that is not part 
     * of this slice since this causes \e undefined \e behaviour.
     * 
     * Use ShapeEnumeration<D>::contains(index) to check whether this slice contains the given node.
     * 
     * <b>complexity: </b>logarithmic in the number of slice-nodes
     * 
     * \param[in] index multi-index of a node in this slice
     * \return position of the specified node
     */
    std::size_t find(const MultiIndex& index) const
    {
        std::size_t ordinal;
        
        if (!try_find(index, ordinal))
            throw std::invalid_argument("slice does not contain multi-index");
        
        return ordinal;
    }
    
    std::array<std::size_t,D> find_backward_neighbours(const MultiIndex& _index) const
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
    
    bool operator==(const ShapeSlice& that) const
    {
        return table_ == that.table_;
    }
    
    bool operator!=(const ShapeSlice& that) const
    {
        return table_ != that.table_;
    }
};

/**
 * \ingroup ShapeEnum
 * 
 * \brief A shape enumeration is a complete, ordered listing of all 
 * nodes that are part of the shape.
 * 
 * Since many algorithms operating on wavepackets use recursive formulas, 
 * shape enumerations are divided into a set of \e slices to simplify data structures.
 * 
 * The \f$ s \f$-th slice contains all multi-indices \f$ \boldsymbol{k} \in \mathfrak{K} \f$ 
 * that satisfy \f$ \displaystyle\sum_{d=1}^{D} k_d = s \f$.
 * 
 * To determine, to which slice a multi-index belongs, use:
 * \code{.cpp}
 * #include <numeric>
 * int islice = std::accumulate(index.begin(), index.end(), int(0));
 * \endcode
 */
template<dim_t D, class MultiIndex>
class ShapeEnum
{
private:
    ShapeSlice<D,MultiIndex> lower_;
    ShapeSlice<D,MultiIndex> upper_;
    
    std::vector< ShapeSlice<D, MultiIndex> > slices_;
    std::size_t n_entries_;
    MultiIndex limits_;
    
public:
    ShapeEnum() = default;
    
    ShapeEnum(std::vector< ShapeSlice<D, MultiIndex> >&& slices,
              std::size_t n_entries,
              MultiIndex limits)
        : lower_(0)
        , upper_(n_entries)
        , slices_(std::move(slices))
        , n_entries_(n_entries)
        , limits_(limits)
    { }
    
    ShapeEnum(const ShapeEnum& that) = default;
    
    ShapeEnum(ShapeEnum&& that)
        : lower_(std::move(that.lower_))
        , upper_(std::move(that.upper_))
        , slices_(std::move(that.slices_))
        , n_entries_(that.n_entries_)
        , limits_(that.limits_)
    { }
    
    ShapeEnum &operator=(const ShapeEnum& that) = default;
    
    ShapeEnum &operator=(ShapeEnum&& that)
    {
        n_entries_ = that.n_entries_;
        limits_ = that.limits_;
        lower_ = std::move(that.lower_);
        upper_ = std::move(that.upper_);
        slices_ = std::move(that.slices_);
        return *this;
    }
    
    /**
     * \brief returns a reference to a slice
     * 
     * This function does not fail if an 'invalid' slice index is passed. 
     * If slice index is negative then this function returns an empty slice with offset 0.
     * If slice index is equals or larger then the number of (non-empty) slices the this function
     * returns an empty slice with offset set equals to total number of entries in all slices.
     * 
     * \param[in] islice ordinal of the desired slice
     * \return reference to slice
     */
    const ShapeSlice<D, MultiIndex>& slice(int islice) const
    {
        if (islice < 0)
            return lower_;
        else if (islice >= (int)slices_.size())
            return upper_;
        else
            return slices_[islice];
    }
    
    ShapeSlice<D, MultiIndex>& slice(int islice)
    {
        if (islice < 0)
            return lower_;
        else if (islice >= (int)slices_.size())
            return upper_;
        else
            return slices_[islice];
    }
    
    const std::vector< ShapeSlice<D, MultiIndex> >& slices() const
    {
        return slices_;
    }
    
    std::size_t n_entries() const
    {
        return n_entries_;
    }
    
    /**
     * \return number of non-empty (probably) slices
     */
    int n_slices() const
    {
        return (int)slices_.size();
    }
    
    const MultiIndex& limits() const
    {
        return limits_;
    }
    
    int limit(dim_t axis) const
    {
        return limits_[axis];
    }
    
    bool operator==(const ShapeEnum &that) const
    {
        if (n_entries_ != that.n_entries_)
            return false;
        
        for (int islice = n_entries_-1; islice >= 0; islice--) {
            if (slice(islice) != that.slice(islice))
                return false;
        }
        
        return true;
    }
    
    bool operator!=(const ShapeEnum &that) const
    {
        return !operator==(that);
    }
};

}

#endif