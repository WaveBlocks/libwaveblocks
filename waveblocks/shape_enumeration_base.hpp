#ifndef WAVEBLOCKS_SHAPE_ENUMERATION_HPP
#define WAVEBLOCKS_SHAPE_ENUMERATION_HPP

#include <iterator>
#include <array>
#include <vector>

#include "basic_types.hpp"

namespace waveblocks {

template<dim_t D>
class ShapeSlice
{
public:
    virtual ~ShapeSlice() { }
    
    /**
     * \return number of nodes in all previous slices
     */
    virtual std::size_t offset() const = 0;
    
    /**
     * \return index/ordinal of this slice
     */
    virtual std::size_t slice_index() const = 0;
    
    /**
     * \return number of nodes in this slice
     */
    virtual std::size_t size() const = 0;
    
    /**
     * \brief Returns the multi-index of the node at position <i>ordinal</i>.
     * 
     * Notice that the first node in the slice has ordinal 0 (not 1 or offset()).
     * 
     * Portable programs should never call this function with an argument that is <i>out-of-range</i>,
     * since this causes <i>undefined behaviour</i>.
     * 
     * <b>complexity: </b>logarithmic in the number of slice-nodes
     * 
     * \param[in] ordinal position of a node in this slice
     * \return multi-index of the specified node
     */
    virtual std::array<int,D> operator[](std::size_t ordinal) const = 0;
    
    /**
     * \brief Returns the position of the node with multi-index <i>index</i>.
     * 
     * Notice that the first node in the slice has position 0 (not 1 or offset()).
     * 
     * Portable programs should never call this function with an node that is not part 
     * of this slice since this causes <i>undefined behaviour</i>.
     * 
     * Use ShapeEnumeration<D>::contains(index) to check whether this slice contains the given node.
     * 
     * <b>complexity: </b>logarithmic in the number of slice-nodes
     * 
     * \param[in] index multi-index of a node in this slice
     * \return position of the specified node
     */
    virtual std::size_t find(const std::array<int,D> &index) const = 0;
    
    /**
     * \brief Retrieves ordinals of all backward neighbours of a given node very efficiently.
     * 
     * Notice that this method assumes that the given node <b>is part of the shape</b>.
     * Therefore this method does not need to perform any contains()-checks since it 
     * knows that all backward neighbours exist, except when the given node contains
     * some zero entries. In the latter case, this method returns an undefined
     * ordinal.
     * 
     * If the given node is not part of the shape, then the behaviour is undefined.
     * 
     * \param[in] index node that <b>is part of the shape</b>
     * \return For each backward neighbour its ordinal. An ordinal is undefined if its node does not exist.
     */
    virtual std::array<std::size_t,D> find_backward_neighbours(const std::array<int,D> &index) const = 0;
    
    /**
     * \brief const_iterator over a slice to support foreach statements
     */
    class Iterator : public std::iterator<std::random_access_iterator_tag, std::array<int,D>>
    {
    private:
        const ShapeSlice *slice_;
        std::size_t ientry_;
        
    public:
        Iterator(const ShapeSlice *slice, std::size_t ientry)
        : slice_(slice)
        , ientry_(ientry)
        { }
        
        // --- Iterator ---
        
        Iterator() = default;
        
        Iterator(const Iterator &other)
            : slice_(other.slice_)
            , ientry_(other.ientry_)
        { }
        
        Iterator &operator=(const Iterator &other) const
        {
            slice_ = other.slice_;
            ientry_ = other.ientry_;
            return *this;
        }
        
        ~Iterator() = default;
        
        // --- Forward Iterator ---
        
        bool operator!=(const Iterator &other) const
        {
            return ientry_ != other.ientry_;
        }
        
        bool operator==(const Iterator &other) const
        {
            return ientry_ == other.ientry_;
        }
        
        std::array<int,D> operator*() const
        {
            return (*slice_)[ientry_];
        }
        
        Iterator &operator++()
        {
            ++ientry_;
            return *this;
        }
        
        Iterator &operator++(int)
        {
            ++ientry_;
            return *this;
        }
        
        // --- Bidirectional Iterator ---
        
        Iterator &operator--()
        {
            --ientry_;
            return *this;
        }
        
        Iterator &operator--(int)
        {
            --ientry_;
            return *this;
        }
        
        
        // --- Random Access Iterator ---
        
        Iterator &operator+=(std::ptrdiff_t n)
        {
            ientry_ += n;
            return *this;
        }
        
        Iterator &operator-=(std::ptrdiff_t n)
        {
            ientry_ -= n;
            return *this;
        }
        
        Iterator operator+(std::ptrdiff_t n) const
        {
            return Iterator(slice_, ientry_+n);
        }
        
        Iterator operator-(std::ptrdiff_t n) const
        {
            return Iterator(slice_, ientry_-n);
        }
        
        std::array<int,D> operator[](std::ptrdiff_t n) const
        {
            return (*slice_)[ientry_+n];
        }
        
        bool operator<(const Iterator &that) const
        {
            return ientry_ < that.ientry_;
        }
        
        bool operator>(const Iterator &that) const
        {
            return ientry_ > that.ientry_;
        }
        
        bool operator<=(const Iterator &that) const
        {
            return ientry_ <= that.ientry_;
        }
        
        bool operator>=(const Iterator &that) const
        {
            return ientry_ >= that.ientry_;
        }
    };
    
    Iterator begin() const
    {
        return Iterator(this, 0);
    }
    
    Iterator end() const
    {
        return Iterator(this, size());
    }
};

template<dim_t D>
class ShapeEnumeration
{
public:
    /**
     * \return number of (ideally non-empty) slices
     */
    virtual std::size_t n_slices() const = 0;
    
    /**
     * \param[in] islice index of requested slice
     * \return requested slice or empty slice if (islice >= #slices)
     */
    virtual const ShapeSlice<D>& slice(std::size_t islice) const = 0;
    
    /**
     * \return forwards
     */
    virtual bool contains(const std::array<int,D> &index) const = 0;
    
    /**
     * \brief retrieves the index of the outmost node in a given axis
     * 
     * \return 
     */
    virtual int bbox(dim_t axis) const = 0;
    
    /**
     * \brief returns number of nodes that are part of the shape
     * 
     * \return number of nodes
     */
    virtual std::size_t size() const
    {
        std::size_t sum = 0;
        for (std::size_t islice = 0; islice < n_slices(); islice++) {
            sum += slice(islice).size();
        }
        return sum;
    }
    
    /**
     * \brief range that contains all slices
     */
    struct Slices
    {
    private:
        const ShapeEnumeration *ref_;
        
    public:
        Slices(const ShapeEnumeration *ref)
            : ref_(ref)
        { }
        
        /**
         * 
         * 
         * \param[in] islice index of requested slice
         * \return reference to requested slice
         */
        const ShapeSlice<D>& operator[](std::size_t islice) const
        {
            return ref_->slice(islice);
        }
        
        /**
         * \return number of slices
         */
        std::size_t size() const
        {
            return ref_->n_slices();
        }
        
        struct Iterator : std::iterator<std::forward_iterator_tag, ShapeSlice<D> >
        {
        private:
            const ShapeEnumeration *ref_;
            std::size_t islice_;
            
        public:
            Iterator(const ShapeEnumeration *ref, std::size_t islice)
                : ref_(ref)
                , islice_(islice)
            { }
            
            // --- Iterator ---
        
            Iterator() = default;
            
            Iterator(const Iterator &other)
                : ref_(other.ref_)
                , islice_(other.islice_)
            { }
            
            Iterator &operator=(const Iterator &other) const
            {
                ref_ = other.ref_;
                islice_ = other.islice_;
                return *this;
            }
            
            ~Iterator() = default;
            
            // --- Forward Iterator ---
            
            bool operator!=(const Iterator &other) const
            {
                return islice_ != other.islice_;
            }
            
            bool operator==(const Iterator &other) const
            {
                return islice_ == other.islice_;
            }
            
            const ShapeSlice<D>& operator*() const
            {
                return ref_->slice(islice_);
            }
            
            const ShapeSlice<D>* operator->() const
            {
                return & ref_->slice(islice_);
            }
            
            Iterator &operator++()
            {
                ++islice_;
                return *this;
            }
            
            Iterator &operator++(int)
            {
                ++islice_;
                return *this;
            }
        };
        
        Iterator begin() const
        {
            return Iterator{ref_, 0};
        }
        
        Iterator end() const
        {
            return Iterator{ref_, ref_->n_slices()};
        }
    };
    
    Slices slices() const
    {
        return Slices{this};
    }
};

}

#endif