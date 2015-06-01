#ifndef WAVEBLOCKS_SHAPE_ENUMERATION_BASE
#define WAVEBLOCKS_SHAPE_ENUMERATION_BASE

#include <array>
#include <string>

#include "basic_types.hpp"

namespace waveblocks {

/**
 * \brief Base class for all shape enumeration implementations.
 * 
 * The purpose of this class is to unify access to all implementations by hiding the
 * real implementation details behind virtual functions.
 * 
 * Since translating multi-indices to ordinals is quite tricky if not impossible,
 * implementations will use a dictionary or lookup-table to perform this conversions.
 * If a shape contains several million nodes, such data structures consume quite a lot of memory, 
 * so implementations will somehow compress multi-indices to save memory. All multi-indices 
 * are passed as as a tuple of integers (std::array<int,D>). Implementations will transform it
 * into the internal format.
 * 
 * Instantiation
 * -------------
 * \code{.cpp}
 * const dim_t D = 4;
 * TinyMultiIndex<std::size_t,D> MultiIndex;
 * typedef HyperbolicCutShape<D> S;
 * 
 * S shape(7.0);
 * 
 * ShapeEnumeration<D> *enumeration = new SlicedShapeEnumeration< D, MultiIndex, S >(shape);
 * \endcode
 * 
 * \tparam D number of multi-index dimensions
 * 
 * \see TinyMultiIndex A compressed multi-index type that represents all multi-indices using a single integer.
 */
template<dim_t D>
class ShapeEnumeration
{
public:
    /**
     * \brief A slice contains all 
     * 
     * The <b>s</b>-th slice contains all nodes <b>k</b> with the property: \f$ \sum_{i=1}^{D} k_i = s \f$
     */
    class Slice
    {
    public:
        virtual ~Slice() { }
        
        /**
         * \return number of nodes in all previous slices
         */
        virtual std::size_t offset() const = 0;
        
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
         * Portable programs should never call this function with an argument 
         * that is inexistant in this slice since this causes <i>undefined behaviour</i>.
         * 
         * Use contains(index) to check whether this slice contains this node.
         * 
         * <b>complexity: </b>logarithmic in the number of slice-nodes
         * 
         * \param[in] index multi-index of a node in this slice
         * \return position of the specified node
         */
        virtual std::size_t find(const std::array<int,D> &index) const = 0;
        
        /**
         * \brief const_iterator over a slice to support foreach statements
         */
        class Iterator
        {
        private:
            const Slice *ref_;
            std::size_t ientry_;
            
        public:
            Iterator(const Slice *ref, std::size_t ientry)
            : ref_(ref)
            , ientry_(ientry)
            { }
            
            Iterator(const Iterator &other)
                : ref_(other.ref_)
                , ientry_(other.ientry_)
            { }
            
            Iterator &operator=(const Iterator &other)
            {
                ref_ = other.ref_;
                ientry_ = other.ientry_;
                return *this;
            }
            
            Iterator &operator++()
            {
                ++ientry_;
                return *this;
            }
            
            std::array<int,D> operator*() const
            {
                return (*ref_)[ientry_];
            }
            
            bool operator==(const Iterator &other) const
            {
                return ientry_ == other.ientry_ && ref_ == other.ref_;
            }
            
            bool operator!=(const Iterator &other) const
            {
                return !operator==(other);
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
    
    virtual ~ShapeEnumeration() { }
    
    /**
     * \return number of nodes
     */
    virtual std::size_t size() const = 0;
    
    virtual std::array<int,D> operator[](std::size_t ordinal) const = 0;
    
    virtual std::size_t find(const std::array<int,D> &index) const = 0;
    
    /**
     * \param[in] islice slice-index (first slice has index 0)
     * \return reference to a slice
     */
    virtual const Slice &slice(std::size_t islice) const = 0;
    
    virtual std::size_t count_slices() const = 0;
    
    /**
     * \brief Checks whether this slice contains a node with multi-index <i>index</i>.
     * 
     * To proof that a shape contains a specific node, it is sufficient to consult 
     * the shape's member function limit(axis, index). However querying the ordinal of a 
     * node however is expensive since it usually requires a dictionary lookup.
     * Therefore contains(<i>index</i>) is very cheap compared to find(<i>index</i>).
     * 
     * <b>complexity: </b> constant-time
     * 
     * \param index multi-index of the node
     * \return true if the shape contains the node; false otherwise
     */
    virtual bool contains(const std::array<int,D> &index) const = 0;
    
    virtual std::string description() const = 0;
    
    /**
     * \return smallest tuple L s.t for every k element of K: k_d <= L_d
     */
    virtual int bbox(dim_t axis) const = 0;
    
    /**
     * \brief const_iterator over the whole shape to support foreach statements
     */
    class Iterator
    {
    private:
        const ShapeEnumeration *ref_;
        std::size_t ientry_;
        
    public:
        Iterator(const ShapeEnumeration *ref, std::size_t ientry)
            : ref_(ref)
            , ientry_(ientry)
        { }
        
        Iterator(const Iterator &other)
            : ref_(other.ref_)
            , ientry_(other.ientry_)
        { }
        
        Iterator &operator=(const Iterator &other)
        {
            ref_ = other.ref_;
            ientry_ = other.ientry_;
            return *this;
        }
        
        Iterator &operator++()
        {
            ++ientry_;
            return *this;
        }
        
        std::array<int,D> operator*() const
        {
            return (*ref_)[ientry_];
        }
        
        bool operator==(const Iterator &other) const
        {
            return ientry_ == other.ientry_ && ref_ == other.ref_;
        }
        
        bool operator!=(const Iterator &other) const
        {
            return !operator==(other);
        }
    };
    
    /**
     * \brief 
     */
    Iterator begin() const
    {
        return Iterator(this, 0);
    }
    
    /**
     * \brief
     */
    Iterator end() const
    {
        return Iterator(this, size());
    }
};

}

#endif