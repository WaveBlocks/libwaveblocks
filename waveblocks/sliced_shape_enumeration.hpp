#ifndef SLICED_SHAPE_ENUMERATION
#define SLICED_SHAPE_ENUMERATION

#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <array>
#include <string>
#include <sstream>

#include "basic_types.hpp"
#include "lexical_shape_enumerator.hpp"

namespace waveblocks {

/**
 * \brief Base class for all shape enumeration implementations.
 * 
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
template<dim_t D, class MultiIndex, class S>
class SlicedShapeEnumeration : public ShapeEnumeration<D>
{
public:
    class Slice;
    
private:
    /**
     * shape description
     */
    S shape_;
    
    /**
     * flag whether a dictionary is used.
     */
    bool use_dict_;
    
    std::vector< MultiIndex > table_;
    
    std::vector< Slice > slices_;
    
public:
    class Slice : public ShapeEnumeration<D>::Slice
    {
        friend class SlicedShapeEnumeration;
        
    private:
        const SlicedShapeEnumeration *enumeration_;
        
        std::size_t start_;
        std::size_t end_;
        
        std::unordered_map< MultiIndex, std::size_t > dict_;
        
        typedef typename std::vector< MultiIndex >::const_iterator Iterator;
        
        Iterator begin() const
        {
            return enumeration_->table_.begin() + start_;
        }
        
        Iterator end() const
        {
            return enumeration_->table_.begin() + end_;
        }
        
    public:
        std::size_t offset() const override
        {
            return start_;
        }
        
        std::size_t size() const override
        {
            return end_ - start_;
        }
        
        std::array<int,D> operator[](std::size_t ientry) const override
        {
            assert (ientry < size());
            
            return static_cast< std::array<int,D> >( *(begin() + ientry) );
        }
        
        std::size_t find(const std::array<int,D> &_index) const override
        {
            MultiIndex index(_index);
            
            if (enumeration_->use_dict_) {
                
                auto it = dict_.find(index);
                if (it == dict_.end())
                    return size(); //shape does not contain node
                else
                    return it->second;
                
            } else {
                std::less< MultiIndex > comp;
                
                auto it = std::lower_bound(begin(), end(), index, comp);
                
                if (*it == index)
                    return it - begin();
                else
                    return size();
            }
        }
    };
    
    /**
     * \brief Takes a shape description and builds a lookup-table for it.
     * 
     * \param[in] shape shape description
     * \param[in] use_dict option whether to use a dictionary (see class description for details)
     */
    SlicedShapeEnumeration(S shape, bool use_dict = false);
    
    std::size_t size() const override
    {
        return table_.size();
    }
    
    std::array<int,D> operator[](std::size_t ordinal) const override
    {
        assert(ordinal < size());
        
        return static_cast< std::array<int,D> >( table_[ordinal] );
    }
    
    const Slice &slice(std::size_t islice) const override
    {
        return slices_[islice];
    }
    
    std::size_t count_slices() const override
    {
        return slices_.size();
    }
    
    std::size_t find(const std::array<int,D> &index) const override
    {
        std::size_t islice = 0;
        for (dim_t i = 0; i < D; i++)
            islice += index[i];
        
        if (islice >= slices_.size())
            return size(); //entry does not exist
        
        std::size_t ientry =  slice(islice).find(index);
        
        if (ientry >= slice(islice).size())
            return size(); //entry does not exist
        else
            return slice(islice).offset() + ientry;
    }
    
    bool contains(const std::array<int,D> &index) const override
    {
        return index[0] <= shape_.template limit< std::array<int,D> >(index, 0);
    }
    
    std::string description() const override
    {
        std::stringstream out;
        out << "ShapeEnumeration {";
        out << "dimension: " << D << ", ";
        out << "shape: " << shape_.description() << ", ";
        out << "#entries: " << size() << ", ";
        out << "#slices: " << count_slices() << "}";
        return out.str();
    }
    
    int bbox(dim_t axis) const override
    {
        return shape_.bbox(axis);
    }
};

template<dim_t D, class MultiIndex, class S>
SlicedShapeEnumeration<D,MultiIndex,S>::SlicedShapeEnumeration(S shape, bool use_dict)
    : shape_(shape)
    , use_dict_(use_dict) //using a dictionary is slower than binary search
    , table_()
{
    //check multi-index type
    for (dim_t d = 0; d < D; d++) {
        if (shape.bbox(d) > MultiIndex::limit(d))
            throw std::runtime_error("multi-index type is not suitable. reason: overflow");
    }
    
    LexicalIndexGenerator<D,MultiIndex,S> it(shape);
    
    //holds for each slice its own vector
    std::vector< std::vector< MultiIndex >* > temp;
    
    //generate nodes and put them into the corresponding slice
    do {
        MultiIndex index = it.index();
        
        std::size_t islice = 0;
        for (dim_t i = 0; i < D; i++)
            islice += index[i];
        
        assert(islice <= temp.size()); //assert that lexical index generator does not jump
        
        if (islice == temp.size()) {
            temp.push_back( new std::vector< MultiIndex >() );
        }
        
        temp[islice]->push_back(index);
    } while (it.forward());
    
    //all slices are now populated
    //copy all nodes into the same vector and store offsets
    slices_ = std::vector<Slice>(temp.size());
    for (std::size_t i = 0; i < temp.size(); i++) {
        slices_[i].enumeration_ = this;
        
        slices_[i].start_ = table_.size();
        table_.insert(table_.end(), temp[i]->begin(), temp[i]->end());
        slices_[i].end_ = table_.size();
        
        //create dictionary
        if (use_dict_) {
            std::size_t ordinal = 0;
            for (auto index : slices_[i]) {
                slices_[i].dict_[index] = ordinal++;
            }
        }
        
        delete temp[i];
        temp[i] = nullptr;
    }
}

}

#endif