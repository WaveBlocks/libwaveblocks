#ifndef WAVEBLOCKS_SHAPE_ENUMERATION_DEFAULT
#define WAVEBLOCKS_SHAPE_ENUMERATION_DEFAULT

#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <array>
#include <string>
#include <sstream>

#include "basic_types.hpp"
#include "lexical_shape_enumerator.hpp"
#include "shape_enumeration_base.hpp"

namespace waveblocks {

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
class DefaultShapeEnumeration : public ShapeEnumeration<D>
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
        friend class DefaultShapeEnumeration;
        
    private:
        const DefaultShapeEnumeration *enumeration_;
        
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
    DefaultShapeEnumeration(S shape, bool use_dict = false);
    
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
DefaultShapeEnumeration<D,MultiIndex,S>::DefaultShapeEnumeration(S shape, bool use_dict)
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