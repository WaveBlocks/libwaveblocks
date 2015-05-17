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

template<dim_t D>
class ShapeEnumeration
{
public:
    class Slice
    {
    public:
        virtual ~Slice() { }
        
        virtual std::size_t offset() const = 0;
        
        virtual std::size_t size() const = 0;
        
        virtual std::array<int,D> operator[](std::size_t ordinal) const = 0;
        
        virtual std::size_t find(const std::array<int,D> &index) const = 0;
        
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
    
    virtual std::size_t size() const = 0;
    
    virtual std::array<int,D> operator[](std::size_t ordinal) const = 0;
    
    virtual std::size_t find(const std::array<int,D> &index) const = 0;
    
    virtual const Slice &slice(std::size_t islice) const = 0;
    
    virtual std::size_t count_slices() const = 0;
    
    virtual bool contains(const std::array<int,D> &index) const = 0;
    
    virtual std::string description() const = 0;
    
    /**
     * \return smallest tuple L s.t for every k element of K: k_d <= L_d
     */
    virtual int bbox(dim_t axis) const = 0;
    
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
    
    Iterator begin() const
    {
        return Iterator(this, 0);
    }
    
    Iterator end() const
    {
        return Iterator(this, size());
    }
};

template<dim_t D, class MultiIndex, class S>
class SlicedShapeEnumeration : public ShapeEnumeration<D>
{
public:
    class Slice;
    
private:
    S shape_;
    
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
    
    SlicedShapeEnumeration(S shape);
    
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
SlicedShapeEnumeration<D,MultiIndex,S>::SlicedShapeEnumeration(S shape)
    : shape_(shape)
    , use_dict_(false) //using a hashmap is slower than binary search
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
    std::size_t ordinal = 0;
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