#ifndef SLICED_SHAPE_ENUMERATION
#define SLICED_SHAPE_ENUMERATION

#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <valarray>

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
        
        virtual MultiIndex<D> operator[](std::size_t ordinal) const = 0;
        
        virtual std::size_t find(MultiIndex<D> index) const = 0;
        
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
            
            MultiIndex<D> operator*() const
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
    
    virtual MultiIndex<D> operator[](std::size_t ordinal) const = 0;
    
    virtual std::size_t find(MultiIndex<D> index) const = 0;
    
    virtual const Slice &slice(std::size_t islice) const = 0;
    
    virtual std::size_t count_slices() const = 0;
    
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
        
        MultiIndex<D> operator*() const
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
    
template<dim_t D, class S>
class SlicedShapeEnumeration : public ShapeEnumeration<D>
{
public:
    class Slice;
    
private:
    S shape_;
    
    bool use_dict_;
    
    std::vector< MultiIndex<D> > table_;
    
    std::valarray<Slice> slices_;
    
public:
    class Slice : public ShapeEnumeration<D>::Slice
    {
        friend class SlicedShapeEnumeration<D,S>;
        
    private:
        const SlicedShapeEnumeration<D,S> *enumeration_;
        
        std::size_t start_;
        std::size_t end_;
        
        std::unordered_map< MultiIndex<D>, std::size_t > dict_;
        
    public:
        typedef typename std::vector<MultiIndex<D>>::const_iterator Iterator;
        
        std::size_t offset() const override
        {
            return start_;
        }
        
        std::size_t size() const override
        {
            return end_ - start_;
        }
        
        Iterator begin() const
        {
            return enumeration_->table_.begin() + start_;
        }
        
        Iterator end() const
        {
            return enumeration_->table_.begin() + end_;
        }
        
        MultiIndex<D> operator[](std::size_t ientry) const override
        {
            assert (ientry < size());
            
            return *(begin() + ientry);
        }
        
        std::size_t find(MultiIndex<D> index) const override
        {
            if (enumeration_->use_dict_) {
                
                auto it = dict_.find(index);
                if (it == dict_.end())
                    return size(); //shape does not contain node
                else
                    return it->second;
                
            } else {
                std::less< MultiIndex<D> > comp;
                
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
    
    MultiIndex<D> operator[](std::size_t ordinal) const override
    {
        assert(ordinal < size());
        
        return table_[ordinal];
    }
    
    typedef typename std::vector<MultiIndex<D>>::const_iterator Iterator;
    
    Iterator begin() const
    {
        return table_.begin();
    }
    
    Iterator end() const
    {
        return table_.end();
    }
    
    const Slice &slice(std::size_t islice) const override
    {
        return slices_[islice];
    }
    
    std::size_t count_slices() const override
    {
        return slices_.size();
    }
    
    std::size_t find(MultiIndex<D> index) const override
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
};

template<dim_t D, class S>
SlicedShapeEnumeration<D,S>::SlicedShapeEnumeration(S shape)
    : shape_(shape)
    , use_dict_(true)
    , table_()
{
    LexicalIndexGenerator<D,S> it(shape);
    
    //holds for each slice its own vector
    std::vector< std::vector< MultiIndex<D> >* > temp;
    
    //generate nodes and put them into the corresponding slice
    std::size_t ordinal = 0;
    do {
        MultiIndex<D> index = it.index();
        
        std::size_t islice = 0;
        for (dim_t i = 0; i < D; i++)
            islice += index[i];
        
        assert(islice <= temp.size()); //assert that lexical index generator does not jump
        
        if (islice == temp.size()) {
            temp.push_back( new std::vector< MultiIndex<D> >() );
        }
        
        temp[islice]->push_back(index);
    } while (it.forward());
    
    //all slices are now populated
    //copy all nodes into the same vector and store offsets
    slices_ = std::valarray<Slice>(temp.size());
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