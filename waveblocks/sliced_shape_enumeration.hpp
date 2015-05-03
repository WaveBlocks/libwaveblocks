#ifndef SLICED_SHAPE_ENUMERATION
#define SLICED_SHAPE_ENUMERATION

#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_map>

#include "lexical_shape_enumerator.hpp"

namespace waveblocks {

template<dim_t D, class S>
class SlicedShapeEnumeration
{
private:
    S shape_;
    
    bool use_dict_;
    
    std::vector< MultiIndex<D> > table_;
    
    /**
     * stores offset of every slice
     * size: #slices + 1
     * first entry: 0
     * last entry: table_.size()
     */
    std::vector<std::size_t> offsets_;
    
    std::vector< std::unordered_map< MultiIndex<D>, std::size_t > > dict_;
    
public:
    class Slice
    {
    private:
        const SlicedShapeEnumeration<D,S> *enumeration_;
        std::size_t islice_;
        
    public:
        Slice(const SlicedShapeEnumeration<D,S> *enumeration, std::size_t islice)
            : enumeration_(enumeration)
            , islice_(islice)
        { }
        
        Slice(const Slice &that)
            : enumeration_(that.enumeration_)
            , islice_(that.islice_)
        { }
        
        Slice &operator=(const Slice &that)
        {
            enumeration_ = that.enumeration_;
            islice_ = that.islice_;
            
            return *this;
        }
        
        typedef typename std::vector<MultiIndex<D>>::const_iterator Iterator;
        
        std::size_t offset() const
        {
            return enumeration_->offsets_[islice_];
        }
        
        std::size_t size() const
        {
            return enumeration_->offsets_[islice_+1] - enumeration_->offsets_[islice_];
        }
        
        Iterator begin() const
        {
            return enumeration_->table_.begin() + enumeration_->offsets_[islice_];
        }
        
        Iterator end() const
        {
            return enumeration_->table_.begin() + enumeration_->offsets_[islice_+1];
        }
        
        MultiIndex<D> operator[](std::size_t ientry) const
        {
            assert (ientry < size());
            
            return *(begin() + ientry);
        }
        
        std::size_t find(MultiIndex<D> index) const
        {
            if (enumeration_->use_dict_) {
                
                auto & slice_dict = enumeration_->dict_.at(islice_);
                
                auto it = slice_dict.find(index);
                if (it == slice_dict.end())
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
    
    class Slices
    {
    private:
        const SlicedShapeEnumeration<D,S> *enumeration_;
        
    public:
        Slices(const SlicedShapeEnumeration<D,S> *enumeration)
            : enumeration_(enumeration)
        { }
        
        Slice operator[](std::size_t islice)
        {
            return Slice(enumeration_, islice);
        }
        
        std::size_t count() const
        {
            return enumeration_->offsets_.size() - 1;
        }
    };
    
    SlicedShapeEnumeration(S shape);
    
    std::size_t size() const
    {
        return table_.size();
    }
    
    MultiIndex<D> operator[](std::size_t ordinal) const
    {
        assert(ordinal < size());
        
        return table_.at(ordinal);
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
    
    Slice slice(std::size_t islice) const
    {
        return Slice(this, islice);
    }
    
    Slices slices() const
    {
        return Slices(this);
    }
    
    std::size_t find(MultiIndex<D> index) const
    {
        std::size_t islice = 0;
        for (dim_t i = 0; i < D; i++)
            islice += index[i];
        
        if (islice >= slices().count())
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
    , use_dict_(false)
    , table_()
    , offsets_()
    , dict_()
{
    LexicalIndexGenerator<D,S> it(shape);
    
    std::vector< std::vector<MultiIndex<D>> > slices;
    std::size_t ordinal = 0;
    do {
        MultiIndex<D> index = it.index();
        
        std::size_t islice = 0;
        for (dim_t i = 0; i < D; i++)
            islice += index[i];
        
        assert(islice <= slices.size()); //assert that lexical index generator does not jump
        
        if (islice == slices.size()) {
            slices.emplace_back();
        }
        
        slices[islice].push_back(index);
    } while (it.forward());
    
    //copy all slices into same vector and store offsets
    for (std::size_t i = 0; i < slices.size(); i++) {
        offsets_.push_back(table_.size());
        table_.insert(table_.end(), slices[i].begin(), slices[i].end());
        
        slices[i].clear();
        slices[i].shrink_to_fit(); //force memory release
    }
    offsets_.push_back(table_.size());
    
    //create dictionary
    if (use_dict_) {
        for (std::size_t islice = 0; islice < slices.size(); islice++) {
            dict_.emplace_back();
            
            auto & slicedict = dict_.back();
            
            std::size_t ordinal = 0;
            for (auto index : slice(islice)) {
                slicedict[index] = ordinal++;
            }
        }
    }
}

}

#endif