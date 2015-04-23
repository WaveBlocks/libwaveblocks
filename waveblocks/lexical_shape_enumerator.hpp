#ifndef WAVEBLOCKS_LEXICAL_SHAPE_ENUMERATOR_HPP
#define WAVEBLOCKS_LEXICAL_SHAPE_ENUMERATOR_HPP

#include <cassert>
#include <memory>
#include <string>
#include <sstream>
#include <iostream>

#include "multi_index.hpp"

namespace waveblocks {

template<dim_t D, class S>
class LexicalIndexGenerator
{
private:
    S shape_;
    
    MultiIndex<D> index_;
    
public:
    LexicalIndexGenerator(S shape)
        : shape_(shape)
        , index_()
    { }
    
    LexicalIndexGenerator(S shape, std::size_t ordinal, MultiIndex<D> index)
        : shape_(shape)
        , index_(index)
    { }
    
    LexicalIndexGenerator(const LexicalIndexGenerator &that)
        : shape_(that.shape_)
        , index_(that.index_)
    { }
    
    LexicalIndexGenerator &operator=(const LexicalIndexGenerator &that)
    {
        shape_ = that.shape_;
        index_ = that.index_;
        
        return *this;
    }
    
    MultiIndex<D> index() const
    {
        return index_;
    }
    
    bool backward()
    {
        //find first non-zero entry
        std::size_t axis = D-1;
        while (index_[axis] == 0) {
            if (axis == 0)
                return false;
            else
                --axis;
        }
        
        index_[axis++] -= 1;
        for (; axis < D; axis++) {
            index_[axis] = shape_.getSurface(axis, index_);
        }
        
        return true;
    }
    
    bool forward()
    {
        std::size_t axis = D-1;
        while ((int)index_[axis] == shape_.getSurface(axis, index_)) {
            index_[axis] = 0;
            if (axis == 0)
                return false; //end of chain
            else
                --axis;
        }
        
        index_[axis] += 1;
        
        return true;
    }
};

template<dim_t D, class S>
class LexicalShapeEnumeration
{
private:
    S shape_;
    std::shared_ptr<std::vector<MultiIndex<D>>> table_;
    
public:
    typedef typename std::vector<MultiIndex<D>>::const_iterator Iterator;
    
    LexicalShapeEnumeration(S shape)
        : shape_(shape)
        , table_(std::make_shared<std::vector<MultiIndex<D>>>())
    {
        LexicalIndexGenerator<D,S> gen(shape);
        do {
            table_->push_back(gen.index());
        } while (gen.forward());
    }
    
    LexicalShapeEnumeration(const LexicalShapeEnumeration &that) 
        : shape_(that.shape_)
        , table_(that.table_)
    { }
    
    LexicalShapeEnumeration &operator=(const LexicalShapeEnumeration &that)
    {
        shape_ = that.shape_;
        table_ = that.table_;
        return *this;
    }
    
    std::size_t size() const
    {
        return table_->size();
    }
    
    Iterator begin() const
    {
        return table_->begin();
    }
    
    Iterator end() const
    {
        return table_->end();
    }
    
    std::size_t find(MultiIndex<D> index) const
    {
        LexicalMultiIndexCompare<D> comp;
        
        auto it = std::lower_bound(table_->begin(), table_->end(), index, comp);
        
        if (*it == index)
            return it - table_->begin();
        else
            return table_->size(); //not found
    }
    
    MultiIndex<D> operator[](std::size_t ordinal) const
    {
        assert (ordinal < table_.size());
        
        return table_->at(ordinal);
    }
};

}

#endif