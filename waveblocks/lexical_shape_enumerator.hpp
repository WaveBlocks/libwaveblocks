#ifndef WAVEBLOCKS_LEXICAL_SHAPE_ENUMERATOR_HPP
#define WAVEBLOCKS_LEXICAL_SHAPE_ENUMERATOR_HPP

#include <cassert>

#include "multi_index.hpp"

namespace waveblocks {

template<std::size_t D, class S>
class LexicalShapeIterator
{
private:
    S shape_;
    
    std::size_t ordinal_;
    MultiIndex<D> index_;
    
public:
    LexicalShapeIterator(S shape) : shape_(shape), ordinal_(), index_()
    {
        
    }
    
    LexicalShapeIterator(S shape, std::size_t ordinal, MultiIndex<D> index)
        : shape_(shape), ordinal_(ordinal), index_(index)
    {
        
    }
    
    LexicalShapeIterator(const LexicalShapeIterator &that)
        : shape_(that.shape_), ordinal_(that.ordinal_), index_(that.index_)
    {
        
    }
    
    LexicalShapeIterator &operator=(const LexicalShapeIterator &that)
    {
        shape_ = that.shape_;
        ordinal_ = that.ordinal_;
        index_ = that.index_;
        
        return *this;
    }
    
    MultiIndex<D> getMultiIndex() const
    {
        return index_;
    }
    
    std::size_t getOrdinal() const
    {
        return ordinal_;
    }
    
    bool operator==(const LexicalShapeIterator<D,S> &that) const
    {
        //dont check index equality as index-ordinal mapping is bijective
        //dont check shape equality to improve performance
        return ordinal_ == that.ordinal_;
    }
    
    bool operator!=(const LexicalShapeIterator<D,S> &that) const
    {
        return this->operator==(that);
    }
    
    bool advance()
    {
        ++ordinal_;
        
        std::size_t axis = 0;
        while ((int)index_[axis] == shape_.getSurface(axis, index_)) {
            index_[axis++] = 0;
            if (axis == D) {
                //end of chain
                return false;
            }
        }
        
        index_[axis] += 1;
        
        return true;
    }
    
    /**
     * return value corresponds to (A - B)
     * 
     * 
     * returns 0 when first is equals second
     * returns positive value when first is greater than second
     * returns negative value when first is smaller than second
     */
    static int compare(MultiIndex<D> a, MultiIndex<D> b)
    {
        for (std::size_t i = D-1; i >= 0; i--) {
            int diff = a[i] - b[i];
            if (diff != 0)
                return diff;
        }
        return 0;
    }
};

template<std::size_t D, class S>
class LexicalShapeEnumeration
{
private:
    S shape_;
    std::shared_ptr<std::vector<MultiIndex<D>>> table_;
    
    std::size_t size_;
    std::size_t stride_;
    
public:
    LexicalShapeEnumeration(S shape, std::size_t stride) : 
            shape_(shape), 
            table_(std::make_shared(std::vector<MultiIndex<D>>{})), 
            stride_(stride == 0 ? 1 : stride)
    {
        LexicalShapeIterator<D,S> it(shape);
        do {
            if (it.getOrdinal()%stride_ == 0)
                table_->push_back(it.getMultiIndex());
        } while (it.advance());
        size_ = it.getOrdinal();
    }
    
    LexicalShapeEnumeration(const LexicalShapeEnumeration &that) : 
            shape_(that.shape_), 
            table_(that.table_), 
            size_(that.size_), 
            stride_(that.stride_)
    {
        
    }
    
    LexicalShapeEnumeration &operator=(const LexicalShapeEnumeration &that)
    {
        shape_ = that.shape_;
        table_ = that.table_;
        size_ = that.size_;
        stride_ = that.stride_;
        return *this;
    }
    
    std::size_t size() const
    {
        return size_;
    }
    
    LexicalShapeIterator<D,S> begin() const
    {
        return LexicalShapeIterator<D,S>(shape_, MultiIndex<D>{}, 0);
    }
    
    LexicalShapeIterator<D,S> end() const
    {
        return LexicalShapeIterator<D,S>(shape_, MultiIndex<D>{}, size_);
    }
    
    LexicalShapeIterator<D,S> find(MultiIndex<D> index) const
    {
        //binary-search for biggest entry that is smaller or equals requested multi-index
        
        //Note:
        // table[a] <= index
        // table[b] > index
        std::size_t a = 0, b = table_->size();
        while (b > a+1) {
            std::size_t mid = (b + a)/2;
            
            if (LexicalShapeIterator<D,S>::compare(table_->at(mid), index) > 0) {
                //selected entry is bigger than requested entry
                b = mid;
            } else {
                a = mid;
            }
        }
        
        //
        LexicalShapeIterator<D,S> it(shape_, a*stride_, table_->at(a));
        
        //Note: this loop breaks if iterator reaches end
        while (it.getMultiIndex() != index && it.advance()) {}
    }
    
    LexicalShapeIterator<D,S> at(std::size_t ordinal) const
    {
        if (ordinal >= size_)
            return end();
        
        std::size_t major = ordinal/stride_;
        
        LexicalShapeIterator<D,S> it(shape_, table_->at(major), major*stride_);
        
        std::size_t minor = ordinal%stride_;
        
        while (minor-- > 0)
        {
            bool result = it.advance();
            assert (result);
        }
    }
};

}

#endif