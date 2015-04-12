#ifndef WAVEBLOCKS_LEXICAL_SHAPE_ENUMERATION
#define WAVEBLOCKS_LEXICAL_SHAPE_ENUMERATION

#include <memory>
#include <cassert>

namespace waveblocks {

template<dim_t D>
class LexicalShapeEnumeration
{
private:
    const Shape<D> *shape_;
    
    std::size_t stride_;
    std::size_t size_;
    std::shared_ptr<std::vector<MultiIndex<D>>> table_;
    
public:
    class Iterator
    {
    private:
        const LexicalShapeEnumeration *enumeration_;
        
        std::size_t ordinal_;
        MultiIndex<D> index_;
        
    public:
        Iterator(const LexicalShapeEnumeration *enumeration_, const MultiIndex<D> &index, std::size_t ordinal) : 
                enumeration_(enumeration_), 
                ordinal_(ordinal), 
                index_(index) {}
        
        Iterator(const Iterator &that) : 
                enumeration_(that.enumeration_), 
                ordinal_(that.ordinal_), 
                index_(that.index_) {}
        
        Iterator &operator=(const Iterator &that)
        {
            enumeration_ = that.enumeration_;
            ordinal_ = that.ordinal_;
            index_ = that.index_;
            return *this;
        }
        
        bool operator==(const Iterator &that) const
        {
            //dont compare multi-index since the mapping of ordinal to multi-index is bijective
            return enumeration_ == that.enumeration_ && ordinal_ == that.ordinal_;
        }
        
        bool operator!=(const Iterator &that) const
        {
            return !operator==(that);
        }
        
        std::size_t getOrdinal() const
        {
            return ordinal_;
        }
        
        MultiIndex<D> getMultiIndex() const
        {
            return index_;
        }
        
        bool getForwardNeighbour(dim_t axis, Iterator &it) const
        {
            int surface = enumeration_->shape_->surface(axis, index_);
            if (index_[axis] == surface)
                return false;
            
            it = *this;
            
            it.index_[axis] += 1;
            
            if (axis == 0)
                it.ordinal_ += 1;
            else if (axis == 1)
                it.ordinal += surface + 1;
            else
                enumeration_->find(it.index_, it.ordinal_);
            
            return true;
        }
        
        bool getBackwardNeighbour(dim_t axis, Iterator &it) const
        {
            int surface = enumeration_->shape_->surface(axis, index_);
            if (index_[axis] == 0)
                return false;
            
            it = *this;
            
            it.index_[axis] -= 1;
            
            if (axis == 0)
                it.ordinal_[0] -= 1;
            else if (axis == 1)
                it.ordinal -= surface + 1;
            else
                enumeration_->find(it.index_, it.ordinal_);
            
            return true;
        }
        
        bool getNext(Iterator &it) const
        {
            it = *this;
            
            dim_t axis = 0;
            while ((int)it.index_[axis] == enumeration_->shape_->surface(axis, it.index_)) {
                if (axis == D-1)
                    return false;
                it.index_[axis++] = 0;
            }
            
            it.index_[axis] += 1;
            
            ++it.ordinal_;
            
            return true;
        }
    };
    
    LexicalShapeEnumeration(const Shape<D> *shape, std::size_t stride) : 
            shape_(shape), stride_(stride == 0 ? 1 : stride), table_(std::make_shared<std::vector<MultiIndex<D>>>())
    {
        //create tableau
        Iterator it = begin();
        do {
            if (it.getOrdinal()%stride == 0) {
                table_->push_back(it.getMultiIndex());
            }
        } while (it.getNext(it));
        
        size_ = 1 + it.getOrdinal();
    }
    
    LexicalShapeEnumeration(const LexicalShapeEnumeration &that) : 
            shape_(that.shape_), stride_(that.stride_), table_(that.table_), size_(that.size_) {}
    
    LexicalShapeEnumeration &operator=(const LexicalShapeEnumeration &that)
    {
        shape_ = that.shape_;
        stride_ = that.stride_;
        table_ = that.table_;
        size_ = that.size_;
        
        return *this;
    }
    
    std::size_t size() const
    {
        return size_;
    }
    
    Iterator begin() const
    {
        return Iterator(this, MultiIndex<D>{}, 0);
    }
    
    Iterator end() const
    {
        return Iterator(this, MultiIndex<D>{}, size());
    }
    
    /**
     * computes a - b
     */
    int compare(const MultiIndex<D> &a, const MultiIndex<D> &b) const
    {
        for (dim_t i = D-1; i >= 0; i--) {
            int diff = a[i] - b[i];
            if (diff != 0)
                return diff;
        }
        return 0;
    }
    
    bool contains(const MultiIndex<D> &index) const
    {
        for (dim_t i = 0; i < D; i++)
            if (index[i] > shape_.surface(i, index))
                return false;
        return true;
    }
    
    bool find(const MultiIndex<D> &index, Iterator &it) const
    {
        if (!contains(index))
            return false;
        
        //binary-search for biggest entry that is smaller or equals requested multi-index
        MultiIndex<D> *table = table_->data();
        
        //Note:
        // table[a] <= index
        // table[b] > index
        std::size_t a = 0, b = table_->size();
        while (b > a+1) {
            std::size_t mid = (b + a)/2;
            
            if (compare(table[mid], index) > 0) {
                //selected entry is bigger than requested entry
                b = mid;
            } else {
                a = mid;
            }
        }
        
        //
        it = Iterator(this, table[a], a*stride_);
        
        while (it.getMultiIndex() != index) {
            bool result = it.getNext(it);
            assert(result);
        }
        
        return true;
    }
    
    bool at(std::size_t ordinal, Iterator &it) const
    {
        if (ordinal >= size_)
            return false;
        
        std::size_t major = ordinal/stride_;
        
        it = Iterator(*this, (*table_)[major], major*stride_);
        
        std::size_t remain = ordinal - major*stride_;
        while (remain-- > 0) {
            bool result = it.getNext(it);
            assert(result);
        }
        
        return true;
    }
};

}

#endif