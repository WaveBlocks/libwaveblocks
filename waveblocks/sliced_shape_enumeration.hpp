#ifndef SLICED_SHAPE_ENUMERATION
#define SLICED_SHAPE_ENUMERATION

#include <vector>

#include "lexical_shape_enumerator.hpp"

namespace waveblocks {

template<std::size_t D, class S>
class SlicedShapeEnumeration
{
public:
    class Slice {
        friend class SlicedShapeEnumeration;
        
    public:
        typedef typename std::vector<MultiIndex<D>>::const_iterator Iterator;
        
        Slice();
        
        std::size_t index() const;
        std::size_t offset() const;
        std::size_t size() const;
        MultiIndex<D> at(std::size_t ordinal) const;
        std::size_t find(MultiIndex<D> index) const;
        SlicedShapeEnumeration<D,S>::Slice::Iterator begin() const;
        SlicedShapeEnumeration<D,S>::Slice::Iterator end() const;
        
    private:
        std::vector<MultiIndex<D>> table_;
        std::size_t offset_;
        std::size_t index_;
    };
    
    SlicedShapeEnumeration(S shape);
    
    std::size_t count() const;
    const SlicedShapeEnumeration<D,S>::Slice &operator[](std::size_t slice);
    
    
    typedef typename std::vector<Slice>::const_iterator Iterator;
    
    SlicedShapeEnumeration<D,S>::Iterator begin() const;
    SlicedShapeEnumeration<D,S>::Iterator end() const;
    
private:
    S shape_;
    std::vector<Slice> slices_;
};

template<std::size_t D, class S>
SlicedShapeEnumeration<D,S>::SlicedShapeEnumeration(S shape)
    : shape_(shape)
    , slices_()
{
    LexicalShapeIterator<D,S> it(shape);
    do {
        MultiIndex<D> index = it.getMultiIndex();
        
        std::size_t slice = 0;
        for (std::size_t i = 0; i < D; i++)
            slice += index[i];
        
        if (slice >= slices_.size())
            slices_.resize(slice+1);
        
        slices_[slice].table_.push_back(index);
    } while (it.advance());
    
    //compute offsets
    std::size_t offset = 0;
    for (std::size_t i = 0; i < slices_.size(); i++) {
        slices_[i].index_ = i;
        slices_[i].offset_ = offset;
        offset += slices_[i].table_.size();
    }
}

template<std::size_t D, class S>
typename SlicedShapeEnumeration<D,S>::Iterator SlicedShapeEnumeration<D,S>::begin() const
{
    return slices_.begin();
}

template<std::size_t D, class S>
typename SlicedShapeEnumeration<D,S>::Iterator SlicedShapeEnumeration<D,S>::end() const
{
    return slices_.end();
}

template<std::size_t D, class S>
SlicedShapeEnumeration<D,S>::Slice::Slice()
    : table_()
    , offset_()
    , index_()
{ }

template<std::size_t D, class S>
std::size_t SlicedShapeEnumeration<D,S>::count() const
{
    return slices_.size();
}

template<std::size_t D, class S>
std::size_t SlicedShapeEnumeration<D,S>::Slice::offset() const
{
    return offset_;
}

template<std::size_t D, class S>
std::size_t SlicedShapeEnumeration<D,S>::Slice::size() const
{
    return table_.size();
}

template<std::size_t D, class S>
std::size_t SlicedShapeEnumeration<D,S>::Slice::index() const
{
    return index_;
}

template<std::size_t D, class S>
typename SlicedShapeEnumeration<D,S>::Slice::Iterator SlicedShapeEnumeration<D,S>::Slice::begin() const
{
    return table_.begin();
}

template<std::size_t D, class S>
typename SlicedShapeEnumeration<D,S>::Slice::Iterator SlicedShapeEnumeration<D,S>::Slice::end() const
{
    return table_.end();
}

template<std::size_t D, class S>
MultiIndex<D> SlicedShapeEnumeration<D,S>::Slice::at(std::size_t ordinal) const
{
    return table_[ordinal];
}

template<std::size_t D, class S>
std::size_t SlicedShapeEnumeration<D,S>::Slice::find(MultiIndex<D> index) const
{
    //do binary search
    std::size_t a = 0, b = table_.size();
    while (b > a+1) {
        std::size_t mid = (b + a)/2;
        if (LexicalShapeIterator<D,S>::compare(table_.at(mid), index) > 0) {
            //table[mid] > index
            b = mid;
        } else {
            a = mid;
        }
    }
    
    if (table_[a] != index)
        return table_.size(); //index not found
    else
        return a;
}

template<std::size_t D, class S>
const typename SlicedShapeEnumeration<D,S>::Slice &SlicedShapeEnumeration<D,S>::operator[](std::size_t slice)
{
    return slices_[slice];
}

}

#endif