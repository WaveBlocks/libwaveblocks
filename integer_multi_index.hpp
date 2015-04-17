#ifndef UINT_MULTI_INDEX
#define UINT_MULTI_INDEX

#include <iostream>
#include <vector>
#include <initializer_list>

template<class UINT, dim_t D, std::size_t BITS>
class IntegerMultiIndex
{
private:
    UINT values_;
    
public:
    class Entry
    {
    private:
        UINT &values_;
        std::size_t index_;
        
        inline UINT offset() const
        {
            return BITS*index_;
        }
        
        inline UINT mask() const
        {
            return (1<<BITS)-1;
        }
        
        inline UINT get() const
        {
            return (values_ >> offset()) & mask();
        }
        
        inline void set(UINT value)
        {
            values_ &= ~(mask() << offset());
            values_ |= (value & mask()) << offset();
        }
        
    public:
        Entry(UINT &values, std::size_t index) : values_(values), index_(index) {}
        
        Entry &operator=(UINT value)
        {
            set(value);
            return *this;
        }
        
        Entry &operator=(const Entry &entry)
        {
            set(entry.get());
            return *this;
        }
        
        Entry &operator+=(UINT value)
        {
            value = get() + value;
            set(value);
            return *this;
        }
        
        Entry &operator-=(UINT value)
        {
            value = get() - value;
            set(value);
            return *this;
        }
        
        Entry &operator*=(UINT value)
        {
            value = get() * value;
            set(value);
            return *this;
        }
        
        Entry &operator/=(UINT value)
        {
            value = get() / value;
            set(value);
            return *this;
        }
        
        Entry &operator%=(UINT value)
        {
            value = get() % value;
            set(value);
            return *this;
        }
        
        operator UINT() const
        {
            return get();
        }
    };
    
    IntegerMultiIndex() : values_(0) {}
    IntegerMultiIndex(const IntegerMultiIndex &that) : values_(that.values_) {}
    
    IntegerMultiIndex &operator=(const IntegerMultiIndex &that)
    {
        values_ = that.values_;
        
        return *this;
    }
    
    UINT operator[](int index) const
    {
        return (values_ >> BITS*index) & ( (1<<BITS)-1 );
    }
    
    Entry operator[](int index)
    {
        return Entry(values_, index);
    }
    
    IntegerMultiIndex(std::initializer_list<UINT> list)
    {
        int i = 0;
        for (typename std::initializer_list<UINT>::iterator it = list.begin(); it != list.end() && i < D; it++)
            operator[](i++) = *it;
    }
    
    
    bool operator<(const IntegerMultiIndex &that)
    {
        return values_ < that.values_;
    }
    
    bool operator>(const IntegerMultiIndex &that)
    {
        return values_ > that.values_;
    }
    
    bool operator<=(const IntegerMultiIndex &that)
    {
        return values_ <= that.values_;
    }
    
    bool operator>=(const IntegerMultiIndex &that)
    {
        return values_ >= that.values_;
    }
    
    bool operator==(const IntegerMultiIndex &that)
    {
        return values_ == that.values_;
    }
    
    bool operator!=(const IntegerMultiIndex &that)
    {
        return values_ != that.values_;
    }
};

#endif