#pragma once

#include <iostream>
#include <array>
#include <cassert>
#include <stdexcept>
#include <string>
#include <initializer_list>
#include <limits>

#include "basic_types.hpp"
#include "stdarray2stream.hpp"


namespace waveblocks {

/**
 * \brief Represents a multi-index using a single integer.
 *
 * This implementation splits an integer into same sized parts.
 * For example, using a 64bit integer to represent a 10-dimensional multi-index:
 * The largest index it can represent is 63 (6 bits).
 * This means when using this class, special care must be taken to prevent _overflows_.
 * Thus, code that uses multi-indices should ensure that a multi-index type is viable:
 *
 * \code
 * MultiIndex mindex;
 * for (dim_t d = 0; d < D; d++) {
 *    mindex[d] = largest_possible_index;
 *    if (mindex[d] != largest_possible_index)
 *       throw std::runtime_error("multi-index type is not compatible!");
 * }
 * \endcode
 *
 * If TinyMultiIndex is not enough for you. You will have to implement your own type.
 * A custom implementation
 * must possess the same semantics as std::array<int,D>.
 * Furthermore it has to specialize std::less, that performs lexical index comparison
 * beginning on the first index.
 * And it has to specialize std::equal_to and std::hash to
 * enable use of multi-indices as hashtable keys.
 *
 * \tparam UINT Option to select, which integer type is used.
 * \tparam D Dimensionality of multi-index i.e number of entries.
 */
template<class UINT, dim_t D>
class TinyMultiIndex
{
    friend struct std::less< waveblocks::TinyMultiIndex<UINT,D> >;
    friend struct std::hash< waveblocks::TinyMultiIndex<UINT,D> >;
    friend struct std::equal_to< waveblocks::TinyMultiIndex<UINT,D> >;

private:

    /**
     * for a given axis: returns the largest value that this implementation is able to store
     */
    static int limit(dim_t axis)
    {
        (void) axis; //supress -Wunused-parameter

        UINT l = (UINT(1)<<BITS_PER_ENTRY)-1;
        if (l > std::numeric_limits<int>::max())
            return std::numeric_limits<int>::max();
        else
            return (int)l;
    }

public:
    UINT values_ = 0;

    const static std::size_t BITS_PER_ENTRY = (8*sizeof(UINT))/D;



    class Entry
    {
    private:
        UINT &values_;

        //care: real index is [index_ - (D-1)]
        std::size_t index_;

        inline UINT offset() const
        {
            return BITS_PER_ENTRY*index_;
        }

        inline UINT mask() const
        {
            return (UINT(1)<<BITS_PER_ENTRY)-1;
        }

        inline int get() const
        {
            return (values_ >> offset()) & mask();
        }

        inline void set(int value)
        {
            assert (value >= 0 && UINT(value) <= mask());

            values_ &= ~(mask() << offset());
            values_ |= (UINT(value) & mask()) << offset();
        }

    public:
        Entry(UINT &values, std::size_t index) : values_(values), index_(index) {}

        Entry &operator=(int value)
        {
            set(value);
            return *this;
        }

        Entry &operator=(const Entry &entry)
        {
            set(entry.get());
            return *this;
        }

        Entry &operator+=(int value)
        {
            value = get() + value;
            set(value);
            return *this;
        }

        Entry &operator-=(int value)
        {
            value = get() - value;
            set(value);
            return *this;
        }

        Entry &operator*=(int value)
        {
            value = get() * value;
            set(value);
            return *this;
        }

        Entry &operator/=(int value)
        {
            value = get() / value;
            set(value);
            return *this;
        }

        Entry &operator%=(int value)
        {
            value = get() % value;
            set(value);
            return *this;
        }

        operator int() const
        {
            return get();
        }
    };

    TinyMultiIndex()
        : values_(0)
    { }

    TinyMultiIndex(const TinyMultiIndex &that)
        : values_(that.values_)
    { }

    TinyMultiIndex(const std::array<int,D> &that)
    {
        for (dim_t d = 0; d < D; d++)
            operator[](d) = that[d];
    }

    TinyMultiIndex(std::initializer_list<int> list)
    {
        dim_t axis = 0;
        for (typename std::initializer_list<int>::iterator it = list.begin(); it != list.end() && axis < D; it++) {
            if (*it > limit(axis))
                throw std::range_error("this multi-index implementation is unable to store a larger value than " + std::to_string(limit(axis)));

            operator[](axis++) = *it;
        }
    }

    TinyMultiIndex &operator=(const TinyMultiIndex &that)
    {
        values_ = that.values_;

        return *this;
    }

    int operator[](dim_t index) const
    {
        return (values_ >> BITS_PER_ENTRY*((D-1)-index)) & ( (UINT(1)<<BITS_PER_ENTRY)-1 );
    }

    Entry operator[](dim_t index)
    {
        return Entry(values_, (D-1)-index);
    }

    bool operator==(const TinyMultiIndex &that) const
    {
        return values_ == that.values_;
    }

    bool operator!=(const TinyMultiIndex &that) const
    {
        return values_ != that.values_;
    }

    operator std::array<int,D>() const
    {
        std::array<int,D> copy;
        for (dim_t d = 0; d < D; d++) {
            copy[d] = operator[](d);
        }
        return copy;
    }
};

template<class UINT, dim_t D>
std::ostream &operator<<(std::ostream &out, const TinyMultiIndex<UINT, D> &index)
{
    std::cout << "(";
    for (dim_t i = 0; i < D-1; i++)
        std::cout << index[i] << ", ";
    if (D != 0)
        std::cout << index[D-1];
    std::cout << ")";
    return out;
}

}

namespace std {
    /**
     * \cond HIDDEN_SYMBOLS
     * Provides less functor (compare) for STL containers (notable std::map).
     * Specializes generic std::less<T>.
     * \endcond
     */
    template<class UINT, waveblocks::dim_t D>
    struct less< waveblocks::TinyMultiIndex<UINT,D> >
    {
    private:
        typedef waveblocks::TinyMultiIndex<UINT,D> MultiIndex;

    public:
        typedef MultiIndex first_argument_type;
        typedef MultiIndex second_argument_type;
        typedef bool result_type;

        bool operator()(const MultiIndex &first, const MultiIndex &second) const
        {
            return first.values_ < second.values_;
        }
    };

    /**
     * \cond HIDDEN_SYMBOLS
     * Provides hash functor for STL containers (notable std::unordered_map).
     * Specializes generic std::hash<T>.
     * \endcond
     */
    template<class UINT, waveblocks::dim_t D>
    struct hash< waveblocks::TinyMultiIndex<UINT,D> >
    {
    private:
        typedef waveblocks::TinyMultiIndex<UINT,D> MultiIndex;

    public:
        std::size_t operator()(const MultiIndex &index) const
        {
            return index.values_;
        }
    };

    /**
     * \cond HIDDEN_SYMBOLS
     * Provides equality functor for STL containers (notable std::unordered_map).
     * Specializes generic std::equal_to<T>.
     * \endcond
     */
    template<class UINT, waveblocks::dim_t D>
    struct equal_to< waveblocks::TinyMultiIndex<UINT,D> >
    {
    private:
        typedef waveblocks::TinyMultiIndex<UINT,D> MultiIndex;

    public:
        typedef MultiIndex first_argument_type;
        typedef MultiIndex second_argument_type;
        typedef bool result_type;

        bool operator()(const MultiIndex &first, const MultiIndex &second) const
        {
            return first.values_ == second.values_;
        }
    };
}
