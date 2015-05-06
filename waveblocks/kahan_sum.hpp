#ifndef WAVEBLOCKS_KAHAN_SUM_HPP
#define WAVEBLOCKS_KAHAN_SUM_HPP

namespace waveblocks
{

/**
 * The Kahan's algorithm achieves O(1) error growth for summing N numbers.
 * 
 * This implementation works for complex numbers too.
 */
template<class T>
class KahanSum
{
private:
    T sum_;
    
    /**
     * A running compensation for lost low-order bits.
     */
    T c_;
    
public:
    /**
     * Zero initialize
     */
    KahanSum()
        : sum_()
        , c_()
    { }
    
    /**
     * Copy constructor
     */
    KahanSum(const KahanSum<T> &that)
        : sum_(that.sum_)
        , c_(that.c_)
    { }
    
    /**
     * Assignment operator
     */
    KahanSum &operator=(const KahanSum<T> &that)
    {
        c_ = that.c_;
        sum_ = that.sum_;
        return *this;
    }
    
    /**
     * Add a number
     */
    KahanSum &operator+=(T input)
    {
        T y = input - c_;
        T t = sum_ + y;
        c_ = (t - sum_) - y;
        sum_ = t;
        return *this;
    }
    
    /**
     * Retrieve accumulated sum
     */
    T operator()() const
    {
        return sum_;
    }
};

}

#endif