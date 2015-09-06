#ifndef WAVEBLOCKS_KAHAN_SUM_HPP
#define WAVEBLOCKS_KAHAN_SUM_HPP

namespace waveblocks
{

/**
 * \brief The Kahan's algorithm achieves O(1) error growth for summing N numbers.
 * 
 * This implementation works for every type that satisfies following requirements:
 *  - copy constructor
 *  - copy assignment operator
 *  - binary plus operator that performs element-wise addition
 *  - binary minus operator that performs element-wise subtraction
 * 
 * \tparam T type of summands
 */
template<class T>
class KahanSum
{
private:
    /**
     * Accumulated sum.
     */
    T sum_;
    
    /**
     * A running compensation for lost low-order bits.
     */
    T c_;
    
    /**
     * Preallocated temporary variables
     */
    T y_;
    T t_;
    
public:
    /**
     * \brief compiler-generated default constructor
     * 
     * Be aware that this constructor fails on summand types that cannot 
     * provide a default constructor like dynamically sized arrays/matrices.
     * 
     * \see KahanSum::KahanSum(const T&)
     */
    KahanSum() = default;
    
    /**
     * \brief A zero initializing constructor for types that cannot provide a default constructor.
     * 
     * A noteable example is a dynamically sized matrix.
     * 
     * \param[in] zero initialized zero value as a template
     */
    KahanSum(const T &zero)
        : sum_(zero)
        , c_(zero)
    { }
    
    /**
     * \brief compiler-generated copy constructor
     */
    KahanSum(const KahanSum<T> &that) = default;
    
    /**
     * \brief compiler-generated copy assignment operator
     */
    KahanSum &operator=(const KahanSum<T> &that) = default;
    
    /**
     * \brief adds a number 
     * 
     * \param[in] summand summand
     */
    KahanSum &operator+=(const T &summand)
    {
        y_ = summand - c_;
        t_ = sum_ + y_;
        c_ = (t_ - sum_) - y_;
        sum_ = t_;
        return *this;
    }
    
    /**
     * \brief retrieves accumulated sum.
     * 
     * \return accumulated sum
     */
    const T &operator()() const
    {
        return sum_;
    }
};

}

#endif