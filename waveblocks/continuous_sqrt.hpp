#ifndef WAVEBLOCKS_CONTINUOUS_SQRT
#define WAVEBLOCKS_CONTINUOUS_SQRT

#include <complex>
#include <stdexcept>

#include "math_util.hpp"

namespace waveblocks {

/**
 * \brief This class deals with the issue, that the square root of complex numbers is not unique.
 * 
 * The equation \f$ z^2 = r \exp{(i\phi)} \f$ has two solutions, namely
 * \f$ z_1=\sqrt{r} \exp{\left(i\frac{\phi}{2}\right)} \f$ and 
 * \f$ z_2=\sqrt{r} \exp{\left(i(\frac{\phi}{2}+\pi)\right)} \f$.
 * 
 * This class chooses the solution, that is nearest to the solution of the previous computation ( = reference solution).
 * Then this class overrides the stored reference solution with the current solution.
 * 
 * The distance between the two complex numbers is determined by the angle-distance.
 * 
 * \tparam T Type of both the real and imaginary components of the complex number.
 */
template<class T>
class ContinuousSqrt
{
private:
    /**
     * stored reference solution
     */
    std::complex<T> sqrt_;
    
    /**
     * argument (angle) of reference solution
     * domain = [-pi;pi]
     */
    T state_;
    
    /**
     * false iff a reference solution is stored
     */
    bool empty_;
    
public:
    /**
     * \brief Delayes initialization of the stored reference solution.
     * 
     * The next call to operator()() yields the principal square root.
     */
    ContinuousSqrt()
        : sqrt_()
        , state_()
        , empty_(true)
    { }
    
    /**
     * \brief Initializes the stored reference solution to a chosen value.
     * 
     * \param sqrt The initial reference solution.
     */
    ContinuousSqrt(std::complex<T> sqrt)
        : sqrt_(sqrt)
        , state_(std::arg(sqrt))
        , empty_(false)
    { }
    
    /**
     * Chooses the square root angle (aka argument) that continuates the reference angle the best.
     * Throws an exception if the deviation above an accepted value (by default >45°) as this
     * strongly indicates a problem in higher level code (for example a too large timestep).
     * \param[in] ref The angle of the reference square root. domain = \f$ [-\pi;\pi] \f$
     * \param[in] arg The angle of the computed square root. domain = \f$ [-\pi;\pi] \f$
     * \return The angle of the continuating square root. domain = \f$ [-\pi;\pi] \f$
     */
    static T continuate(T ref, T arg)
    {
        const T PI = pi<T>();
        const T RANGE = 0.25*PI; // 0.5*pi allows all inputs
        
        //determine, how long one needs to rotate
        //the reference angle counter-clock-wise to hit the angle of the 1st root
        T rot = arg - ref; // domain = [-2*pi;2*pi]
        
        // force rotation into domain [-pi;pi]
        if (rot >= PI)
            rot -= 2.0*PI;
        else if (rot < -PI)
            rot += 2.0*PI;
        
        if (rot > -RANGE && rot < RANGE) {
            return arg;
        }
        else if (rot > PI-RANGE || rot < -PI+RANGE) {
            if (arg > 0.0)
                return arg - PI;
            else
                return arg + PI;
        }
        else {
            throw std::runtime_error("continuous_sqrt: too large deviation between computed square root and reference solution");
        }
    }
    
    /**
     * \brief Solves the quadratic equation \f$ z^2 = c \f$, chooses the solution \f$ \hat{z} \f$ that best continuates 
     * the prior result \f$ z_0 \f$ and updates the reference solution (\f$ z_o \gets \hat{z} \f$).
     * 
     * \param input The right-hand-side \f$ c \f$.
     * \return The best solution \f$ \hat{z} \f$.
     */
    std::complex<T> operator()(std::complex<T> input)
    {
        if (empty_) {
            state_ = 0.5*std::arg(input); // choose principal solution
            empty_ = false;
        } else {
            state_ = continuate(state_, 0.5*std::arg(input) );
        }
        
        sqrt_ = std::polar(std::sqrt(std::abs(input)), state_);
        
        return sqrt_;
    }
    
    /**
     * \brief Retrieve the stored reference solution.
     */
    std::complex<T> operator()() const
    {
//         if (empty_)
//             throw std::runtime_error("continuous_sqrt: no reference solution stored");
//         else
        
        return sqrt_;
    }
};

}

#endif