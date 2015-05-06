#ifndef WAVEBLOCKS_CONTINUOUS_SQRT
#define WAVEBLOCKS_CONTINUOUS_SQRT

#include <complex>
#include <stdexcept>

namespace waveblocks {

template<class T>
class ContinuousSqrt
{
private:
    /**
     * stored square root
     */
    std::complex<T> sqrt_;
    
    /**
     * argument (angle) of stored square root
     * domain = [0;2*pi)
     */
    T state_;
    
public:
    ContinuousSqrt()
        : sqrt_()
        , state_()
    {}
    
    ContinuousSqrt(const ContinuousSqrt &that) 
        : sqrt_(that.sqrt_)
        , state_(that.state_)
    {}
    
    ContinuousSqrt &operator=(const ContinuousSqrt &that)
    {
        sqrt_ = that.sqrt_;
        state_ = that.state_;
        return *this;
    }
    
    /**
     * retrieve stored square root
     */
    std::complex<T> operator()() const
    {
        return sqrt_;
    }
    
    /**
     * Chooses the square root angle (aka argument) that continuates the reference angle the best.
     * Throws an exception if the deviation above an accepted value (by default >45°) as this
     * strongly indicates a problem in higher level code (for example a too large timestep).
     * \param[in] ref angle of reference root. domain = [0;2*pi)
     * \param[in] arg angle of first(!) root. domain = [0;pi)
     * \return angle of continuating root. domain = [0;2*pi)
     */
    static T continuate(T ref, T arg)
    {
        const T pi = 3.14159265359;
        const T RANGE = 0.25*pi; // 0.5*pi allows all inputs
        
        //determine, how long one needs to 
        //rotate the reference angle counter-clock-wise to hit the angle of the 1st root
        T rot = arg - ref; // domain = (-2pi;pi)
        
        // force rotation into domain [0;2*pi)
        if (rot < 0.0)
            rot += 2.0*pi;
        
        if (rot < RANGE || rot > 2.0*pi-RANGE) {
            return arg;
        }
        else if (rot > pi-RANGE && rot < pi+RANGE) {
            return arg + pi;
        }
        else {
            throw std::runtime_error("too large step");
        }
    }
    
    /**
     * update stored square root
     */
    std::complex<T> operator()(std::complex<T> input)
    {
        const T pi = 3.14159265359;
        const T range = 0.25*pi; // 0.5*pi allows all inputs
        
        state_ = continuate(state_, 0.5*std::arg(input) );
        
        sqrt_ = std::polar(std::sqrt(std::abs(input)), state_);
        
        return sqrt_;
    }
};

}

#endif