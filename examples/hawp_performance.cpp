#include <Eigen/Core>

#include <iostream>
#include <chrono>

#include "waveblocks/tiny_multi_index.hpp"

#include "waveblocks/shape_commons.hpp"

#include "waveblocks/shape_enum.hpp"
#include "waveblocks/shape_enumerator.hpp"

#include "waveblocks/hawp_commons.hpp"

using namespace waveblocks;

class Timer
{
public:
    double millis() const
    {
        return std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1,1000> > >(stop_ - start_).count();
    }
    
    double seconds() const
    {
        return std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1,1> > >(stop_ - start_).count();
    }
    
    void start()
    {
        start_ = std::chrono::high_resolution_clock::now();
    }
    
    void stop()
    {
        stop_ = std::chrono::high_resolution_clock::now();
    }
    
private:
    std::chrono::high_resolution_clock::time_point start_;
    std::chrono::high_resolution_clock::time_point stop_;
};

template<dim_t D>
class PerformanceTest
{
public:
    typedef TinyMultiIndex<std::size_t, D> MultiIndex;
    
    static void run()
    {
        Timer timer;
        Eigen::IOFormat CleanFmt(4, 0, ", ", "\n   ", "[", "]");
        
        ScalarHaWp<D,MultiIndex> wavepacket;
        
        std::cout << "Initialize wavepacket {" << std::endl;
        
        wavepacket.eps() = 0.9;
        wavepacket.parameters().p = RMatrix<D,1>::Random();
        wavepacket.parameters().q = RMatrix<D,1>::Random();
        wavepacket.parameters().P = CMatrix<D,D>::Random();
        wavepacket.parameters().Q = CMatrix<D,D>::Random();
        
        std::array<int,D> limits;
        std::fill(limits.begin(), limits.end(), 4);
        int sparsity = 1 << D;
        
        LimitedHyperbolicCutShape<D> shape(sparsity, limits);
        
        ShapeEnumerator<D,MultiIndex> enumerator;
        
        wavepacket.shape() = enumerator.enumerate(shape);
        wavepacket.coefficients() = std::vector<complex_t>(wavepacket.shape()->n_entries(), 1.0/wavepacket.shape()->n_entries());
        
        std::cout << "   dimensionality: " << D << std::endl;
        std::cout << "   shape: " << shape << std::endl;
        std::cout << "   number of basis functions: " << wavepacket.shape()->n_entries() << std::endl;
        std::cout << "}" << std::endl;
        
        
        
        timer.start();
        
        std::cout << "Evaluate wavepacket {" << std::endl;
        RMatrix<D, Eigen::Dynamic> grid(D,1);
        for (int i = 0; i < 1; i++) {
            grid(i,0) = -1.0 + 2.0*(i-1)/D;
        }
        
        CMatrix<1, Eigen::Dynamic> result = wavepacket.evaluate(grid);

        timer.stop();
        
        std::cout << "   value: " << result.format(CleanFmt) << std::endl;
        std::cout << "   time:  " << timer.millis() << " [ms]" << std::endl;
        std::cout << "}" << std::endl;
    }
    
private:
    
};

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    
    PerformanceTest<5>::run();
    PerformanceTest<8>::run();
    PerformanceTest<11>::run();
    PerformanceTest<14>::run();
    
    return 0;
}