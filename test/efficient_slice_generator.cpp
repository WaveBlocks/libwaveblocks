#include <iostream>
#include <array>
#include <unordered_map>

#include "waveblocks/multi_index.hpp"
#include "waveblocks/hyperbolic_shape.hpp"

using namespace waveblocks;

template<int D, class S, int I>
struct Loop
{
    Loop<D,S,I+1> inner;

    int total;
    const S &shape;

    Loop(const S &shape, int total)
        : inner(shape, total)
        , total(total)
        , shape(shape)
    { }

    void operator()(MultiIndex<D> &index, int sum) const
    {
        int limit = shape.getSurface(I-1,index);

        for (int ki = 0; ki <= std::min(limit, total - sum); ki++) {
            index[I-1] = ki;
            inner(index, sum + ki);
        }
        index[I-1] = 0;
    }
};

template<int D, class S>
struct Loop<D,S,D>
{
    int total;
    const S &shape;

    Loop(const S &shape, int total)
        : total(total), shape(shape)
    { }

    void operator()(MultiIndex<D> &index, int sum) const
    {
        int kd = total - sum;
        if (kd <= shape.getSurface(D-1,index)) {
            index[D-1] = kd;
            std::cout << index << std::endl;
            index[D-1] = 0;
        }
    }
};

int main()
{
    const int D = 4;

    typedef HyperbolicCutShape<D> S;

    S shape(7.0);

    Loop<D,S,1> loop(shape, 2);

    MultiIndex<D> index {};
    loop(index, 0);

    return 0;
}
