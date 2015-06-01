#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <functional>
#include <cassert>

template<class T, std::size_t S>
std::ostream &operator<<(std::ostream &out, const std::array<T,S> &index)
{
    std::cout << "(";
    for (std::size_t i = 0; i < S-1; i++)
        std::cout << index[i] << ", ";
    if (S != 0)
        std::cout << index[S-1];
    std::cout << ")";
    return out;
}

template<int D>
void enumerate(std::array<int,D> &index, int depth, int acc_sum, int tot_sum)
{
    
    if (depth == D-1) {
        index[depth] = tot_sum - acc_sum;
        std::cout << index << std::endl;
        index[depth] = 0;
    }
    else if (acc_sum == tot_sum) {
        std::cout << index << std::endl;
    }
    else {
        for (int i = 0; i <= tot_sum - acc_sum; i++) {
            index[depth] = i;
            enumerate<D>(index, depth+1, acc_sum+i, tot_sum);
        }
        index[depth] = 0;
    }
}

    
//     NodeSize build_(const std::vector<NodeSize>& data, std::size_t begin, std::array<int,D> &index, int depth, int acc_sum, int tot_sum)
//     {
//         if (depth == D-1) {
//             index[depth] = tot_sum - acc_sum;
//             acc_sum = tot_sum;
//         }
//         
//         NodeSize size;
//         if (acc_sum == tot_sum) {
//             // write leaf
//             size.length = sizeof(NodeSize);
//             size.n_leafs = 1;
//         }
//         else if (empty(index, depth)) {
//             // write empty subtree
//             size.length = sizeof(NodeSize);
//             size.n_leafs = 0;
//         }
//         else {
//             // write node & subtrees
//             size.length = sizeof(NodeSize);
//             size.n_leafs = 0;
//             
//             //write subtrees
//             for (int i = 0; i <= tot_sum - acc_sum; i++) {
//                 index[depth] = i;
//                 size += build_(data, begin + size.length, index, depth+1, acc_sum+i, tot_sum);
//             }
//         }
//         data[begin] = size;
//         
//         index[depth] = 0;
//     }

template<int D>
struct HyperbolicCutShape
{
private:
    int K_;
    
public:
    HyperbolicCutShape(int K)
        : K_(K)
    { }
    
    struct stack_entry_type
    {
        int product = 1;
    };
    
    inline bool empty(const std::array<int,D> &index, int dim, int value, stack_entry_type &top, int remain) const
    {
        // update accumulated product
        top.product *= value+1;
        
        // 
        int lowerbound = top.product * (remain+1);
        
        //std::cout << dim << ": " << index << " : " << top.product << " lowerbound: " << lowerbound << std::endl;
        
        return lowerbound > K_;
    }
    
    inline bool contains(const std::array<int,D> &index)
    {
        int prod = 1;
        for (int i = 0; i < D; i++) {
            prod *= index[i]+1;
        }
        return prod < K_;
    }
};

template<int D>
struct HyperCubicShape
{
private:
    std::array<int,D> bbox_;
    
public:
    struct stack_entry_type { };
    
    inline bool empty(const std::array<int,D> &index, int dim, int value, stack_entry_type &top, int remain) const
    {
        return value > bbox_[dim];
    }
};

template<int D, class S>
struct ExtendedShape
{
private:
    S nested_;
    
public:
    ExtendedShape(const S& nested)
        : nested_(nested)
    { }
    
    struct stack_entry_type
    {
        typename S::stack_entry_type previous;
        typename S::stack_entry_type current;
        
        stack_entry_type()
            : previous()
            , current()
        { }
    };
    
    inline bool empty(const std::array<int,D> &index, int dim, int value, stack_entry_type &top, int remain) const
    {
        bool empty_prev = nested_.empty(index, dim, value-1, top.previous, remain);
        bool empty_curr = value == 0 ? false : nested_.empty(index, dim, value-1, top.previous, remain);
        
        if (empty_curr) {
            // check whether neighbour in previous slice is feasible
            if (value > 0)
                return empty_prev;
            else
                return true;
        }
        else {
            assert(!empty_prev);
            
            // nested shape contains at least one entry => extended shape too
            return true;
        }
    }
    
    inline bool contains(std::array<int,D> index)
    {
        //check whether nested shape already contains node
        if (nested_.contains(index))
            return true;
        
        //check whether nested shape contains at least one neighbor
        for (int i = 0; i < D; i++) {
            if (index[i] > 0) {
                index[i] -= 1;
                if (nested_.contains(index))
                    return true;
                index[i] += 1;
            }
        }
        
        return false;
    }
};

template<int D, class S>
void enumerate(const S& shape,
               const std::function< void(const std::array<int,D>&) > &func,
               typename S::stack_entry_type user,
               std::array<int,D> &index,
               int depth, int remain)
{
    if (depth == D-1) {
        index[depth] = remain;
        remain = 0;
    }
    
    if (remain == 0) {
        // leaf
        func(index);
    }
    else {
        for (int i = 0; i <= remain; i++) {
            index[depth] = i;
            
            typename S::stack_entry_type top = user;
            if (!shape.empty(index, depth, i, top, remain-i)) {
                enumerate<D,S>(shape, func, top, index, depth+1, remain-i);
            } else {
                std::cout << "filtered: " << depth << " " << index << std::endl;
            }
        }
    }
    
    index[depth] = 0;
}

template<int D, class S>
void enumerate(const S& shape, const std::function< void(const std::array<int,D>&) > &func, int slice)
{
    std::array<int,D> index{};
    
    typename S::stack_entry_type user{};
    
    enumerate<D,S>(shape, func, user, index, 0, slice);
}

int main()
{
    //SparseManhattanTree<3> tree(5);
    
    const int D = 5;
    
    typedef HyperbolicCutShape<D> S;
    typedef ExtendedShape<D, HyperbolicCutShape<D> > ES;
    
    double K = 11.0;
    
    S shape(K);
    
    ES eshape(shape);
    
    enumerate<D,S>(shape, [&shape,K](const std::array<int,D>& index) { 
        int prod = 1;
        for (int i = 0; i < D; i++)
            prod *= index[i]+1;
        
        if (!shape.contains(index))
            std::cout << "fail: " << index << std::endl;
    }, 5);
    
    return 0;
}