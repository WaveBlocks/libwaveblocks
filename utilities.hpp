#pragma once
#include<vector>
#include<functional>
#include"types.hpp"

template<class I, int N, int M>
class Adaptor;

template<class I,int N, int M>
class AdaptorIterator {
private:
	Adaptor<I,N,M> g;
	int i;
	const Adaptor<I,N,M>& adaptor;
	
public:
	AdaptorIterator(const Adaptor<I,N,N>& adaptor, int i) : adaptor(adaptor), i(i) {}
	bool operator!=(AdaptorIterator<I,N,M> other) {return other.i != i;}
	
	void operator++() {++i;}
	GVector<I,N>& operator*() {
		return g[i];
		
	}
	
};

template<class I, int N, int M>
class Adaptor {
private:
	const GMatrix<I,N,M>& matrix;
	
public:
	size_t size(){return M;}

	Adaptor(size_t size) {}
	Adaptor(const GMatrix<I,N,M>& matrix) : matrix(matrix) {}
	
	GVector<I,N>& operator[] (int i) {
		return matrix.block<N,1>(i,0);
	}
	
	
	AdaptorIterator<I,N,M> begin() {
		return AdaptorIterator<I,N,M>(*this, 0);
	}
	AdaptorIterator<I,N,M> end() {
		return AdaptorIterator<I,N,M>(*this, M);
	}
		
	
};

template<int N, int M>
struct AdaptorHelper {
	template<class I>
	using type = Adaptor<I,N,M>;
};


// Function evaluations
template<class A, class R, template<typename ...> class G = std::vector, template<typename ...> class F = std::function>
G<R> evaluate_function_in_grid(F<R(A)> f, G<A> g)
{
    G<R> result(g.size());
    auto it = result.begin();
    for (auto& arg: g) {
        *it = f(arg);
        ++it;
    }
    return result;
}

template<int N, template<typename, int> class V, class A, class R, template<typename ...> class F = std::function>
V<R,N> evaluate_function_vector(V<F<R(A)>,N> mf, A arg)
{
    V<R,N> m;
    for (int i = 0; i < N; ++i){
		m[i] = mf[i](arg);
    }
	
	return m;
}

template<int N, template<typename, int> class V, class A, class R, template<typename ...> class G = std::vector, template<typename ...> class F = std::function>
G<V<R,N>> evaluate_function_vector_in_grid(V<F<R(A)>,N> mf, G<A> g)
{
	G<V<R,N>> result(g.size());
	auto it = result.begin();
	for (auto& arg: g) {
		*it = evaluate_function_vector<N,V,A,R,F>(mf, arg);
		++it;
	}
	return result;
}

template<int N, template<typename, int, int> class M, class A, class R, template<typename ...> class F = std::function>
M<R,N,N> evaluate_function_matrix(M<F<R(A)>,N,N> mf, A arg)
{
    M<R,N,N> m;
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j) {			
            m(i,j) = mf(i,j)(arg);
        }
    }
	
	return m;
}

template<int N, template<typename, int, int> class M, class A, class R, template<typename ...> class G = std::vector, template<typename ...> class F = std::function>
G<M<R,N,N>> evaluate_function_matrix_in_grid(M<F<R(A)>,N,N> mf, G<A> g)
{
	G<M<R,N,N>> result(g.size());
	auto it = result.begin();
	for (auto& arg: g) {
		*it = evaluate_function_matrix<N,M,A,R,F>(mf, arg);
		++it;
	}
	return result;
}
