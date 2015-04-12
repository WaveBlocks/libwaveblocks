#ifndef WAVEBLOCKS_HYPERCUBIC_SHAPE
#define WAVEBLOCKS_HYPERCUBIC_SHAPE

namespace WaveBlocks {

template<dim_t D>
class HyperCubicShape
{
private:
    /**
     * largest multi-index
     */
    const MultiIndex<D> limits_;
    
public:
    /**
     * An enumeration is by design immutable, therefore the same instance
     * can be shared across multiple threads.
     */
    class Enumeration
    {
    private:
        const HyperCubicShape<D> &shape_;
        
    public:
        class Iterator
        {
        private:
            const HyperCubicShape<D> &shape_;
            const MultiIndex<D> index_;
            
        public:
            Iterator(const HyperCubicShape<D> shape, MultiIndex<D> index) : shape_(shape), index_(index) {}
            Iterator(const Iterator &that) : shape_(that.shape_), index_(that.index_) {}
            
            MultiIndex<D> getMultiIndex() const
            {
                return index_;
            }
            
            std::size_t getOrdinal() const
            {
                std::size_t ordinal = 0;
                std::size_t area = 1;
                for (dim_t d = 0; d < D; d++) {
                    ordinal += area*index_[d];
                    area *= 1+shape_.limits_[d];
                }
                return ordinal;
            }
            
            bool getForwardNeighbour(dim_t axis, Iterator &it) const
            {
                if (index_[axis] == 0)
                    return false;
                
                MultiIndex<D> neighbour(index_);
                neighbour[axis] += 1;
                it = Iterator(shape_, neighbour);
                return true;
            }
            
            Iterator getBackwardNeighbour(dim_t axis) const
            {
                if (index_[axis] == 0)
                    return false;
                
                MultiIndex<D> neighbour(index_);
                neighbour[axis] -= 1;
                return Iterator(shape_, neighbour);
            }
        };
        
        Enumeration(const HyperCubicShape<D> &shape) : shape_(shape) {}
        Enumeration(const Enumeration &that) : shape_(that.shape_) {}
        
        /**
        * 
        */
        std::size_t size() const
        {
            std::size_t size = 1;
            for (dim_t d = 0; d < D; d++)
                size *= shape_.limits_[d];
            return size;
        }
        
        /**
        * get first multi-index
        */
        Iterator first() const
        {
            return Iterator(shape_, MultiIndex<D>());
        }
        
        /**
        * 
        */
        Iterator last() const
        {
            return Iterator(shape_, shape_.limits_);
        }
        
        /**
        * get multi-index of ordinal
        * behaviour is undefined if 'ordinal >= this->size()'
        * 
        * worst case complexity: amortized log(size)
        */
        Iterator at(std::size_t ordinal) const
        {
            assert(ordinal < size);
            
            MultiIndex<D> index;
            std::size_t area = 1;
            for (dim_t d = 0; d < D; d++) {
                index[d] = ordinal%area;
                ordinal -= ordinal%area;
                area *= 1+shape_.limits_[d];
            }
            return Iterator(shape_, index);
        }
        
        /**
        * check whether this enumeration contains a multi-index
        * 
        * worst case complexity: amortized log(size)
        */
        bool contains(const MultiIndex &index) const
        {
            for (dim_t d = 0; d < D; d++)
                if (index[d] > shape_.box_[d])
                    return false;
            return true;
        }
        
        /**
        * lookup ordinal of multi-index
        * return value is this->end() if 'this->contains(index) == false'
        * 
        * worst case complexity: amortized log(size)
        */
        bool find(const MultiIndex &index, Iterator &it) const
        {
            if (!contains(index))
                return false;
            
            return Iterator(shape_, index);
        }
    };
};

};

#endif