#ifndef WAVEBLOCKS_POLYMORPHOUS_SHAPE
#define WAVEBLOCKS_POLYMORPHOUS_SHAPE

namespace waveblocks {

template<dim_t D>
class BasicShape
{
public:
    virtual int surface(dim_t axis, const MultiIndex<D> &coordinate) const = 0;
    virtual MultiIndex<D> getLimits() const = 0;
};

/*
template<dim_t D>
class PolymorphousShape
{
private:
    const BasicShape<D> *shape_;
    
public:
    PolymorphousShape(const BasicShape *shape) : shape_(shape) {}
    PolymorphousShape(const PolymorphousShape &that) : shape_(that.shape_) {}
    
    PolymorphousShape &operator=(const PolymorphousShape &that)
    {
        shape_ = that.shape_;
        return *this;
    }
    
    int surface(dim_t axis, const MultiIndex<D> &coordinate) const
    {
        return shape_->surface(axis, coordinate);
    }
    
    MultiIndex<D> getLimits() const
    {
        return shape_->getLimits();
    }
};*/

}

#endif