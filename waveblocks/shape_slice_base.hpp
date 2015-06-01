#ifndef WAVEBLOCKS_SHAPE_SLICE
#define WAVEBLOCKS_SHAPE_SLICE

template<dim_t D>
class ShapeSlice
{
public:
    virtual ~Slice() { }
    
    /**
     * \return number of nodes in all previous slices
     */
    virtual std::size_t offset() const = 0;
    
    /**
     * \return number of nodes in this slice
     */
    virtual std::size_t size() const = 0;
    
    /**
     * \brief Returns the multi-index of the node at position <i>ordinal</i>.
     * 
     * Notice that the first node in the slice has ordinal 0 (not 1 or offset()).
     * 
     * Portable programs should never call this function with an argument that is <i>out-of-range</i>,
     * since this causes <i>undefined behaviour</i>.
     * 
     * <b>complexity: </b>logarithmic in the number of slice-nodes
     * 
     * \param[in] ordinal position of a node in this slice
     * \return multi-index of the specified node
     */
    virtual std::array<int,D> operator[](std::size_t ordinal) const = 0;
    
    /**
     * \brief Returns the position of the node with multi-index <i>index</i>.
     * 
     * Notice that the first node in the slice has position 0 (not 1 or offset()).
     * 
     * Portable programs should never call this function with an argument 
     * that is inexistant in this slice since this causes <i>undefined behaviour</i>.
     * 
     * Use contains(index) to check whether this slice contains this node.
     * 
     * <b>complexity: </b>logarithmic in the number of slice-nodes
     * 
     * \param[in] index multi-index of a node in this slice
     * \return position of the specified node
     */
    virtual std::size_t find(const std::array<int,D> &index) const = 0;
    
    /**
     * \brief const_iterator over a slice to support foreach statements
     */
    class Iterator
    {
    private:
        const Slice *ref_;
        std::size_t ientry_;
        
    public:
        Iterator(const Slice *ref, std::size_t ientry)
        : ref_(ref)
        , ientry_(ientry)
        { }
        
        Iterator(const Iterator &other)
            : ref_(other.ref_)
            , ientry_(other.ientry_)
        { }
        
        Iterator &operator=(const Iterator &other)
        {
            ref_ = other.ref_;
            ientry_ = other.ientry_;
            return *this;
        }
        
        Iterator &operator++()
        {
            ++ientry_;
            return *this;
        }
        
        std::array<int,D> operator*() const
        {
            return (*ref_)[ientry_];
        }
        
        bool operator==(const Iterator &other) const
        {
            return ientry_ == other.ientry_ && ref_ == other.ref_;
        }
        
        bool operator!=(const Iterator &other) const
        {
            return !operator==(other);
        }
    };
    
    Iterator begin() const
    {
        return Iterator(this, 0);
    }
    
    Iterator end() const
    {
        return Iterator(this, size());
    }
};

#endif