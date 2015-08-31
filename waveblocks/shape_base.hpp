#ifndef WAVEBLOCKS_ABSTRACT_SHAPE_HPP
#define WAVEBLOCKS_ABSTRACT_SHAPE_HPP

#include <iostream>

namespace waveblocks
{
    /**
     * \brief Abstract class that represents a shape object.
     */
    template<dim_t D>
    class AbstractShape
    {
    public:
        virtual ~AbstractShape() { }
        
        /**
         * \brief Returns upper shape limit in one direction.
         * 
         * Returns for a given \e axis \f$ j \f$ a (preferably) as small as possible \e limit \f$ L_j \f$
         * such that:
         * \f[
         * \forall \boldsymbol{k} \in \mathfrak{K} \,\colon\; k_j \leq L_j
         * \f]
         * \param axis \f$ j \f$
         * \return \f$ L_j \f$
         */
        virtual int bbox(dim_t axis) const = 0;
        
        /**
         * \brief Returns shape limit in one direction starting from a base node.
         * 
         * Returns for a given \e axis \f$ j \f$ and a given \e base \e node \f$ \boldsymbol{n} \f$
         * the largest element \f$ k^\star \f$ that satisfies: 
         * \f[ \boldsymbol{k} \in \mathfrak{K}, \;
         * k_i =
         *    \begin{cases}
         *       n_i,& i \neq j\\
         *       k^\star, & i = j
         *    \end{cases}
         * \f]
         * 
         * A \e base \e node is a node whose \f$ j \f$ -th entry is zero.
         * 
         * Notice that \f$ n_j \f$ does not influence return value \f$ k^\star \f$. 
         * Therefore \f$ n_j \f$ can be any value since it is simply ignored.
         * 
         * \param base_node base point \f$ \boldsymbol{n} \f$. \f$ n_j \f$ can be arbitrary.
         * \param axis \f$ j \f$
         * \return \f$ k^\star \f$
         */
        virtual int limit(int const* base_node, dim_t axis) const = 0;
        
        virtual void print(std::ostream & out) const = 0;
    };
    
    template<dim_t D>
    std::ostream & operator<<(std::ostream & out, AbstractShape<D> const& shape)
    {
        shape.print(out);
        return out;
    }
}



#endif