#pragma once

#include <iostream>


namespace waveblocks {
    namespace wavepackets {
        namespace shapes {
            /**
             * \brief Subclasses provide a description of a basis shape.
             *
             * A \f$ D \f$-dimensional basis shape \f$ \mathfrak{K} \f$ is a set
             * of _unordered_ \f$ D \f$ dimensional integer tuples (aka _node_).
             *
             * Subclasses provide an description of a basis shape
             * \f$ \mathfrak{K} \subset \mathbb{N}_0^D\f$.
             * It describes, which nodes \f$ \underline{k} \in \mathbb{N}_0^D \f$
             * are part of the shape.
             * Instances are passed to the ShapeEnumerator,
             * which converts the information to an ShapeEnum.
             *
             * Keep in mind, that basis shapes must fulfill the fundamental property
             * \f[
             * \underline{k} \in \mathfrak{K} \Rightarrow \forall
             * \underline{k}-\underline{e}^d \in \mathfrak{K} \;\forall d \in \{d \;|\;k_d \geq 1\}
             * \f]
             * where \f$ \underline{e}^d \f$ is the unit vector in direction \f$ d \f$.
             * That means, if an arbitrary node is part of the basis shape, then all nodes
             * in the backward cone are part of the shape too.
             *
             * \tparam D basis shape dimensionality
             */
            template<dim_t D>
            class AbstractShape
            {
            public:
                virtual ~AbstractShape() { }

                /**
                 * \brief Retrieves the length of the minimum bounding box in one direction.
                 *
                 * The minimum bounding box is given by
                 * \f[
                 * L_{\alpha}=\max_{k_{\alpha}}\left\{\underline{k} \in \mathfrak{K}\right\}
                 * \f]
                 *
                 * \param axis The direction \f$ \alpha \f$.
                 * \return Length of the bbox.
                 */
                virtual int bbox(dim_t axis) const = 0;

                /**
                 * \brief Evaluates one surface function on a base node.
                 *
                 * The surface function to direction \f$ \alpha \f$ is given by
                 *
                 * \f[
                 * s_{\alpha}(\underline{n})=\max_{k_{\alpha}}
                 * \left\{\underline{k} \in \mathfrak{K} \;|\;
                 * k_d = n_d \; \forall d \neq \alpha
                 * \right\}
                 * \f]
                 *
                 * Notice that the \f$ \alpha \f$-th entry of \f$ \underline{n} \f$
                 * does not influence return value.
                 * It can be of any value since it is simply ignored.
                 *
                 * \param base_node The basis node \f$ \underline{n} \f$. It contains D indices.
                 * \param axis The direction \f$ \alpha \f$.
                 * \return Value of the surface function.
                 */
                virtual int limit(int const* base_node, dim_t axis) const = 0;

                /**
                 * \brief Prints a pretty description of the shape.
                 *
                 * \param out The output stream.
                 */
                virtual void print(std::ostream & out) const = 0;
            };

            template<dim_t D>
            std::ostream & operator<<(std::ostream & out, AbstractShape<D> const& shape)
            {
                shape.print(out);
                return out;
            }
        }
    }
}
