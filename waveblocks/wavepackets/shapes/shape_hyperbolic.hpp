#pragma once

#include <cmath>
#include <string>
#include <sstream>

#include "../../basic_types.hpp"

#include "shape_base.hpp"


namespace waveblocks {
    namespace wavepackets {
        namespace shapes {
            /**
             * \brief This class implements the hyperbolic cut shape.
             *
             * This class implements the hyperbolic cut basis shape, which is a special
             * type of a sparse basis shape.
             * The hyperbolic cut shape in \f$ D \f$ dimensions with _sparsity_ \f$S\f$
             * is defined as the set
             *
             * \f[
             * \mathfrak{K}(D,S) := \left\{(k_1,\dots,k_D) \in \mathbb{N}_0^D \mid
             *      \displaystyle\prod_{d=1}^{D} (1+k_d) \leq S \right\}
             * \f]
             *
             * \tparam D basis shape dimensionality
             */
            template<dim_t D>
            class HyperbolicCutShape : public AbstractShape<D>
            {
            private:
                int S_;

            public:
                /**
                 * \brief General constructor to set the sparsity parameter \f$ S \f$.
                 *
                 * \param S The sparsity parameter \f$ S \f$.
                 */
                HyperbolicCutShape(int S) : S_(S) {}

                virtual int bbox(dim_t axis) const override
                {
                    (void)(axis); // unused
                    return S_ - 1;
                }

                virtual int limit(int const* base_node, dim_t axis) const override
                {
                    double s = S_;

                    for (dim_t i = 0; i < D; i++) {
                        if (i != axis) {
                            s /= 1 + base_node[i];
                        }
                    }

                    return (int)s - 1;
                }

                virtual void print(std::ostream & out) const override
                {
                    out << "HyperbolicCutShape{sparsity: " << S_ << "}";
                }
            };

            /**
             * \brief This class implements the limited hyperbolic cut shape.
             *
             * This class implements the limited hyperbolic cut basis shape which is a special
             * type of a sparse basis shape.
             * The limited hyperbolic cut shape in \f$ D \f$ dimensions with _sparsity_ \f$S\f$ and
             * _limits_ \f$ \boldsymbol{K} = (K_1,\ldots,K_D) \f$ is defined as the set
             *
             * \f[
             * \mathfrak{K}(D,S,\boldsymbol{K}) := \left\{(k_1,\dots,k_D) \in \mathbb{N}_0^D \mid
             *      0 \leq k_d < K_d \; \land
             *      \displaystyle\prod_{d=1}^{D} (1+k_d) \leq S \right\}
             * \f]
             *
             * It is an intersection of the hyperbolic cut shape with a hypercubic shape.
             *
             * \tparam D basis shape dimensionality
             */
            template<dim_t D>
            class LimitedHyperbolicCutShape : public AbstractShape<D>
            {
            private:
                int S_;
                std::array<int,D> limits_;

            public:
                /**
                 * \brief General constructor to define sparsity parameter and limits.
                 *
                 * \param S The sparsity parameter \f$ S \f$.
                 * \param limits Tuple of all limits \f$ \boldsymbol{K} \f$.
                 */
                LimitedHyperbolicCutShape(int S, const std::array<int,D> &limits)
                    : S_(S)
                    , limits_(limits)
                { }

                /**
                 * \brief Specialized constructor to set all limits \f$ K_d \f$ to the same value \f$ K^\star \f$.
                 *
                 * \param S The sparsity parameter \f$ S \f$.
                 * \param size The limit \f$ K^\star \f$.
                 */
                LimitedHyperbolicCutShape(int S, int size)
                    : S_(S)
                {
                    for (std::size_t d = 0; d < D; d++)
                        limits_[d] = size;
                }

                /**
                 * \brief General constructor to define sparsity parameter and limits.
                 *
                 * \param S The sparsity parameter \f$ S \f$.
                 * \param list List of all limits \f$ \boldsymbol{K} \f$.
                 */
                LimitedHyperbolicCutShape(int S, std::initializer_list<int> list)
                    : S_(S)
                {
                    int deflt = 0;
                    std::size_t i = 0;
                    for (int e : list) {
                        limits_[i++] = deflt = e;
                    }
                    //fill remaining elements with last value of initializer list
                    while (i < D) {
                        limits_[i++] = deflt;
                    }
                }

                virtual int bbox(dim_t axis) const override
                {
                    return std::min( limits_[axis]-1, S_ - 1);
                }

                virtual int limit(int const* base_node, dim_t axis) const override
                {
                    double s = S_;

                    for (dim_t i = 0; i < D; i++) {
                        if (i != axis) {
                            if (base_node[i] >= limits_[i])
                                return -1;
                            else
                                s /= 1 + base_node[i];
                        }
                    }

                    return std::min((int)limits_[axis]-1, (int)s - 1);
                }

                virtual void print(std::ostream & out) const override
                {
                    out << "LimitedHyperbolicCutShape{ sparsity: " << S_ << ", limits (exclusive): [";
                    for (dim_t i = 0; i < D-1; i++) {
                        out << limits_[i] << ",";
                    }
                    out << limits_[D-1] << "]";
                    out << "}";
                }
            };
        }
    }
}
