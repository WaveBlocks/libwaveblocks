#pragma once

#include <stdexcept>
#include <memory>

#include "hawp_paramset.hpp"
#include "hawp_evaluator.hpp"
#include "shapes/shape_extension_cache.hpp"


namespace waveblocks {
    namespace wavepackets {
        /**
         * \brief Abstract superclass that represents a set
         * of basis function to a scalar Hagedorn wavepacket.
         *
         * It provides read-only access to
         * the scaling parameter \f$ \varepsilon \f$,
         * the Hagedorn parameter set \f$ \Pi \f$ and
         * the basis shape \f$ \mathfrak{K} \f$.
         *
         * Therefore it is able to evaluate its basis functions \f$ \phi_k \f$ on quadrature points \f$ x \f$.
         */
        template<dim_t D, class MultiIndex>
        class AbstractScalarHaWpBasis
        {
        public:
            /**
             * \brief Retrieves the semi-classical scaling parameter
             * \f$ \varepsilon \f$ of the wavepacket.
             */
            virtual real_t eps() const = 0;

            /**
             * \brief Grants read-only access to the Hagedorn parameter set
             * \f$ \Pi \f$ of the wavepacket.
             */
            virtual HaWpParamSet<D> const& parameters() const = 0;

            /**
             * \brief Retrieves the basis shape \f$ \mathfrak{K} \f$
             * of the wavepacket.
             */
            virtual shapes::ShapeEnumSharedPtr<D, MultiIndex> shape() const = 0;

            template<int N> HaWpEvaluator<D,MultiIndex,N>
            create_evaluator(CMatrix<D,N> const& grid) const
            {
                return {eps(), &parameters(), shape().get(), grid};
            }

            /**
             * \brief Evaluates all basis functions \f$ \{\phi_k\} \f$ on complex grid nodes \f$ x \in \gamma \f$.
             *
             * \param grid
             * Complex grid nodes / quadrature points \f$ \gamma \f$.
             * Complex matrix with shape (dimensionality, number of grid nodes).
             * \return Complex 2D-array with shape (basis shape size, number of grid nodes)
             */
            template<int N> HaWpBasisVector<N>
            evaluate_basis(CMatrix<D,N> const& grid) const
            {
                return create_evaluator(grid).all();
            }

            /**
             * \brief Evaluates all basis functions \f$ \{\phi_k\} \f$ on real grid nodes \f$ x \in \gamma \f$.
             *
             * \param rgrid
             * Real grid nodes / quadrature points \f$ \gamma \f$.
             * Real matrix with shape (dimensionality, number of grid nodes).
             * \return Complex 2D-array with shape (basis shape size, number of grid nodes)
             *
             * \tparam N
             * Number of quadrature points.
             * Don't choose Eigen::Dynamic. It works, but performance is bad.
             */
            template<int N> HaWpBasisVector<N>
            evaluate_basis(RMatrix<D,N> const& rgrid) const
            {
                CMatrix<D,N> cgrid = rgrid.template cast <complex_t>();
                return evaluate_basis(cgrid);
            }

            //     virtual HaWpBasisVector<Eigen::Dynamic> evaluate_basis(ComplexGrid<D,Eigen::Dynamic> const& grid) const
            //     {
            //         return create_evaluator<Eigen::Dynamic>(grid).all();
            //     }

            /**
             * \brief Computes the extension \f$ \mathfrak{K}_{ext} \f$ of the stored
             * basis shape \f$ \mathfrak{K} \f$.
             *
             *
             * Computing an extended shape is expensive. Therefore this function
             * caches computed extensions.
             *
             * \e Thread-Safety: The stored pointer to the cached shape extension is not guarded by a mutex.
             * Therefore race conditions may occur when calling this function concurrently.
             *
             * \return Shared pointer to the extended shape.
             */
            shapes::ShapeEnumSharedPtr<D,MultiIndex>
            extended_shape() const
            {
                return shape_extension_cache_.get_extended_shape( this->shape() );
            }

        private:
            shapes::ShapeExtensionCache<D,MultiIndex> shape_extension_cache_;
        };


        /**
         * \brief Abstract superclass that represents a scalar (1-component) hagedorn wavepacket.
         *
         * A subclass provides read-only access to
         * the scaling parameter \f$ \varepsilon \f$,
         * the Hagedorn parameter set \f$ \Pi \f$,
         * the basis shape \f$ \mathfrak{K} \f$ and
         * the coefficients \f$ c \f$.
         *
         * Therefore it is able to evaluate itself (\f$ \Phi(x) \f$) on grid points \f$ x \f$.
         *
         * \tparam D wavepacket dimensionality
         * \tparam MultiIndex type to represent a multi-index
         */
        template<dim_t D, class MultiIndex>
        class AbstractScalarHaWp : public AbstractScalarHaWpBasis<D,MultiIndex>
        {
        public:
            virtual real_t eps() const = 0;
            virtual HaWpParamSet<D> const& parameters() const = 0;
            virtual shapes::ShapeEnumSharedPtr<D, MultiIndex> shape() const = 0;

            /**
             * \brief Grants read-only access to the coefficients \f$ \{c_k\} \f$
             * for all \f$ k \in \mathfrak{K} \f$ of this wavepacket.
             */
            virtual Coefficients const& coefficients() const = 0;

            /**
             * \brief Evaluates this wavepacket \f$ \Phi(x) \f$ at complex grid nodes \f$ x \in \gamma \f$.
             *
             * Notice that this function does not include the prefactor
             * \f$ \frac{1}{\sqrt{det(Q)}} \f$ nor the global phase
             * \f$ \exp{(\frac{iS}{\varepsilon^2})} \f$.
             *
             * \param grid
             * Complex grid nodes / quadrature points \f$ \gamma \f$.
             * Complex matrix with shape (dimensionality, number of grid nodes).
             * \return Complex matrix with shape (1, number of grid nodes)
             */
            template<int N> CArray<1,N>
            evaluate(CMatrix<D,N> const& grid) const
            {
                if (this->shape()->n_entries() != (std::size_t)coefficients().size())
                    throw std::runtime_error("shape.size() != coefficients.size()");

                return this->template create_evaluator<N>(grid).reduce(coefficients());
            }

            /**
             * \brief Evaluates this wavepacket \f$ \Phi(x) \f$ at real grid nodes \f$ x \in \gamma \f$.
             *
             * Notice that this function does not include the prefactor
             * \f$ \frac{1}{\sqrt{det(Q)}} \f$ nor the global phase
             * \f$ \exp{(\frac{iS}{\varepsilon^2})} \f$.
             *
             * \param rgrid
             * Real grid nodes / quadrature points \f$ \gamma \f$.
             * Real matrix with shape (dimensionality, number of grid nodes).
             * \return Complex matrix with shape (1, number of grid nodes)
             */
            template<int N> CArray<1,N>
            evaluate(RMatrix<D,N> const& rgrid) const
            {
                CMatrix<D,N> cgrid = rgrid.template cast <complex_t>();
                return evaluate(cgrid);
            }

            //     virtual HaWpBasisVector<Eigen::Dynamic> evaluate(ComplexGrid<D,Eigen::Dynamic> const& grid) const
            //     {
            //         return create_evaluator<Eigen::Dynamic>().reduce(coefficients());
            //     }

            /**
             * \brief Computes the prefactor \f$ \frac{1}{\sqrt{det(Q)}} \f$.
             */
            complex_t prefactor() const
            {
                return real_t(1) / this->parameters().sdQ();
            }

            /**
             * \brief Computes the global phase factor \f$ \exp{(\frac{i S}{\varepsilon^2})} \f$.
             */
            complex_t phasefactor() const
            {
                return std::exp(complex_t(0,1) * this->parameters().S() / eps() / eps());
            }
        };


        /**
         * \brief Concrete implementation of a scalar Hagedorn wavepacket.
         */
        template<dim_t D, class MultiIndex>
        class ScalarHaWp : public AbstractScalarHaWp<D, MultiIndex>
        {
        public:
            /**
             * \brief Grants writeable access to the semi-classical scaling parameter
             * \f$ \varepsilon \f$ of the wavepacket.
             */
            real_t & eps()
            {
                return eps_;
            }

            real_t eps() const override
            {
                return eps_;
            }

            /**
             * \brief Grants writeable access to the Hagedorn parameter set
             * \f$ \Pi \f$ of the wavepacket.
             */
            HaWpParamSet<D> & parameters()
            {
                return parameters_;
            }

            HaWpParamSet<D> const& parameters() const override
            {
                return parameters_;
            }

            /**
             * \brief Grants access to the basis shape
             * \f$ \mathfrak{K} \f$ of the wavepacket.
             *
             * \return
             * Reference to the shape enumeration pointer.
             * You can assign a new pointer to it!
             */
            shapes::ShapeEnumSharedPtr<D, MultiIndex> & shape()
            {
                return shape_;
            }

            shapes::ShapeEnumSharedPtr<D, MultiIndex> shape() const override
            {
                return shape_;
            }

            /**
             * \brief Grants writeable access to the coefficients \f$ c \f$
             * of the wavepacket.
             */
            Coefficients & coefficients()
            {
                return coefficients_;
            }

            Coefficients const& coefficients() const override
            {
                return coefficients_;
            }

        private:
            real_t eps_;
            HaWpParamSet<D> parameters_;
            shapes::ShapeEnumSharedPtr<D, MultiIndex> shape_;
            Coefficients coefficients_;
        };


        /**
         * \brief Represents a homogeneous Hagedorn wavepacket \f$ \Psi \f$
         * with \f$ N \f$ components \f$ \Phi_n \f$.
         * All components share the same Hagedorn parameterset \f$ \Pi \f$
         * and scaling parameter \f$ \varepsilon \f$.
         *
         * The number of components is determined at runtime.
         *
         * \tparam D wavepacket dimensionality
         * \tparam MultiIndex
         */
        template<dim_t D, class MultiIndex>
        class HomogeneousHaWp
        {
        public:
            /**
             * \brief Represents a component of a homogeneous wavepacket.
             *
             * Such a component is a full-fledged scalar Hagedorn wavepacket
             * that shares the Hagedorn parameter \f$ \Pi \f$ and scaling
             * parameter \f$ \varepsilon \f$ with other components.
             */
            class Component : public AbstractScalarHaWp<D,MultiIndex>
            {
            public:
                Component(HomogeneousHaWp const* const owner)
                    : owner_(owner)
                    , coefficients_()
                { }

                Component(Component&& that)
                    : owner_(that.owner_)
                    , shape_(that.shape_)
                    , coefficients_(std::move(that.coefficients_))
                { }

                Component(Component const& that)
                    : owner_(that.owner_)
                    , shape_(that.shape_)
                    , coefficients_(that.coefficients_)
                { }

                Component & operator=(Component&& that)
                {
                    shape_ = that.shape_;
                    coefficients_ = std::move(that.coefficients_);
                    return *this;
                }

                Component & operator=(Component const& that)
                {
                    shape_ = that.shape_;
                    coefficients_ = that.coefficients_;
                    return *this;
                }

                real_t eps() const override
                {
                    return owner_->eps();
                }

                HaWpParamSet<D> const& parameters() const override
                {
                    return owner_->parameters();
                }

                /**
                 * \brief Grants access to the basis shape
                 * \f$ \mathfrak{K} \f$ of the wavepacket.
                 *
                 * \return
                 * Reference to the shape enumeration pointer.
                 * You can assign a new pointer to it!
                 */
                shapes::ShapeEnumSharedPtr<D, MultiIndex> & shape()
                {
                    return shape_;
                }

                shapes::ShapeEnumSharedPtr<D, MultiIndex> shape() const override
                {
                    return shape_;
                }

                /**
                 * \brief Grants writeable access to the coefficients of this
                 * wavepacket component.
                 */
                Coefficients & coefficients()
                {
                    return coefficients_;
                }

                Coefficients const& coefficients() const override
                {
                    return coefficients_;
                }

            private:
                HomogeneousHaWp const* const owner_;

                shapes::ShapeEnumSharedPtr<D, MultiIndex> shape_;
                Coefficients coefficients_;
            };

            HomogeneousHaWp(std::size_t n)
                : eps_(0.0)
                , parameters_()
                , components_(n, Component(this))
            { }

            /**
             * \brief Grants access to the semi-classical scaling parameter
             * \f$ \varepsilon \f$ of the wavepacket.
             */
            real_t & eps()
            {
                return eps_;
            }

            /**
             * \brief Retrieves the semi-classical scaling parameter
             * \f$ \varepsilon \f$ of the wavepacket.
             */
            real_t eps() const
            {
                return eps_;
            }

            /**
             * \brief Grants writeable access to the Hagedorn parameter set
             * \f$ \Pi \f$ of the wavepacket.
             */
            HaWpParamSet<D> & parameters()
            {
                return parameters_;
            }

            /**
             * \brief Grants read-only access to the Hagedorn parameter set
             * \f$ \Pi \f$ of the wavepacket.
             */
            HaWpParamSet<D> const& parameters() const
            {
                return parameters_;
            }

            /**
             * \brief Grants writeable access to all components
             * \f$ \{\Phi_n\} \f$ of this wavepacket.
             */
            std::vector<Component> & components()
            {
                return components_;
            }

            /**
             * \brief Grants read-only access to all components
             * \f$ \{\Phi_n\} \f$ of this wavepacket.
             */
            std::vector<Component> const& components() const
            {
                return components_;
            }

            /**
             * \brief Grants writeable access to the \f$ n \f$-th component
             * \f$ \Phi_n \f$.
             *
             * \param n The index \f$ n \f$ of the requested component.
             * \return Reference to the requested component.
             */
            Component & component(std::size_t n)
            {
                return components_[n];
            }

            /**
             * \brief Grants read-only access to the \f$ n \f$-th component
             * \f$ \Phi_n \f$.
             *
             * \param n The index \f$ n \f$ of the requested component.
             * \return Reference to the requested component.
             */
            Component const& component(std::size_t n) const
            {
                return components_[n];
            }

            /**
             * \brief Grants writeable access to the \f$ n \f$-th component
             * \f$ \Phi_n \f$.
             *
             * \param n The index \f$ n \f$ of the requested component.
             * \return Reference to the requested component.
             */
            Component & operator[](std::size_t n)
            {
                return component(n);
            }

            /**
             * \brief Grants read-only access to the \f$ n \f$-th component
             * \f$ \Phi_n \f$.
             *
             * \param n The index \f$ n \f$ of the requested component.
             * \return Reference to the requested component.
             */
            Component const& operator[](std::size_t n) const
            {
                return component(n);
            }

            /**
             * \brief Returns the number of components.
             */
            std::size_t n_components() const
            {
                return components_.size();
            }

            /**
             * \brief Computes the union of basis shapes of all components.
             *
             * _Thread Safety:_ This function caches the basis shape union.
             * The cache is not protected by a mutex. This concurrent access may
             * introduce race conditions.
             */
            shapes::ShapeEnumSharedPtr<D,MultiIndex> union_shape() const
            {
                // check cache status
                bool rebuild_cache = false;

                if (union_cache_snapshot_.size() != n_components()) {
                    rebuild_cache = true;
                    union_cache_snapshot_.resize(n_components());
                }

                for (std::size_t n = 0; n < n_components() && !rebuild_cache; n++) {
                    if (union_cache_snapshot_[n] != component(n).shape().get())
                        rebuild_cache = true;
                }

                // rebuild cache if necessary
                if (rebuild_cache) {
                    std::vector< ShapeEnum<D,MultiIndex> const* > list(n_components());
                    for (std::size_t c = 0; c < n_components(); c++) {
                        list[c] = union_cache_snapshot_[c] = component(c).shape().get();
                    }
                    cached_shape_union_ = std::make_shared< ShapeEnum<D,MultiIndex> >(shapes::shape_enum::strict_union(list));
                }

                return cached_shape_union_;
            }

            /**
             * \brief Evaluate the value of all components at once.
             *
             * Evaluates \f$ \Psi(x) = \{\Phi_i(x)\} \f$,
             * where \f$ x \f$ is is a complex quadrature point.
             *
             * _Thread Safety:_ Some intermediate results are cached. The cache is not
             * protected by a mutex. Thus concurrent access may introduce race conditions.
             *
             * \param grid
             * Complex quadrature points.
             * Complex matrix of shape (dimensionality, number of quadrature points)
             * \return
             * Complex matrix of shape (number of components, number of quadrature points)
             *
             * \tparam N
             * Number of quadrature points.
             * Don't use Eigen::Dynamic. It works, but performance is bad.
             */
            template<int N>
            CArray<Eigen::Dynamic,N> evaluate(CMatrix<D,N> const& grid) const
            {
                ScalarHaWp<D,MultiIndex> unionwp;

                unionwp.eps() = eps();
                unionwp.parameters() = parameters();
                unionwp.shape() = union_shape();

                std::vector< ShapeEnum<D,MultiIndex>* > shape_list(n_components());
                std::vector< complex_t const* > coeffs_list(n_components());

                for (std::size_t n = 0; n < n_components(); n++) {
                    shape_list[n] = component(n).shape().get();
                    coeffs_list[n] = component(n).coefficients().data();
                }

                return unionwp.template create_evaluator<N>(grid).vector_reduce(shape_list.data(), coeffs_list.data(), n_components());
            }

            /**
             * \brief Evaluates the value of all components at once.
             *
             * Evaluates \f$ \Psi(x) = \{\Phi_i(x)\} \f$,
             * where \f$ x \f$ is is a real quadrature point.
             *
             * \param rgrid
             * Real quadrature points.
             * Real matrix of shape (dimensionality, number of quadrature points)
             * \return
             * Complex matrix of shape (number of components, number of quadrature points)
             *
             * \tparam N
             * Number of quadrature points.
             * Don't use Eigen::Dynamic. It works, but performance is bad.
             */
            template<int N>
            CArray<Eigen::Dynamic,N> evaluate(RMatrix<D,N> const& rgrid) const
            {
                CMatrix<D,N> cgrid = rgrid.template cast<complex_t>();
                return evaluate(cgrid);
            }

        private:
            real_t eps_;
            HaWpParamSet<D> parameters_;
            std::vector<Component> components_;

            mutable std::vector< ShapeEnum<D,MultiIndex>* > union_cache_snapshot_;
            mutable shapes::ShapeEnumSharedPtr<D,MultiIndex> cached_shape_union_;
        }; // class HomogeneousHaWp


        /**
         * \brief Represents an inhomogeneous Hagedorn wavepacket \f$ \Psi \f$
         * with \f$ N \f$ components \f$ \Phi_n \f$.
         * All components have a different set of Hagedorn parameters \f$ \Pi \f$,
         * basis shapes \f$ \mathfrak{K} \f$ and coefficients \f$ c \f$.
         *
         * The number of components is determined at runtime.
         *
         * \tparam D wavepacket dimensionality
         * \tparam MultiIndex
         */
        template<dim_t D, class MultiIndex>
        class InhomogeneousHaWp
        {
        public:
            /**
             * \brief Represents a component \f$ \Phi_n \f$ of an
             * inhomogeneous wavepacket \f$ \Psi \f$.
             *
             * Such a component is a full-fledged scalar Hagedorn wavepacket
             * that shares only the scaling parameter \f$ \varepsilon \f$ with other
             * components.
             */
            class Component : public AbstractScalarHaWp<D,MultiIndex>
            {
            public:
                Component(InhomogeneousHaWp const* const owner)
                    : owner_(owner)
                    , shape_()
                    , coefficients_()
                { }

                Component(Component&& that)
                    : owner_(that.owner_)
                    , parameters_(std::move(that.parameters_))
                    , shape_(that.shape_)
                    , coefficients_(std::move(that.coefficients_))
                { }

                Component(Component const& that)
                    : owner_(that.owner_)
                    , parameters_(that.parameters_)
                    , shape_(that.shape_)
                    , coefficients_(that.coefficients_)
                { }

                Component & operator=(Component&& that)
                {
                    parameters_ = std::move(that.parameters_);
                    shape_ = that.shape_;
                    coefficients_ = std::move(that.coefficients_);
                    return *this;
                }

                Component & operator=(Component const& that)
                {
                    parameters_ = that.parameters_;
                    shape_ = that.shape_;
                    coefficients_ = that.coefficients_;
                    return *this;
                }

                /**
                 * \brief Forwards the scaling parameter \f$ \varepsilon \f$
                 * of the owning inhomogeneous wavepacket.
                 */
                real_t eps() const override
                {
                    return owner_->eps();
                }

                /**
                 * \brief Grants writeable access to the Hagedorn parameter set \f$ \Pi \f$.
                 */
                HaWpParamSet<D> & parameters()
                {
                    return parameters_;
                }


                HaWpParamSet<D> const& parameters() const override
                {
                    return parameters_;
                }

                /**
                 * \brief Grants access to the basis shape
                 * \f$ \mathfrak{K} \f$ of the wavepacket.
                 *
                 * \return
                 * Reference to the shape enumeration pointer.
                 * You can assign a new pointer to it!
                 */
                shapes::ShapeEnumSharedPtr<D, MultiIndex> & shape()
                {
                    return shape_;
                }

                shapes::ShapeEnumSharedPtr<D, MultiIndex> shape() const override
                {
                    return shape_;
                }

                /**
                 * \brief Grants writeable access to the coefficients of this
                 * wavepacket component.
                 */
                Coefficients & coefficients()
                {
                    return coefficients_;
                }

                Coefficients const& coefficients() const override
                {
                    return coefficients_;
                }

            private:
                InhomogeneousHaWp const* const owner_;

                HaWpParamSet<D> parameters_;
                shapes::ShapeEnumSharedPtr<D, MultiIndex> shape_;
                Coefficients coefficients_;
            };

            InhomogeneousHaWp(std::size_t n)
                : eps_()
                , components_(n, Component(this))
            { }

            /**
             * \brief Grants access to the semi-classical scaling parameter
             * \f$ \varepsilon \f$ of the wavepacket.
             */
            real_t & eps()
            {
                return eps_;
            }

            /**
             * \brief Retrieves the semi-classical scaling parameter
             * \f$ \varepsilon \f$ of the wavepacket.
             */
            real_t eps() const
            {
                return eps_;
            }

            /**
             * \brief Grants writeable access to all components
             * \f$ \{\Phi_n\} \f$ of this wavepacket.
             */
            std::vector<Component> & components()
            {
                return components_;
            }

            /**
             * \brief Grants read-only access to all components
             * \f$ \{\Phi_n\} \f$ of this wavepacket.
             */
            std::vector<Component> const& components() const
            {
                return components_;
            }

            /**
             * \brief Grants writeable access to the \f$ n \f$-th component
             * \f$ \Phi_n \f$.
             *
             * \param n The index \f$ n \f$ of the requested component.
             * \return Reference to the requested component.
             */
            Component & component(std::size_t n)
            {
                return components_[n];
            }

            /**
             * \brief Grants read-only access to the \f$ n \f$-th component
             * \f$ \Phi_n \f$.
             *
             * \param n The index \f$ n \f$ of the requested component.
             * \return Reference to the requested component.
             */
            Component const& component(std::size_t n) const
            {
                return components_[n];
            }

            /**
             * \brief Grants writeable access to the \f$ n \f$-th component
             * \f$ \Phi_n \f$.
             *
             * \param n The index \f$ n \f$ of the requested component.
             * \return Reference to the requested component.
             */
            Component & operator[](std::size_t n)
            {
                return component(n);
            }

            /**
             * \brief Grants read-only access to the \f$ n \f$-th component
             * \f$ \Phi_n \f$.
             *
             * \param n The index \f$ n \f$ of the requested component.
             * \return Reference to the requested component.
             */
            Component const& operator[](std::size_t n) const
            {
                return component(n);
            }

            /**
             * \brief Returns the number of components.
             */
            std::size_t n_components() const
            {
                return components_.size();
            }

            /**
             * \brief Evaluate the value of all components at once.
             *
             * Evaluates \f$ \Psi(x) = \{\Phi_i(x)\} \f$,
             * where \f$ x \f$ is is a complex quadrature point.
             *
             * \param grid
             * Complex quadrature points.
             * Complex matrix of shape (dimensionality, number of quadrature points)
             * \return
             * Complex matrix of shape (number of components, number of quadrature points)
             *
             * \tparam N
             * Number of quadrature points.
             * Don't use Eigen::Dynamic. It works, but performance is bad.
             */
            template<int N>
            CArray<Eigen::Dynamic,N> evaluate(CMatrix<D,N> const& grid) const
            {
                CArray<Eigen::Dynamic,N> result(n_components(),grid.cols());

                for (std::size_t c = 0; c < n_components(); c++) {
                    result.row(c) = component(c).evaluate(grid);
                }

                return result;
            }

            /**
             * \brief Evaluates the value of all components at once.
             *
             * Evaluates \f$ \Psi(x) = \{\Phi_i(x)\} \f$,
             * where \f$ x \f$ is is a real quadrature point.
             *
             * \param rgrid
             * Real quadrature points.
             * Real matrix of shape (dimensionality, number of quadrature points)
             * \return
             * Complex matrix of shape (number of components, number of quadrature points)
             *
             * \tparam N
             * Number of quadrature points.
             * Don't use Eigen::Dynamic. It works, but performance is bad.
             */
            template<int N>
            CArray<Eigen::Dynamic,N> evaluate(RMatrix<D,N> const& rgrid) const
            {
                CMatrix<D,N> cgrid = rgrid.template cast<complex_t>();
                return evaluate(cgrid);
            }

        private:
            real_t eps_;
            std::vector<Component> components_;
        };
    }
}
