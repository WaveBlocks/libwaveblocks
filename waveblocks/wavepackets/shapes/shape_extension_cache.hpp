#pragma once

#include <memory>

#include "shape_enum_union.hpp"
#include "shape_enum_extended.hpp"


namespace waveblocks {
    namespace wavepackets {
        namespace shapes {

            template<dim_t D, class MultiIndex>
            using ShapeEnumSharedPtr = std::shared_ptr< ShapeEnum<D, MultiIndex> >;

            // template<dim_t D>
            // class AbstractScalarWavepacket
            // {
            // public:
            //     virtual Eigen::Array<complex_t,1,Eigen::Dynamic> evaluate(ComplexGrid<D,Eigen::Dynamic> const& grid) const = 0;
            //     virtual HaWpBasisVector<Eigen::Dynamic> evaluate_basis(ComplexGrid<D,Eigen::Dynamic> const& grid) const = 0;
            // };


            template<dim_t D, class MultiIndex>
            class ShapeExtensionCache
            {
            public:
                /**
                 * \brief shape reference to shape to check actuality of cache
                 * \return
                 */
                ShapeEnumSharedPtr<D,MultiIndex> get_extended_shape(std::shared_ptr< ShapeEnum<D,MultiIndex> > shape) const
                {
                    if (shape.get() != cached_extended_shape_source_)
                        update_extended_shape(shape); // this function is not thread-safe

                    return cached_extended_shape_;
                }

                /**
                 * \brief Manually sets extended shape
                 *
                 * \param shape source of new extended shape
                 * \param extension new extended shape
                 */
                void set_extended_shape(std::shared_ptr< ShapeEnum<D,MultiIndex> > shape,
                                        std::shared_ptr< ShapeEnum<D,MultiIndex> > extension)
                {
                    cached_extended_shape_source_ = shape.get();
                    cached_extended_shape_ = extension;
                }

                /**
                 * \brief Recomputes extended shape if source shape changed.
                 *
                 * \param shape new source shape
                 */
                void update_extended_shape(std::shared_ptr< ShapeEnum<D,MultiIndex> > shape) const
                {
                    if (shape.get() != cached_extended_shape_source_) {
                        cached_extended_shape_source_ = shape.get();
                        cached_extended_shape_ = std::make_shared< ShapeEnum<D,MultiIndex> >(shapes::shape_enum::extend(shape.get()));
                    }
                }

            private:
                mutable ShapeEnumSharedPtr<D,MultiIndex> cached_extended_shape_;
                mutable ShapeEnum<D,MultiIndex> * cached_extended_shape_source_; // source of the cached extended shape
            };
        }
    }
}
