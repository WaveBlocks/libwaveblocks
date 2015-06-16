
namespace waveblocks {

/**
 * \defgroup ShapeEnum Shape Enumeration
 * \ingroup WaveBlocksND
 * 
 * \brief Provides methods to enumerate basic shapes and manipulate enumerations.
 * 
 * ## Definition ##
 * A \f$ D \f$-dimensional shape enumeration \f$ \mathfrak{K} \f$ is a 
 * set of _ordered_ D-dimensional integer-tuples (aka _node_).
 * 
 * ## Rationale ##
 * A \ref BasicShape "shape" just tells you whether it contains a specific node. But
 * for most algorithms, you need to associate coefficients with shape nodes. One way to
 * to that is using a dictionary. But it is simpler to enumerate all nodes in a shape.
 * This way you can keep those coefficients in an array, ordered according to the enumeration. 
 * 
 * ## Implementation ##
 * Many algorithms, notable evaluation of a _hagedorn wavepacket_, use recursive formulas in the form 
 * \f$
 * c_{\underline{k}} = f(c_{\underline{k}-\underline{e}^1}, \ldots, c_{\underline{k}-\underline{e}^D}) 
 * \f$
 * where \f$ c_{\underline{k}} \f$ is a value associated with the node \f$ \underline{k} \f$
 * and where \f$ \underline{e}^d \f$ is the unit vector in direction \f$ d \f$.
 * To simplify such algorithms, the class ShapeEnum organizes a shape into _slices_. 
 * The \f$ s \f$-th slice of a shape \f$ \mathfrak{K} \f$ contains all nodes \f$ 
 * \underline{k} \in \mathfrak{K} \f$ that satisfy \f$ \sum_{d=1}^{D} k_d = s \f$.
 * 
 * Nodes in the same slice are ordered lexically. This ordering enables simple and efficient
 * union- and intersect-operations on shape enumerations.
 * 
 * ## Usage ##
 * 
 * \code{.cpp}
 * #include "waveblocks/shape_hyperbolic.hpp"
 * 
 * #include "waveblocks/tiny_multi_index.hpp"
 * #include "waveblocks/shape_enumerator.hpp"
 * #include "waveblocks/shape_enum.hpp"
 * 
 * using namespace waveblocks;
 * \endcode
 * 
 * Create a suitable \ref BasicShape "shape" .
 * \code{.cpp}
 * const dim_t D = 5;
 * LimitedHyperbolicCutShape<D> shape(7.0, {2,2,4,4,4});
 * \endcode
 * 
 * Select an appropriate data type to represent an integer tuple.
 * \code{.cpp}
 * typedef TinyMultiIndex<std::size_t,D> MultiIndex;
 * \endcode
 * Pass it to a _shape enumerator_ and you get an _shape enumeration_:
 * \code{.cpp}
 * ShapeEnumerator<D, MultiIndex> enumerator;
 * ShapeEnum<D, MultiIndex> enumeration = enumerator.generate(shape);
 * \endcode
 * Select a slice:
 * \code 
 * const ShapeSlice<D, MultiIndex>& slice = enumeration.slice(slice_index);
 * \endcode
 * Get \f$ i \f$-th node of current slice:
 * \code{.cpp}
 * MultiIndex index = slice[i];
 * \endcode
 * Get position (aka ordinal) of a node \f$ \underline{k} \f$, if you know
 * that this node _is part of the shape_, ...
 * \code{.cpp}
 * std::size_t ordinal = slice.find(k);
 * \endcode
 * ..., if not, use:
 * \code{.cpp}
 * std::size_t ordinal;
 * if (slice.try_find(k, ordinal)) {}
 * \endcode
 * If you need to know the ordinals of all backward neighbours of a parent node, use
 * find_backward_neighbours().
 * Ensure that you consult the slice of the backward neighbours (and not the slice of the parent node).
 * Pay attention that this method has an important restriction: The shape must contain 
 * the parent node.
 * \code{.cpp}
 * std::array<std::size_t,D> ordinals = slice.find_backward_neighbours(node);
 * \endcode
 */

}