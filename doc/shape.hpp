
namespace waveblocks {

/**
 * \defgroup BasicShape Shape Definition
 * \ingroup WaveBlocksND
 * 
 * \brief Provides definitions for commonly used shapes.
 * 
 * ## Definition ##
 * A \f$ D \f$-dimensional shape \f$ \mathfrak{K} \f$ is a set of _unordered_ D-dimensional integer-tuples (aka _node_).
 * A shape is suitable for our needs if it satisfies:
 * \f[
 * k \in \mathfrak{K} \Rightarrow \forall k^\star \preceq k \colon k^\star \in \mathfrak{K}
 * \f]
 * That means, if an arbitrary node is part of the shape, then all nodes in the backward cone are part of the shape too.
 * 
 * ## Implementation ##
 * To define a custom shape definition, your class must expose two methods:
 * 
 * __Surface__
 * \code{.cpp}
 * template<class MultiIndex>
 * int limit(MultiIndex base_node, dim_t axis) const;
 * \endcode
 * 
 * > Returns for a given \e axis \f$ j \f$ and a given \e base \e node \f$ \underline{n} \f$
 * > the largest element \f$ k^\star \f$ that satisfies \f$ \underline{k}(k^\star) \in \mathfrak{K} \f$: 
 * > \f[ k_i(k^\star) =
 *    \begin{cases}
 *       n_i,& i \neq j\\
 *       k^\star, & i = j
 *    \end{cases}
 * \f]
 * > 
 * > A \e base \e node is a node whose \f$ j \f$ -th entry is zero.
 * > 
 * > Notice that \f$ n_j \f$ does not influence return value \f$ k^\star \f$. 
 * > Therefore \f$ n_j \f$ can be any value since it is simply ignored.
 * 
 * __Bounding Volume__
 * \code{.cpp}
 * int bbox(dim_t axis) const;
 * \endcode
 * 
 * > Returns for a given \e axis \f$ j \f$ a (preferably) as small as possible \e limit \f$ L_j \f$
 * > such that:
 * > \f[
 * \forall \underline{k} \in \mathfrak{K} \,\colon\; k_j \leq L_j
 * \f]
 * 
 * ## Usage ##
 * A shape instance itself is quite useless. It can tell you whether it contains a specific node:
 * \code{.cpp}
 * bool contains = test_node[d] <= shape.limit(test_node,d); // d can be any axis
 * \endcode
 * A \ref ShapeEnum "Shape Enumeration" is extremely more useful since it orders
 * all nodes of a shape.
 */

/**
 * \defgroup ShapeExtension Shape Extension
 * \ingroup WaveBlocksND
 * 
 * \brief Provides methods to extend shapes.
 * 
 * __Def:__ Given a basis shape \f$ \mathfrak{K} \f$, the shape extension
 * \f$ \overline{\mathfrak{K}} \f$ is defined by:
 * \f[ \overline{\mathfrak{K}} := \mathfrak{K} \cup 
 * \left\{\underline{k}' \colon \underline{k}' = \underline{k} + \underline{e_d} 
 * \forall d \in \{1,\ldots,D\} \forall \underline{k} \in \mathfrak{K}\right\} \f]
 * 
 * __Usage:__ The shape extension is used to compute the gradient of a hagedorn wavepacket.
 * 
 * ## Usage ##
 * 
 * ### Creation at compile-time ###
 * 
 * \code{.cpp}
 * const dim_t D = 5;
 * 
 * typedef HyperCubicShape<D> S;
 * 
 * S shape({4,4,4,2,2});
 * 
 * ExtendedShape<D, S> shape_extension(shape);
 * 
 * typedef TinyMultiIndex<std::size_t, D> MultiIndex;
 * 
 * ShapeEnumerator<D,MultiIndex> shape_enumerator;
 * ShapeEnum<D,MultiIndex> enum_extension = enumerator.generate(shape_extension);
 * \endcode
 * 
 * ### Creation at runtime  ###
 * 
 * \code{.cpp}
 * // Not implemented yet
 * \endcode
 */

/**
 * 
 * \defgroup ShapeUnion Shape Union
 * \ingroup WaveBlocksND
 * 
 * \brief Provides methods to create unions of shape definitions and shape enumerations.
 * 
 * __Def:__ The union of the shapes \f$ \mathfrak{K}_1,\ldots,\mathfrak{K}_N \f$
 * is defined by:
 * \f[
 * union(\mathfrak{K}_1,\ldots,\mathfrak{K}_N) := \left\{ 
 *      k \mid \exists i \in \{1, \ldots, N \} \colon k \in \mathfrak{K}_i
 * \right\}
 * \f]
 * 
 * __Usage:__ Efficient evaluation of homogeneous wavepackets.
 * 
 * ## Usage ##
 * 
 * ### Creation at compile-time ###
 * 
 * \code{.cpp}
 * const dim_t D = 5;
 * 
 * typedef HyperCubicShape<D> S1;
 * typedef LimitedHyperbolicCutShape<D> S2;
 * 
 * S1 shape1({4,4,4,2,2});
 * S2 shape2(9.0, {4,4,4,4,4});
 * SupersetShape<D, S1, S2> shape_union(shape1, shape2);
 * 
 * typedef TinyMultiIndex<std::size_t, D> MultiIndex;
 * 
 * ShapeEnumerator<D,MultiIndex> enumerator;
 * ShapeEnum<D,MultiIndex> enum1 = enumerator.generate(shape1);
 * ShapeEnum<D,MultiIndex> enum2 = enumerator.generate(shape2);
 * ShapeEnum<D,MultiIndex> enum_union = enumerator.generate(shape_union);
 * \endcode
 * 
 * ### Creation at runtime  ###
 * 
 * \code{.cpp}
 * const dim_t D = 5;
 * 
 * typedef HyperCubicShape<D> S1;
 * typedef LimitedHyperbolicCutShape<D> S2;
 * 
 * S1 shape1({4,4,4,2,2});
 * S2 shape2(9.0, {4,4,4,4,4});
 * 
 * typedef TinyMultiIndex<std::size_t, D> MultiIndex;
 * 
 * ShapeEnumerator<D,MultiIndex> enumerator;
 * ShapeEnum<D,MultiIndex> enum1 = enumerator.generate(shape1);
 * ShapeEnum<D,MultiIndex> enum2 = enumerator.generate(shape2);
 * ShapeEnum<D,MultiIndex> enum_union = shape_enum::strict_union({&enum1, &enum2});
 * \endcode
 */

}