
namespace waveblocks {

/**
 * \defgroup BasicShape Shape Definition
 * \ingroup WaveBlocksND
 * 
 * \brief Provides definitions for commonly used shapes.
 * 
 * A shape definition class must expose two methods:
 * 
 * __Surface__
 * \code
 * template<class MultiIndex>
 * int limit(MultiIndex base_node, dim_t axis) const;
 * \endcode
 * 
 * > Returns for a given \e axis \f$ j \f$ and a given \e base \e node \f$ \boldsymbol{n} \f$
 * > the largest element \f$ k^\star \f$ that satisfies: 
 * > \f[ \boldsymbol{k} \in \mathfrak{K}, \;
 * k_i =
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
 * \code
 * int bbox(dim_t axis) const;
 * \endcode
 * 
 * > Returns for a given \e axis \f$ j \f$ a (preferably) as small as possible \e limit \f$ L_j \f$
 * > such that:
 * > \f[
 * \forall \boldsymbol{k} \in \mathfrak{K} \,\colon\; k_j \leq L_j
 * \f]
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
 */

}