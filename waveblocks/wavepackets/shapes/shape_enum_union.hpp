#pragma once

#include <algorithm>
#include <vector>
#include <utility>
#include <initializer_list>

#include "shape_enum.hpp"


namespace waveblocks {
    namespace wavepackets {
        namespace shapes {
            namespace shape_enum {
                /**
                 * \ingroup ShapeEnum
                 *
                 * \brief Creates union of two shape slices.
                 *
                 * \f[ \forall k \colon k \in union(\mathfrak{K}_1,\mathfrak{K}_2) \iff k \in \mathfrak{K}_1 \lor k \in \mathfrak{K}_2 \f]
                 *
                 * \param[in] begin1 input iterator to begin of lattice point vector of first slice
                 * \param[in] end1 input iterator to end of lattice point vector of first slice
                 * \param[in] begin2 input iterator to begin of lattice point vector of second slice
                 * \param[in] end2 input iterator to end of lattice point vector of second slice
                 * \param[out] sink output iterator to store united slice
                 * \param[in] less function object that implements \e less operation on multi-indices
                 * \tparam Input1 type of first input iterator
                 * \tparam Input2 type of second input iterator
                 * \tparam Output type of output iterator
                 * \tparam Compare type of function object that implements \e less operation on multi-indices
                 * \return union of both slices
                 */
                template<class Input1, class Input2, class Output, class Compare>
                Output strict_union(Input1 begin1, Input1 end1, Input2 begin2, Input2 end2, Output sink, Compare less)
                {
                    auto it1 = begin1;
                    auto it2 = begin2;

                    while (it1 != end1 && it2 != end2) {
                        if (*it1 == *it2) {
                            // lhs = rhs
                            *sink = *it1;
                            ++it1;
                            ++it2;
                            ++sink;
                        }
                        else if ( less(*it1, *it2) ) {
                            // lhs < rhs
                            *sink = *it1;
                            ++it1;
                            ++sink;
                        }
                        else {
                            // lhs > rhs
                            *sink = *it2;
                            ++it2;
                            ++sink;
                        }
                    }

                    sink = std::copy(it1, end1, sink);
                    sink = std::copy(it2, end2, sink);

                    return sink;
                }

                template<class MultiIndex>
                struct _strict_union__heap_entry
                {
                public:
                    MultiIndex value;
                    std::size_t source;

                    bool operator<(const _strict_union__heap_entry& lhs) const
                    {
                        // this operator is used for heap construction
                        // the heap puts the largest value on top but we want the smallest value on top because
                        // the resulting enumeration must be sorted in ascending order
                        // therefore this operator does the opposite, namely 'LHS >= RHS' instead of 'LHS < RHS'
                        std::less<MultiIndex> less;
                        return less(lhs.value, value);
                    }
                };

                /**
                 *
                 */
                template<class MultiIndex>
                std::vector<MultiIndex> strict_union(std::vector< typename std::vector<MultiIndex>::const_iterator > begin,
                                                     std::vector< typename std::vector<MultiIndex>::const_iterator > end)
                {
                    assert( begin.size() == end.size() );

                    typedef _strict_union__heap_entry<MultiIndex> HeapEntry;

                    // get size of largest source
                    std::size_t minsize = 0;
                    for (std::size_t i = 0; i < begin.size(); i++) {
                        minsize = std::max(minsize, static_cast<std::size_t>(end[i] - begin[i]));
                    }

                    std::vector<MultiIndex> superset;
                    std::vector< HeapEntry > heap; // tuple (multi-index, source)

                    // initialize heap
                    for (std::size_t i = 0; i < begin.size(); i++) {
                        if (begin[i] != end[i]) {
                            heap.push_back( HeapEntry{*(begin[i]), i} );
                            std::push_heap(heap.begin(), heap.end());
                            (begin[i])++;
                        }
                    }
                    std::make_heap(heap.begin(), heap.end());

                    // merge multi-indices
                    while (!heap.empty()) {
                        HeapEntry entry = heap.front();
                        std::pop_heap(heap.begin(), heap.end()); heap.pop_back();

                        if (superset.empty() || entry.value != superset.back()) {
                            std::less<MultiIndex> less;

                            (void)less;
                            assert( superset.empty() || less(superset.back(), entry.value) ); // multi-indices must be sorted in ascending order
                            superset.push_back(entry.value);
                        }

                        if ( begin[entry.source] != end[entry.source] ) {
                            heap.push_back( HeapEntry{*(begin[entry.source]), entry.source} );
                            std::push_heap(heap.begin(), heap.end());
                            (begin[entry.source])++;
                        }
                    }

                    return superset;
                }

                template<dim_t D, class MultiIndex>
                ShapeSlice<D, MultiIndex> strict_union(const std::vector< const ShapeSlice<D, MultiIndex>* >& slices, std::size_t union_offset)
                {
                    std::vector< typename std::vector<MultiIndex>::const_iterator > begin(slices.size());
                    std::vector< typename std::vector<MultiIndex>::const_iterator > end(slices.size());

                    for (std::size_t i = 0; i < slices.size(); i++) {
                        begin[i] = slices[i]->cbegin();
                        end[i] = slices[i]->cend();
                    }

                    std::vector<MultiIndex> superset = strict_union<MultiIndex>(begin, end);

                    return {std::move(superset), union_offset};
                }

                template<dim_t D, class MultiIndex>
                ShapeEnum<D, MultiIndex> strict_union(std::vector< ShapeEnum<D, MultiIndex> const* > const& enums)
                {
                    // determine number of slices in union
                    int n_slices = 0;
                    MultiIndex limits{};
                    for (auto _enum : enums) {
                        n_slices = std::max(n_slices, _enum->n_slices());
                        for (dim_t d = 0; d < D; d++) {
                            limits[d] = std::max((int)limits[d], _enum->limit(d));
                        }
                    }

                    // create union
                    std::size_t offset = 0;
                    std::vector< ShapeSlice<D,MultiIndex> > superset(n_slices);
                    for (int islice = 0; islice < n_slices; islice++) {
                        std::vector< const ShapeSlice<D, MultiIndex>* > slices;
                        for (std::size_t source = 0; source < enums.size(); source++) {
                            auto slice = &enums[source]->slice(islice);
                            if (slice->size() != 0) {
                                slices.push_back(slice);
                            }
                        }
                        superset[islice] = strict_union<D,MultiIndex>(slices, offset);
                        offset += superset[islice].size();
                    }

                    return {std::move(superset), offset, limits};
                }

                /**
                 * \ingroup ShapeUnion
                 *
                 * \param enum_list
                 */
                template<dim_t D, class MultiIndex, std::size_t N>
                ShapeEnum<D, MultiIndex> strict_union(const std::array< ShapeEnum<D, MultiIndex>*, N >& enum_list)
                {
                    std::vector< const ShapeEnum<D, MultiIndex>* > enum_vec(enum_list.cbegin(), enum_list.cend());
                    return strict_union<D,MultiIndex>(enum_vec);
                }

                /**
                 * \ingroup ShapeUnion
                 * \param enum_list
                 */
                template<dim_t D, class MultiIndex, std::size_t N>
                ShapeEnum<D, MultiIndex> strict_union(const std::array< const ShapeEnum<D, MultiIndex>*, N >& enum_list)
                {
                    std::vector< const ShapeEnum<D, MultiIndex>* > enum_vec(enum_list.cbegin(), enum_list.cend());
                    return strict_union<D,MultiIndex>(enum_vec);
                }

                /**
                 * \ingroup ShapeUnion
                 * \param enum_list
                 */
                template<dim_t D, class MultiIndex>
                ShapeEnum<D, MultiIndex> strict_union(std::initializer_list< const ShapeEnum<D, MultiIndex>* > enum_list)
                {
                    std::vector< const ShapeEnum<D, MultiIndex>* > enum_vec(enum_list);
                    return strict_union<D,MultiIndex>(enum_vec);
                }
            }
        }
    }
}
