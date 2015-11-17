/* Author: R. Bourquin
 * Copyright: (C) 2014 R. Bourquin
 * License: GNU GPL v2 or above
 *
 * A library of helper functions to search for
 * Kronrod extensions of Gauss quadrature rules.
 */

#pragma once

#include <array>
#include <list>

template<int D> using partition_t = std::array<int, D>;
template<int D> using partitions_t = std::list<partition_t<D>>;
template<int D> using point_t = std::array<int, D>;

template<int D>
int sum(const point_t<D> Z) {
    /* Sum of the components of the point Z.
     */
    int s = 0;
    for(int i=0; i<D; i++) {
        s += Z[i];
    }
    return s;
}

template<int D>
partitions_t<D> partitions(const int K) {
    /*
     * Enumerate integer partitions in anti-lexocographic
     * order for integers up to some limit K. All partitions
     * have exactly D parts, some may be zero.
     *
     * K: Level
     */
    partitions_t<D> partitions;
    partition_t<D> P;

    P.fill(0);
    partitions.push_back(P);

    while(sum<D>(P) <= K) {
        int p0 = P[0];
        bool broke = false;
        for(int i=1; i < D; i++) {
            p0 += P[i];
            if(P[0] <= P[i] + 1) {
                P[i] = 0;
            } else {
                P[0] = p0 - i * (P[i] + 1);
                for(int j=1; j <= i; j++) {
                    P[j] = P[i] + 1;
                }
                partitions.push_back(P);
                broke = true;
                break;
            }
        }
        if(!broke) {
            P[0] = p0 + 1;
            if(sum<D>(P) <= K) {
                partitions.push_back(P);
            }
        }
    }

    return partitions;
}
