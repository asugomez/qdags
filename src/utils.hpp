//
// Created by Asunción Gómez on 19-02-24.
//

#ifndef QDAGS_UTILS_HPP
#define QDAGS_UTILS_HPP

#include<ctime>
#include<bits/stdc++.h>
#include <sdsl/bit_vectors.hpp>


struct qdagWeight {
    uint16_t level; // the level of the node
    uint64_t* roots; // the parent of the subtree of each qdag that conform the tuple.
    double weight; // priority or number of leaves of the tuple
    uint64_t path;  // the bits that encode the path down the leaf in the first qdag.
    bool operator<(const qdagWeight& qd) const {
        return weight < qd.weight;
    }
};

/**
 * Get the coordinates of a point in the grid given the path from the root to the leaf in the quadtree.
 * This function is O(m), with m the number of 1s in the path.
 * @param path see the paper for the definition (page 5).
 * @param l number of bits to define a child on each level
 * @param height height of the quadtree
 * @param pointCoord the coordinates of the point in the grid. We are going to write on this array.
 */
void getCoordinates(uint64_t path, uint16_t l, uint64_t height, uint32_t *pointCoord) {

    uint32_t bits_path = (height + 1) * l;
    path = (uint32_t) path << (31 - bits_path + 1); // move the beginning of the path to the begining of the 32 bits
    uint32_t n_ones = sdsl::bits::cnt(path);
    uint32_t msb, level;
    // we only pass across the 1s to update the coordinates.
    for (uint32_t k = 0; k < n_ones; k++) {
        msb = __builtin_clz(path);
        level = msb/l;
        path &= (((uint32_t) 0xffffffff) >> (msb + 1)); // delete the msb of children
        uint32_t num_point = msb%l;
        pointCoord[num_point] += 1 << (height - level);
    }

}


#endif //QDAGS_UTILS_HPP
