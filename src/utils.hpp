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
    //uint64_t path;  // the bits that encode the path down the leaf in the first qdag.
    uint16_t* coordinates;
    bool operator<(const qdagWeight &qdag) const {
        return weight < qdag.weight; // max heap
    }
};


struct qdagResults{ // we store the output of the join
    uint16_t* coordinates;
    //uint64_t path;
    double weight; // priority or number of leaves of the tuple
    bool operator<(const qdagResults &qdag) const {
        return weight > qdag.weight; // min heap
    }
};

struct orderJoinQdag{
    uint64_t index;
    //uint64_t path;
    uint16_t* coordinates;
    double weight;
    bool operator<(const orderJoinQdag &ojq) const {
        return weight < ojq.weight; // max heap
    }
};
//
///**
// * Get the coordinates of a point in the grid given the path from the root to the leaf in the quadtree.
// * This function is O(m), with m the number of 1s in the path.
// * @param path see the paper for the definition (page 5).
// * @param l number of bits to define a child on each level
// * @param height height of the quadtree
// * @param pointCoord the coordinates of the point in the grid. We are going to write on this array.
// */
//void getCoordinates(uint64_t path, uint16_t nAtt, uint64_t height, uint16_t *pointCoord) {
//    uint64_t n_ones = sdsl::bits::cnt(path);
//    uint64_t msb, level;
//    // we only pass across the 1s to update the coordinates.
//    for (uint64_t k = 0; k < n_ones; k++) {
//        msb = __builtin_clz(path);
//        level = (31-msb)/nAtt;
//        pointCoord[(msb+1)%nAtt] += 1 << level;
//        path &= ((0xffffffff) >> (msb + 1)); // delete the msb of children
//    }
//}

/**
 * Transform the i-th child into a coordinate according to the Morton code.
 * When we descend a level, we compute the new coordinates of the point.
 * @param coordinates
 * @param nAtt
 * @param diff_level
 * @param child
 */
void transformCoordinates(uint16_t* coordinates, uint16_t nAtt, uint64_t diff_level, uint16_t child){
    uint16_t n_ones = sdsl::bits::cnt(child);
    uint64_t msb;
    for(uint16_t i=0; i < n_ones; i++){
        msb = __builtin_clz(child);
        coordinates[(msb+1)%nAtt] += 1 << diff_level; // the last bits correspond to the initial coordinates : [....,z,w,y,x]
        child &= ((0xffffffff) >> (msb + 1));
    }

}


#endif //QDAGS_UTILS_HPP
