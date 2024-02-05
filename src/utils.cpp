//
// Created by Asunción Gómez on 18-01-24.
//
#include <bits/stdc++.h>
#include <sdsl/bit_vectors.hpp>

using namespace sdsl;
using namespace std;


/**
 * Given a bitvector and the number of bits to represent the children of a node, return the coordinates of the point in the grid.
 * Concatenating the labels in the path down to the integer yields the bit-string.
 * For example, for a quadtree we'll represent a node with 2 bits.
 * If the path to a leaf is, for example, 0010, then the coordinates of the leaf are a=1="01" and b=0="00" in the grid.
 * @param bv
 * @param l
 * @return an array of l dimension with the point coordinates in the grid according to the bit path.
 */
uint64_t* getCoordinates(std::bitset<64> bv, uint16_t l){
    uint64_t* pointCoord = new uint64_t[l];
    for(uint16_t i = 0 ; i < l; i++){

    }

    return pointCoord;
}

