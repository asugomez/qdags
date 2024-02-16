//
// Created by Asunción Gómez on 18-01-24.
//
#include <bits/stdc++.h>
//#include <sdsl/bit_vectors.hpp>


using namespace std;

// TODO: otra opción es precomputar tablas con mascaras de bits 100100100, 010010010, etc.
// TODO: analyze the uint16,32,64..
/**
 * Given a bitvector and the number of bits to represent the children of a node, return the coordinates of the point in the grid.
 * Concatenating the labels in the path down to the integer yields the bit-string.
 * For example, for a quadtree we'll represent a node with 2 bits.
 * If the path to a leaf is, for example, 0010, then the coordinates of the leaf are a=1="01" and b=0="00" in the grid.
 * @param path
 * @param l bits number to define a node.
 * @return an array of l dimension with the point coordinates in the grid according to the bit path.
 */
void getCoordinates(uint64_t path, uint64_t l, uint64_t height, uint64_t *pointCoord, uint64_t nAttr) {

    uint64_t diff_level;
    for(uint64_t level=0; level<height; level++) {
        diff_level = height - level;
        uint32_t node_cur_level = (uint32_t) path >> (diff_level * l);
        uint64_t n_ones = bits::cnt((uint64_t) node_cur_level);

        uint64_t msb, num_point;
        for (uint64_t k = 0; k < n_ones; k++) {
            msb = __builtin_ctz(node_cur_level);
            node_cur_level &= ~(1 << msb); // delete the msb of children
            msb %= l;
            num_point = (l - 1) - msb; // the index is the inverse
            pointCoord[num_point] += 1 << diff_level;
        }
    }

}

