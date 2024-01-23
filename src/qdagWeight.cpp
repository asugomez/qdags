//
// Created by Asunción Gómez on 24-12-23.
//
#include <algorithm>
#include "qdags.hpp"

/**
 * Nos indica la posicion de la tupla que nos encontramos del qdag
 * Example: level = -1, node = 0, weight = 30 --> prioridad de la raiz
 */
struct qdagWeight {
    uint16_t level; // the level of the node
    uint64_t node; // the i-th non-empty node (1) of that level (the quadrant)
    uint64_t weight; // priority or number of leaves
    vector<bool> bv;  // the bits that encode the path down the leaf.
    bool operator<(const qdagWeight& qd) const {
        return weight < qd.weight;
    }
};

// TODO: ver como esto ayuda para num_leaves
// TODO: ver como pasar esto a coordenadas en la grilla