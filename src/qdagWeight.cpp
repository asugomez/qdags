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
    int16_t level; // the level of the node
    uint64_t* roots; // the parent of the subtree of each qdag that conform the tuple.
    uint64_t weight; // priority or number of leaves of the tuple
    uint64_t path;  // the bits that encode the path down the leaf in the first qdag.
    bool operator<(const qdagWeight& qd) const {
        return weight < qd.weight;
    }
};

// TODO: ver como esto ayuda para num_leaves
// TODO: ver como pasar esto a coordenadas en la grilla