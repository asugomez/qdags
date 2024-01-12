//
// Created by Asunción Gómez on 24-12-23.
//
#include <algorithm>
#include "qdags.hpp"

struct qdagPri {
    uint16_t level;
    uint64_t node;
    uint64_t pri;
    bool operator<(const qdagPri& qd) const {
        return pri < qd.pri;
    }
};