//
// Created by Asunción Gómez on 19-04-24.
//

#include "../lqdag.hpp"

// ---------- FULL RELATIONAL ALGEBRA ---------- //

// (JOIN, L1(A1), L2(A2)) = (AND, (EXTEND, L1, A1 ∪ A2), (EXTEND, L2, A1 ∪ A2))
lqdag* join(vector<qdag> &Q, bool bounded_result, uint64_t UPPER_BOUND) {
    qdag::att_set A;
    map<uint64_t, uint8_t> attr_map;

    for (uint64_t i = 0; i < Q.size(); i++) {
        uint64_t nAttr = Q[i].nAttr();
        for (uint64_t j = 0; j < nAttr; j++)
            attr_map[Q[i].getAttr(j)] = 1;
    }

    // set of attributs
    for (map<uint64_t, uint8_t>::iterator it = attr_map.begin(); it != attr_map.end(); it++)
        A.push_back(it->first);

    subQuadtreeChild* subQ[Q.size()];
    lqdag* extend_qdag[Q.size()];
    for(uint64_t i = 0; i < Q.size(); i++) {
        subQ[i] = new subQuadtreeChild{&Q[i], 0, 0};
        extend_qdag[i] = new lqdag(FUNCTOR_EXTEND, new lqdag(FUNCTOR_QTREE, subQ[i]), A);
    }
    bool is_odd = Q.size() % 2 == 0;

    lqdag* and_extend[Q.size()/2 + (Q.size() % 2)];
    for(uint64_t i = 0; i < Q.size()/2; i++) {
        and_extend[i] = new lqdag(FUNCTOR_AND, extend_qdag[2*i], extend_qdag[2*i+1]);
    }

    if(!is_odd)
        and_extend[Q.size()/2] = extend_qdag[Q.size()-1];

    // TODO: completar el arbol de ANDs (la altura es log_2(Q.size))

}

// (DIFF, L1(A), L2(A)) = (AND, L1, (NOT, L2))
lqdag* diff(qdag q1, qdag q2, bool bounded_result, uint64_t UPPER_BOUND) {
    subQuadtreeChild* subQ1 = new subQuadtreeChild{&q1, 0, 0};
    subQuadtreeChild* subQ2 = new subQuadtreeChild{&q2, 0, 0};
    lqdag* diff = new lqdag(FUNCTOR_AND, new lqdag(FUNCTOR_QTREE, subQ1), new lqdag(FUNCTOR_NOT, subQ2));
    return diff;
}

// is a virtual lqdag on A whose cells are exactly those that satisfy the predicate θ .
void pred(){

}

void selection(){

}

