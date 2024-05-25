//
// Created by Asunción Gómez on 19-04-24.
//

#include "../lqdag.hpp"

// ---------- FULL RELATIONAL ALGEBRA ---------- //

// (JOIN, L1(A1), L2(A2)) = (AND, (EXTEND, L1, A1 ∪ A2), (EXTEND, L2, A1 ∪ A2))
lqdag* lqdag_join(vector<qdag> &Q, bool bounded_result, uint64_t UPPER_BOUND) {
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

    for (uint64_t i = 0; i < Q.size(); i++) {
        subQ[i] = new subQuadtreeChild{&Q[i], 0, 0};
        extend_qdag[i] = new lqdag(FUNCTOR_EXTEND, new lqdag(FUNCTOR_QTREE, subQ[i]), A);
        if (A.size() == 3)
            Q[i].createTableExtend3();
        else if (A.size() == 4)
            Q[i].createTableExtend4();
        else if (A.size() == 5)
            Q[i].createTableExtend5();
        else {
            cout << "Code only works for queries of up to 5 attributes..." << endl;
            exit(1);
        }
    }

    lqdag* join_result;
    if (Q.size() == 3){
        lqdag* join_r_s = new lqdag(FUNCTOR_AND, extend_qdag[0], extend_qdag[1]);
        join_result = new lqdag(FUNCTOR_AND, join_r_s, extend_qdag[2]);
    }
    else if (Q.size() == 4){
        lqdag* join_r_s = new lqdag(FUNCTOR_AND, extend_qdag[0], extend_qdag[1]);
        lqdag* join_t_u = new lqdag(FUNCTOR_AND, extend_qdag[2], extend_qdag[3]);
        join_result = new lqdag(FUNCTOR_AND, join_r_s, join_t_u);
    }
    else if (Q.size() == 5){
        lqdag* join_r_s = new lqdag(FUNCTOR_AND, extend_qdag[0], extend_qdag[1]);
        lqdag* join_t_u = new lqdag(FUNCTOR_AND, extend_qdag[2], extend_qdag[3]);
        lqdag* join_r_s_t_u = new lqdag(FUNCTOR_AND, join_r_s, join_t_u);
        join_result = new lqdag(FUNCTOR_AND, join_r_s_t_u, extend_qdag[4]);
    }
    else {
        cout << "Code only works for queries of up to 5 attributes..." << endl;
        exit(1);
    }

    return join_result;
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

