//
// Created by Asunción Gómez on 19-04-24.
//

//#include "../lqdag.hpp"

// ---------- FULL RELATIONAL ALGEBRA ---------- //

// (JOIN, L1(A1), L2(A2)) = (AND, (EXTEND, L1, A1 ∪ A2), (EXTEND, L2, A1 ∪ A2))
/**
 * Formula for the join operation
 * @param Q
 * @param UPPER_BOUND
 * @return
 */
lqdag* lqdag_join(vector<qdag> &Q, uint64_t UPPER_BOUND) {
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

/**
 * Completion for the join operation
 * @param Q
 * @param UPPER_BOUND
 * @param results
 * @return
 */
quadtree_formula* compute_lqdag_join(vector<qdag> &Q, uint64_t UPPER_BOUND, uint64_t &results) {
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

    uint64_t p = std::pow(Q[0].getK(), A.size());

    // join result is the formula
    return join_result->completion(p, Q[0].getHeight(), 0, UPPER_BOUND, results);
}

/**
 * Logical disjunction (formula)
 * @param Q
 * @return
 */
lqdag* lqdag_or(vector<qdag> &Q) {
    subQuadtreeChild* subQ[Q.size()];
    lqdag* or_qdag[Q.size()];

    for (uint64_t i = 0; i < Q.size(); i++) {
        subQ[i] = new subQuadtreeChild{&Q[i], 0, 0};
        or_qdag[i] = new lqdag(FUNCTOR_QTREE, subQ[i]);
    }

    lqdag* or_result;
    if (Q.size() == 3){
        lqdag* or_r_s = new lqdag(FUNCTOR_OR, or_qdag[0], or_qdag[1]);
        or_result = new lqdag(FUNCTOR_OR, or_r_s, or_qdag[2]);
    }
    else if (Q.size() == 4){
        lqdag* or_r_s = new lqdag(FUNCTOR_OR, or_qdag[0], or_qdag[1]);
        lqdag* or_t_u = new lqdag(FUNCTOR_OR, or_qdag[2], or_qdag[3]);
        or_result = new lqdag(FUNCTOR_OR, or_r_s, or_t_u);
    }
    else if (Q.size() == 5){
        lqdag* or_r_s = new lqdag(FUNCTOR_OR, or_qdag[0], or_qdag[1]);
        lqdag* or_t_u = new lqdag(FUNCTOR_OR, or_qdag[2], or_qdag[3]);
        lqdag* or_r_s_t_u = new lqdag(FUNCTOR_OR, or_r_s, or_t_u);
        or_result = new lqdag(FUNCTOR_OR, or_r_s_t_u, or_qdag[4]);
    }
    else {
        cout << "Code only works for queries of up to 5 attributes..." << endl;
        exit(1);
    }

    return or_result;
}

lqdag* lqdag_and(vector<qdag> &Q){
    // TODO
}

lqdag* lqdag_not(vector<qdag> &Q){
    // TODO
}


/**
 * Completion of a predicate over a qdag
 * @param Q
 * @param UPPER_BOUND
 * @param results
 * @param pred
 * @return
 */
quadtree_formula* compute_pred_formula(qdag Q, uint64_t UPPER_BOUND, uint64_t &results, predicate* pred){
    qdag::att_set A;
    map<uint64_t, uint8_t> attr_map;

    uint64_t nAttr = Q.nAttr();
    for (uint64_t j = 0; j < nAttr; j++)
        attr_map[Q.getAttr(j)] = 1;

    // set of attributs
    for (map<uint64_t, uint8_t>::iterator it = attr_map.begin(); it != attr_map.end(); it++)
        A.push_back(it->first);

    uint64_t p = std::pow(Q.getK(), A.size());

    subQuadtreeChild* subQ = new subQuadtreeChild{&Q, 0, 0};
    lqdag* lqdag_formula = new lqdag(FUNCTOR_QTREE, subQ);

    uint16_t coordinates[A.size()];
    for(uint16_t i = 0; i < A.size(); i++)
        coordinates[i] = 0;

    return lqdag_formula->completion_with_pred(p,Q.getHeight(),0,0,UPPER_BOUND,results,pred,coordinates);
}


// TODO: completion pred over a more complex formula (two or more relations) R(A,B) and S(B,C) and pred : A=C

// (DIFF, L1(A), L2(A)) = (AND, L1, (NOT, L2))
lqdag* diff(qdag q1, qdag q2, bool bounded_result, uint64_t UPPER_BOUND) {
    subQuadtreeChild* subQ1 = new subQuadtreeChild{&q1, 0, 0};
    subQuadtreeChild* subQ2 = new subQuadtreeChild{&q2, 0, 0};
    lqdag* diff = new lqdag(FUNCTOR_AND, new lqdag(FUNCTOR_QTREE, subQ1), new lqdag(FUNCTOR_NOT, subQ2));
    return diff;
}

quadtree_formula* compute_lqdag_diff(qdag q1, qdag q2, bool bounded_result, uint64_t UPPER_BOUND,uint64_t &results){
    subQuadtreeChild* subQ1 = new subQuadtreeChild{&q1, 0, 0};
    subQuadtreeChild* subQ2 = new subQuadtreeChild{&q2, 0, 0};
    lqdag* diff = new lqdag(FUNCTOR_AND, new lqdag(FUNCTOR_QTREE, subQ1), new lqdag(FUNCTOR_NOT, subQ2));
    uint64_t p = std::pow(q1.getK(), q1.nAttr());
    return diff->completion(p, q1.getHeight(), 0, UPPER_BOUND, results);
}


/**
 * Takes a subexpression F and a predicate θ, which is a logical expression on the attributes of F.
 */
lqdag* selection(){
    // TODO
}

lqdag* projection(){
    // TODO
}

lqdag* semijoin(){

}

lqdag* antijoin(){

}

lqdag* division(){

}

