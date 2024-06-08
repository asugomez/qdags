//
// Created by Asunción Gómez on 19-04-24.
//

//#include "../lqdag.hpp"

// ---------- QUADTREE  ---------- //

quadtree_formula* compute_quadtree(qdag Q, uint64_t UPPER_BOUND, uint64_t &results){
    subQuadtreeChild* subQ = new subQuadtreeChild{&Q, 0, 0};
    lqdag* lqdag_formula = new lqdag(FUNCTOR_QTREE, subQ);
    return lqdag_formula->completion(Q.getKD(), Q.getHeight()-1, 0, UPPER_BOUND, results);
}


/**
 * Completion of a predicate over a qdag
 * @param Q
 * @param UPPER_BOUND
 * @param results
 * @param pred
 * @return
 */
quadtree_formula* compute_pred_quadtree(qdag Q, uint64_t UPPER_BOUND, uint64_t &results, predicate* pred){
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

    uint64_t grid_side = Q.getGridSide();
    uint8_t k = Q.getK();
    uint16_t coordinates[A.size()];
    for(uint16_t i = 0; i < A.size(); i++)
        coordinates[i] = 0;

    return lqdag_formula->completion_with_pred(p, Q.getHeight()-1, 0, grid_side, k, UPPER_BOUND, results, pred, coordinates, A.size());
}

// ---------- JOIN ---------- //

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
    }

    lqdag* join_result;
    if(Q.size() == 2){
        join_result = new lqdag(FUNCTOR_AND, extend_qdag[0], extend_qdag[1]);
    }
    else if (Q.size() == 3){
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
        cout << "Code only implemented for a maximum of 5 qdags..." << endl;
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
    }

    lqdag* join_result;
    if(Q.size() == 2){
        join_result = new lqdag(FUNCTOR_AND, extend_qdag[0], extend_qdag[1]);
    }
    else if (Q.size() == 3){
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
        cout << "Code only implemented for a maximum of 5 qdags..." << endl;
        exit(1);
    }

    uint64_t p = std::pow(Q[0].getK(), A.size());

    // join result is the formula
    return join_result->completion(p, Q[0].getHeight()-1, 0, UPPER_BOUND, results);
}


quadtree_formula* compute_pred_lqdag_join(vector<qdag>&Q, uint64_t UPPER_BOUND, uint64_t &results, predicate* pred) {
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
    }

    lqdag* join_result;
    if(Q.size() == 2){
        join_result = new lqdag(FUNCTOR_AND, extend_qdag[0], extend_qdag[1]);
    }
    else if (Q.size() == 3){
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
        // TODO: fix this
        cout << "Code only implemented for a maximum of 5 qdags..." << endl;
        exit(1);
    }

    uint64_t grid_side = Q[0].getGridSide();
    uint8_t k = Q[0].getK();
    uint64_t p = std::pow(Q[0].getK(), A.size());
    uint16_t coordinates[A.size()];
    for(uint16_t i = 0; i < A.size(); i++)
        coordinates[i] = 0;

    // join result is the formula
    return join_result->completion_with_pred(p, Q[0].getHeight()-1, 0, grid_side, k, UPPER_BOUND, results, pred, coordinates, A.size());
}


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
    return diff->completion(p, q1.getHeight()-1, 0, UPPER_BOUND, results);
}


/**
 * Takes a subexpression F and a predicate θ, which is a logical expression on the attributes of F.
 */
lqdag* selection(lqdag* F, predicate* pred){
    // TODO
}

lqdag* projection(){
    // TODO
}

lqdag* semijoin(){
    // TODO
}

lqdag* antijoin(){
    // TODO
}

lqdag* division(){
    // TODO
}

