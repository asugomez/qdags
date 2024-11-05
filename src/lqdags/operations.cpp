//
// Created by Asunción Gómez on 19-04-24.
//

#include "../lqdag.hpp"
#include "formula_lqdag.hpp"

// ---------- QUADTREE  ---------- //

formula_lqdag* create_formula_qdag_leaf(qdag* Q){
	return new formula_lqdag(FUNCTOR_QTREE, Q);
}

formula_lqdag create_formula_lqdag_leaf(lqdag* L){
	return formula_lqdag(FUNCTOR_LQDAG, L);
}


// ---------- JOIN ---------- //

// (JOIN, L1(A1), L2(A2)) = (AND, (EXTEND, L1, A1 ∪ A2), (EXTEND, L2, A1 ∪ A2))
/**
 * Formula for the join operation
 * @param Q vector with the qdags to join
 * @return
 */
lqdag* lqdag_join(vector<qdag> &Q/* uint64_t &p*/) {
    qdag::att_set A;
    map<uint64_t, uint8_t> attr_map;

    for (uint64_t i = 0; i < Q.size(); i++) {
        uint64_t nAttr = Q.at(i).nAttr();
        for (uint64_t j = 0; j < nAttr; j++)
            attr_map[Q.at(i).getAttr(j)] = 1;
    }

    // set of attributs
    for (map<uint64_t, uint8_t>::iterator it = attr_map.begin(); it != attr_map.end(); it++)
        A.push_back(it->first);

//	uint16_t p = std::pow(Q.at(0).getK(), A.size());

    uint8_t k = Q.at(0).getK();
	uint16_t max_level = Q.at(0).getHeight()-1;
	uint64_t grid_side = Q.at(0).getGridSide();


    formula_lqdag* extend_formula_lqdag[Q.size()];

    for (uint64_t i = 0; i < Q.size(); i++) { // EXTEND((Q_TREE, index) , A)
        extend_formula_lqdag[i] = new formula_lqdag(FUNCTOR_EXTEND, new formula_lqdag(&Q.at(i), i),k, A, max_level, grid_side);
    }

	formula_lqdag* join_formula;

    if(Q.size() == 2){
		join_formula = new formula_lqdag(FUNCTOR_AND, extend_formula_lqdag[0], extend_formula_lqdag[1], k, A.size(), max_level, grid_side);
    }
    else if (Q.size() == 3){
		formula_lqdag* join_r_s = new formula_lqdag(FUNCTOR_AND, extend_formula_lqdag[0], extend_formula_lqdag[1]);
		join_formula = new formula_lqdag(FUNCTOR_AND, join_r_s, extend_formula_lqdag[2], k, A.size(), max_level, grid_side);
    }
    else if (Q.size() == 4){
		formula_lqdag* join_r_s = new formula_lqdag(FUNCTOR_AND, extend_formula_lqdag[0], extend_formula_lqdag[1]);
		formula_lqdag* join_t_u = new formula_lqdag(FUNCTOR_AND, extend_formula_lqdag[2], extend_formula_lqdag[3]);
		join_formula = new formula_lqdag(FUNCTOR_AND, join_r_s, join_t_u, k, A.size(), max_level, grid_side);
    }
    else if (Q.size() == 5){
		formula_lqdag* join_r_s = new formula_lqdag(FUNCTOR_AND, extend_formula_lqdag[0], extend_formula_lqdag[1]);
		formula_lqdag* join_t_u = new formula_lqdag(FUNCTOR_AND, extend_formula_lqdag[2], extend_formula_lqdag[3]);
		formula_lqdag* join_r_s_t_u = new formula_lqdag(FUNCTOR_AND, join_r_s, join_t_u);
		join_formula = new formula_lqdag(FUNCTOR_AND, join_r_s_t_u, extend_formula_lqdag[4], k, A.size(), max_level, grid_side);
    }
    else {
        cout << "Code only implemented for a maximum of 5 qdags..." << endl;
        exit(1);
    }


    return new lqdag{Q, join_formula};
}

lqdag* compute_dfs(lqdag* l, uint64_t UPPER_BOUND, uint64_t &results){
	return l->completion_dfs(l->getMaxLevel(),0,UPPER_BOUND,results);
}

lqdag* compute_dfs_join(vector<qdag> &Q, uint64_t UPPER_BOUND, uint64_t &results){
	lqdag* join = lqdag_join(Q);
	return compute_dfs(join,UPPER_BOUND,results);
}

lqdag* compute_dfs_join_run_time(lqdag* join,uint64_t UPPER_BOUND, uint64_t &results){
	return compute_dfs(join,UPPER_BOUND,results);
}

lqdag* compute_dfs_selection(lqdag* l, predicate* pred, uint64_t UPPER_BOUND, uint64_t &results){
	lqdag* copy_lqdag = new lqdag(l->get_arr_qdags(), l->get_formula(), l->get_position_qdags(), l->get_node_completion_lqdag());
	// TODO: OFT
	copy_lqdag->get_node_completion_lqdag()->val_node = 5;
	// TODO: check original lqdag l didn't change
	return copy_lqdag->completion_selection_dfs(pred,l->getMaxLevel(),0,UPPER_BOUND,results);
}

lqdag* compute_dfs_selection_join(vector<qdag> &Q, predicate* pred, uint64_t UPPER_BOUND, uint64_t &results){
	lqdag* join = lqdag_join(Q);
	return join->completion_selection_dfs(pred,join->getMaxLevel(),0,UPPER_BOUND,results);
}



//lqdag* projection(){
//    // TODO
//}
//
//lqdag* semijoin(){
//    // TODO
//}
//
//lqdag* antijoin(){
//    // TODO
//}
//
//lqdag* division(){
//    // TODO
//}

//// (DIFF, L1(A), L2(A)) = (AND, L1, (NOT, L2))
//lqdag* diff(qdag q1, qdag q2, bool bounded_result, uint64_t UPPER_BOUND) {
//	subQuadtreeChild* subQ1 = new subQuadtreeChild{&q1, 0, 0};
//	subQuadtreeChild* subQ2 = new subQuadtreeChild{&q2, 0, 0};
//	lqdag* diff = new lqdag(FUNCTOR_AND, new lqdag(FUNCTOR_QTREE, subQ1), new lqdag(FUNCTOR_NOT, subQ2));
//	return diff;
//}
//
//quadtree_formula* compute_lqdag_diff(qdag q1, qdag q2, bool bounded_result, uint64_t UPPER_BOUND,uint64_t &results){
//	lqdag* diff_res = diff(q1, q2, bounded_result, UPPER_BOUND);
//	uint64_t p = std::pow(q1.getK(), q1.nAttr());
//	return diff_res->completion(p, q1.getHeight()-1, 0, UPPER_BOUND, results);
//}

