//
// We will not use qdags, only quadtrees
//
#ifndef INCLUDED_LQDAGS
#define INCLUDED_LQDAGS

#include<bits/stdc++.h>
#include "./qdags.hpp"
#include "./utils.hpp"
#include "./lqdags/utils_lqdags.hpp"
#include "./lqdags/formula_lqdag.hpp"


// lqdag as a syntax tree
class lqdag {

public:
	typedef vector<position> position_qdags;


private:
	qdag** arr_qdags; // array of qdags
    formula_lqdag *form_lqdag; // expression, syntax tree, formula of the lqdag (the same for its children)
	position_qdags* pos_qdags; // each qdag will have a space in this vector with its position (level and coordinates)
	node_completion* node_completion_lqdag; // children of the lqdag (computed)

public:
    lqdag() = default;

	lqdag(qdag** arr_qdags, formula_lqdag *form_lqdag) {
		this->arr_qdags = arr_qdags;
		this->form_lqdag = form_lqdag;
		this->pos_qdags = new position_qdags{};
		this->node_completion_lqdag = new node_completion{};
	}

    lqdag(qdag** arr_qdags, formula_lqdag *form_lqdag, position_qdags* pos_qdags) {
		this->arr_qdags = arr_qdags;
		this->form_lqdag = form_lqdag;
		this->pos_qdags = pos_qdags;
		this->node_completion_lqdag = new node_completion{};
	}

	lqdag(qdag** arr_qdags, formula_lqdag *form_lqdag, position_qdags* pos_qdags, node_completion* completion_children_lqdag) {
		this->arr_qdags = arr_qdags;
		this->form_lqdag = form_lqdag;
		this->pos_qdags = pos_qdags;
		this->node_completion_lqdag = completion_children_lqdag;
	}

	uint8_t get_functor(){
		return this->form_lqdag->get_functor();
	}

	formula_lqdag* get_formula(){
		return this->form_lqdag;
	}

	double val_lqdag1(){
		return this->form_lqdag->get_val_lqdag1();
	}

	double val_lqdag2(){
		return this->form_lqdag->get_val_lqdag2();
	}

	void set_val_lqdag1(double newVal){
		this->form_lqdag->set_val_lqdag1(newVal);
	}

	void set_val_lqdag2(double newVal){
		this->form_lqdag->set_val_lqdag2(newVal);
	}

	uint8_t getK() const {
		return this->form_lqdag->getK();
	}

	uint16_t nAttr() const {
		return this->form_lqdag->nAttr();
	}

	uint16_t getMaxLevel() const {
		return this->form_lqdag->getMaxLevel();
	}

	uint64_t getGridSide() const {
		return this->form_lqdag->getGridSide();
	}

	uint64_t n_qdags(){
		return this->pos_qdags->size(); // number of qdags
	}

	/**
	 * Compute the value of the i_qdag-th qdag to the current level and current node of the i-th qdag.
	 * And it returns true if the curr_node has a child in the i-th position.
	 * @param i_qdag the position in the pos_qdags vector.
	 * @param i_child the i-th child of the current node.
	 * @param pos_i_child the position of the i-th child in the next level.
	 * Return value of the i-th child: FULL_LEAF, EMPTY_LEAF, NO_VALUE_LEAF
	 */
	double val_child_qdag(uint64_t i_qdag, uint64_t i_child, uint64_t& pos_i_child){
		uint16_t cur_level = (*this->pos_qdags)[i_qdag].level;
		uint64_t cur_node = (*this->pos_qdags)[i_qdag].cur_node;
		uint16_t max_level = this->form_lqdag->getMaxLevel();
		bool node_exists = true;
		if(cur_level == max_level) {
			pos_i_child = cur_node + i_child;
			node_exists = this->arr_qdags[i_qdag]->get_ith_bit(cur_level, pos_i_child);
			return node_exists ? FULL_LEAF : EMPTY_LEAF;
		} else{
			// ith_child: start position of the ith child of cur_node in the next level (cur_level + 1)
			pos_i_child = this->arr_qdags[i_qdag]->get_child(cur_level, cur_node, i_child, node_exists);
			return node_exists ? NO_VALUE_LEAF : EMPTY_LEAF; // if node exists, we will have to compute its value after
		}
	}


	/**
	 * Get the new positions for each qdags of its children.
	 * @param i the i-th children
	 * @return
	 */
	position_qdags* get_child_pos_qdags(uint16_t i){
		// vector of the children positions for each branch (qdags at the end, leaves)
		vector<uint16_t> cur_pos_nodes_qdags; // cur_pos_nodes_qdags[j] is the i-th child of the j-th qdag
		this->form_lqdag->get_index_i(i, cur_pos_nodes_qdags); // get the index i that corresponds to each qdag
		// new vector with the children coordinates and level
		position_qdags* new_child_pos_nodes_qdags;
		for(uint64_t j = 0; j < this->pos_qdags->size(); j++){ // for each qdag
			// copy this->pos_qdags[j].coordinates
			uint16_t* copy_curr_coordinates = new uint16_t[this->nAttr()];
			uint16_t diff_level = getMaxLevel()-(*pos_qdags)[j].level;
			for(uint16_t k = 0; k < this->nAttr(); k++)
				copy_curr_coordinates[k] = (*pos_qdags)[j].coordinates[k];
			// compute the children coordinates path
			transformCoordinates(copy_curr_coordinates, nAttr(), diff_level, cur_pos_nodes_qdags[j]);
			new_child_pos_nodes_qdags->push_back(position((*pos_qdags)[j].level, copy_curr_coordinates));
		}
		return new_child_pos_nodes_qdags;
	}

	// TODO: check child is a uint16_t and not uint64_t
	/**
	 * Return the child of the lqdag.
	 * This will have the same formula as the parent, but with a different position and children.
	 * We assume that the parent is not empty
	 * @param i
	 * @return
	 */
	 // MALOOO debo ver si es OR
	lqdag* compute_child_lqdag(uint16_t i){
		bool child_exist = true;
		// get index of each child: cur_pos_nodes_qdags has n_qdags() elements
		vector<uint16_t> cur_pos_nodes_qdags(n_qdags()); // cur_pos_nodes_qdags[j] is the i-th child of the j-th qdag
		this->form_lqdag->get_index_i(i, cur_pos_nodes_qdags);
		// check value of position
		// compute the i-th child for each qdag to check if there is a possible child.
		vector<uint64_t> pos_ith_child(n_qdags());
		double val_qdag_j;
		for(uint64_t j = 0; j < n_qdags(); j++){
			val_qdag_j = this->val_child_qdag(j, cur_pos_nodes_qdags[j], pos_ith_child[j]);
			if(val_qdag_j == EMPTY_LEAF){
				child_exist = false;
				break;
			}
		}

		if(!child_exist){
			cur_pos_nodes_qdags.clear();
			pos_ith_child.clear();
			return new lqdag(this->arr_qdags, this->form_lqdag, nullptr, new node_completion{EMPTY_LEAF,0});
		}

		// call to transform coordinates with the position (level, coordinates
		position_qdags* new_pos_qdags = this->get_child_pos_qdags(i);
		node_completion* new_completion_children_lqdag = new node_completion(NO_VALUE_LEAF, this->node_completion_lqdag->n_children);
		lqdag* child_lqdag = new lqdag(this->arr_qdags, this->form_lqdag, new_pos_qdags, new_completion_children_lqdag);
		return child_lqdag;
	}


    // algorithm 8 child(L,i)
    /**
     * We descend recursively by the syntax tree (formula_lqdag). Once we arrive to the leaves we compute the pos_qdags[j]
     * Can only be invoked when value is 1/2 or VALUE_NEED_CHILDREN.
     * The formula will remain the same (we use a pointer to the original formula)
     * @param i
     * @return
     */
    lqdag* get_child_lqdag(formula_lqdag* formula_pos, uint64_t ith_child, uint64_t jth_qdag, position_qdags& m_pos_qdag, double& val_qdag){
        // case base
        switch (this->get_functor()) {
            case FUNCTOR_QTREE: case FUNCTOR_NOT: {
				// computer los valores de pos_qdags y node_copletion y
				uint64_t pos_ith_child = 0;
				val_qdag = val_child_qdag(jth_qdag,ith_child,pos_ith_child);
				m_pos_qdag[jth_qdag].level = (*pos_qdags)[jth_qdag].level+1;
				m_pos_qdag[jth_qdag].cur_node = pos_ith_child;
				uint16_t* copy_curr_coordinates = new uint16_t[this->nAttr()];
				for(uint16_t k = 0; k < this->nAttr(); k++)
					copy_curr_coordinates[k] = (*pos_qdags)[jth_qdag].coordinates[k];
				uint16_t diff_level = getMaxLevel()-(*pos_qdags)[jth_qdag].level;
				// compute the children coordinates path
				transformCoordinates(copy_curr_coordinates, nAttr(), diff_level, ith_child);
				m_pos_qdag[jth_qdag].coordinates = copy_curr_coordinates;
				return new lqdag(this->arr_qdags,
								 this->form_lqdag,
								 &m_pos_qdag,
								 new node_completion{ val_qdag,
													  val_qdag == EMPTY_LEAF ? 0 : this->node_completion_lqdag->n_children });

            }
            case FUNCTOR_AND: {
                double val_l1 = this->val_lqdag1() != NO_VALUE_LEAF ? this->val_lqdag1() : this->form_lqdag->get_formula_lqdag1()->value_lqdag();
				position* pos_qdag_and = new position{};
				if (val_l1== 1)
                    return get_child_lqdag(this->form_lqdag->get_formula_lqdag2(), ith_child, jth_qdag+1, pos_qdag_and, val_qdag);
                double val_l2 = this->val_lqdag2() != NO_VALUE_LEAF ? this->val_lqdag2() : this->form_lqdag->get_formula_lqdag2()->value_lqdag();
                if (val_l2 == 1)
					return get_child_lqdag(this->form_lqdag->get_formula_lqdag1(), ith_child, jth_qdag);
//                lqdag *l = new lqdag(this->arr_qdags, this->form_lqdag, get_child_lqdag(this->form_lqdag->get_formula_lqdag1(), i, count_qdag), get_child_lqdag(this->form_lqdag->get_formula_lqdag2(), i, count_qdag+1));
                return l;
            }
            case FUNCTOR_OR: {
                double val_l1 = this->val_lqdag1() != NO_VALUE_LEAF ? this->val_lqdag1() : this->lqdag1->value_lqdag();
                if (val_l1 == 0)
                    return lqdag2->get_child_lqdag(i);
                double val_l2 = this->val_lqdag2 != NO_VALUE_LEAF ? this->val_lqdag2 : this->lqdag2->value_lqdag();
                if (val_l2 == 0)
                    return lqdag1->get_child_lqdag(i);
                lqdag* l = new lqdag(this->form_lqdag, this->lqdag1->get_child_lqdag(i), this->lqdag2->get_child_lqdag(i));
                return l;
            }
            case FUNCTOR_EXTEND: { // Extend quadtree to attributes att_set
                // i --> i' (hago el mapeo)
                return new lqdag(this->form_lqdag, this->get_child_pos_qdags(M[i]));
            }
            default:
                throw "error: value_lqdag non valid functor";

        }
    }

    /**
     * algorithm 7 value(L).
     * In lqdags we introduce a new idea: full leaves, that denote subgrids full of 1s.
     * We save the value in val_lqdag_1, val_lqdag_2 or subQuatree->value. If we've already computed the value, we return it.
     * @return the value of the root of the lqdag. Can be EMPTY_LEAF, FULL_LEAF, INTERNAL_NODE or VALUE_NEED_CHILDREN.
     */
    double value_lqdag(){
        switch(this->get_functor()){
            case FUNCTOR_QTREE: {
                return this->value_quadtree();
            }
            case FUNCTOR_NOT: {
                return 1 - this->value_quadtree();
            }
            case FUNCTOR_AND: {
                if(this->val_lqdag1 == NO_VALUE_LEAF)
                    this->val_lqdag1 = this->lqdag1->value_lqdag();
                if(this->val_lqdag1 == 0)
                    return EMPTY_LEAF;

                if(this->val_lqdag2 == NO_VALUE_LEAF)
                    this->val_lqdag2 = this->lqdag2->value_lqdag();
                if(this->val_lqdag2 == 0)
                    return EMPTY_LEAF;

                if (this->val_lqdag1 == 1)
                    return this->val_lqdag2;
                if (this->val_lqdag2 == 1)
                    return this->val_lqdag1;

                return VALUE_NEED_CHILDREN; // if both roots have Value 1⁄2, one cannot be sure of the Value of the resulting root until the AND between the children of Q1 and Q2 has been computed
            }
            case FUNCTOR_OR: {
                if(this->val_lqdag1 == NO_VALUE_LEAF)
                    this->val_lqdag1 = this->lqdag1->value_lqdag();
                if(this->val_lqdag1 == 1)
                    return FULL_LEAF;

                if(this->val_lqdag2 == NO_VALUE_LEAF)
                    this->val_lqdag2 = this->lqdag2->value_lqdag();
                if(this->val_lqdag2 == 1)
                    return FULL_LEAF;

                if (this->val_lqdag1 == 0)
                    return this->val_lqdag2;
                if (this->val_lqdag2 == 0)
                    return this->val_lqdag1;

                return VALUE_NEED_CHILDREN; // if both roots have Value 1⁄2, one cannot be sure of the Value of the resulting root until the OR between the children of Q1 and Q2 has been computed
            }
            case FUNCTOR_EXTEND: {
                if (this->val_lqdag1 == NO_VALUE_LEAF)
                    this->val_lqdag1 = this->lqdag1->value_lqdag();
                return this->val_lqdag1;
            }
            default:
                throw "error: value_lqdag non valid functor";
        }
    }

    /**
     * Algorithm 1 Value(Q), adapted for lqdag
     * Value of a qdag
     * In lqdags we introduce a new idea: full leaves, that denote subgrids full of 1s.
     * Leaves can be in a higher level than the last one, when the subgrids are all 0s or all 1s.
     * FULL_LEAF if the qdag represents a full single cell, EMPTY_LEAF if it is an empty leaf, INTERNAL_NODE if is an internal quadtree_formula.
     * @return EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE.
     */
    double value_qdag(){
        if(this->get_functor() == FUNCTOR_QTREE) { // only is called
            uint16_t cur_level = this->pos_qdags->level;
            uint64_t cur_node = this->subQuadtree->node;
            if(this->subQuadtree->value == NO_VALUE_LEAF)
                this->subQuadtree->value = is_a_qdag_leaf(cur_level, cur_node);
            return this->subQuadtree->value;

        } else {
            //cout << "error: value_quadtree called on non-quadtree" << endl;
            throw "error: value_quadtree called on non-quadtree";
        }
    }

    /**
     * This function is called from a (FUNCTOR_QTREE, Q_f). Detect an emptu leaf, full leaf or an internal node.
     * @param level
     * @param node
     * @param k
     * @param k_d
     * If level = 0 and node = 0, we are looking the root of the quadtree.
     * If level > max_level, we are looking for the bit representing a point in the grid.
     * @return whether is an EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE
     */
    double is_a_qdag_leaf(uint16_t level, uint64_t node){
        uint16_t max_level = this->form_lqdag->getMaxLevel();
        if(level > max_level){
            return this->form_lqdag->get_qdag()->get_ith_bit(max_level, node);
        }
        uint8_t node_description = this->form_lqdag->get_qdag()->get_node_last_level(level, node);
        if((node_description | 0) == 0)
            return EMPTY_LEAF;
        uint16_t val_full_ones;
        uint64_t k_d = this->form_lqdag->get_qdag()->getKD();
        switch (k_d) {
            case 2:
                val_full_ones = 3;
                break;
            case 4:
                val_full_ones = 15;
                break;
            default:
                throw "error: invalid k_d";
        }
        if((node_description & val_full_ones) == val_full_ones){ // check the children
            if(level == max_level)
                return FULL_LEAF;
            uint64_t siblings = this->form_lqdag->get_qdag()->rank(level,node);
            for(uint64_t i = 0; i < k_d; i++){
                if(is_a_qdag_leaf(level+1, (siblings+i)*k_d) != 1)
                    return INTERNAL_NODE;
            }
            return FULL_LEAF;
        }
        return INTERNAL_NODE;
    }


    /**
     * Algorithm 9 of paper.
     * The completion Q_f is the quadtree representing the output of the formula F, represented as an lqdag.
     *
     * @param p
     * @param max_level
     * @param cur_level the current level of the quadtree. If cur_level = max_level + 1, we are looking for the bit representing the point in the grid.
     * @param UPPER_BOUND
     * @param results
     * @return
     */
    quadtree_formula* completion(uint64_t p,
                                 uint16_t max_level,
                                 uint16_t cur_level,
                                 uint64_t UPPER_BOUND,
                                 uint64_t &results) {

        if(results >= UPPER_BOUND){ // top k results only
            quadtree_formula* newNode = create_leaf(NO_VALUE_LEAF);
            return newNode;
        }
        double val_lqdag = this->value_lqdag();
        if(val_lqdag == EMPTY_LEAF || val_lqdag == FULL_LEAF){
            quadtree_formula* newNode = create_leaf(val_lqdag);
            if(val_lqdag == FULL_LEAF){
                results += (cur_level > max_level) ? 1 :  pow(p, (max_level - cur_level + 1));
            }
            return newNode;
        }

        double max_value = 0;
        double min_value = 1;
        quadtree_formula **Q_f_children = new quadtree_formula*[p];
        // if all quadtres are empty, return a 0 leaf.
        // C_i ← Completion(Child(L_F,i))
        for(uint64_t i = 0; i < p; i++){
            Q_f_children[i] = this->get_child_lqdag(i)->completion(p, max_level, cur_level+1, UPPER_BOUND, results);
            max_value = max(max_value, Q_f_children[i]->val_leaf); // val_leaf can be EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE
            min_value = min(min_value, Q_f_children[i]->val_leaf);
        }
        // return a leaf
        if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){ // return a leaf
            for(uint64_t i = 0; i < p; i++)
                delete Q_f_children[i];
            delete[] Q_f_children;
            double val_leaf = (max_value == EMPTY_LEAF) ? EMPTY_LEAF : FULL_LEAF;
            quadtree_formula* newNode = create_leaf(val_leaf);
            return newNode;
        }
        else {
            // return an internal node with children
            quadtree_formula* quadtree = new quadtree_formula();
            quadtree->val_leaf = val_lqdag; //INTERNAL_NODE;
            quadtree->children = Q_f_children;
            return quadtree;
        }
    }


    /**
     *
     * @param p
     * @param max_level
     * @param cur_level
     * @param UPPER_BOUND
     * @param results
     * @param pred
     * @param coordinates an array of nAtt coordinates of the node representing the quadrant.
     * @param checkPred once we have evaluated the predicate and if this evaluation returns 1, we do not have to evaluate it again for its children.
     * @return
     */
    quadtree_formula *
    completion_with_pred(uint64_t p, uint16_t max_level, uint16_t cur_level, uint64_t grid_side, uint8_t k, uint64_t UPPER_BOUND, uint64_t &results,
                         predicate *pred, uint16_t *coordinates, uint64_t nAttr) {

        if(results >= UPPER_BOUND){ // top k results only
            quadtree_formula* newNode = create_leaf(EMPTY_LEAF);
            return newNode;
        }
        // if results < UPPER_BOUND
        uint64_t quadrant_side = grid_side/pow(k,cur_level);
        double val_eval_pred = eval_pred(pred, coordinates, quadrant_side, nAttr);
        if(val_eval_pred == 0){
            quadtree_formula* newNode = create_leaf(EMPTY_LEAF);
            return newNode;
        } else if(val_eval_pred == 1){
            return this->completion(p, max_level, cur_level, UPPER_BOUND, results);
        } else {
            double max_value = 0;
            double min_value = 1;
            quadtree_formula **Q_f_children = new quadtree_formula*[p];
            // if all quadtres are empty, return a 0 leaf.
            // C_i ← Completion(Child(L_F,i))
            for(uint64_t i = 0; i < p; i++){
                uint16_t diff_level = max_level-cur_level;
                uint16_t l = (uint16_t) log2(p);
                uint16_t* coordinatesTemp = new uint16_t[l];
                for(uint16_t j = 0; j < l; j++)
                    coordinatesTemp[j] = coordinates[j];
                transformCoordinates(coordinatesTemp, l, diff_level, i);
                quadtree_formula* newNode = this->get_child_lqdag(i)->completion_with_pred(p, max_level, cur_level + 1,
                                                                                           grid_side, k,
                                                                                           UPPER_BOUND, results, pred,
                                                                                           coordinatesTemp, nAttr);
                Q_f_children[i] = newNode;
                max_value = max(max_value, newNode->val_leaf); // val_leaf can be EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE
                min_value = min(min_value, newNode->val_leaf);
            }
            // return a leaf
            if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){ // return a leaf
                for(uint64_t i = 0; i < p; i++)
                    delete Q_f_children[i];
                delete[] Q_f_children;
                double val_leaf = (max_value == EMPTY_LEAF)? EMPTY_LEAF : FULL_LEAF;
                quadtree_formula* newNode = create_leaf(val_leaf);
                return newNode;
            }
            else {
                // return an internal node with children
                quadtree_formula* quadtree = new quadtree_formula;
                quadtree->val_leaf = INTERNAL_NODE;
                quadtree->children = Q_f_children;
                return quadtree;
            }

        }

    }


    /**
     * It builds the quadtree_formula Q_f of the lqdag L_F, as we ask for a particular node.
     * @param p
     * @param max_level
     * @param cur_level
     * @param UPPER_BOUND
     * @param results
     * @param level
     * @param child the i-th child of the current node
     * @param Q_f the quadtree_formula built so far.
     * @return a pointer to the child of the lqdag
     */
    lqdag* lazy_child_completion(uint64_t p,
                                   uint64_t child,
                                   quadtree_formula& Q_f){
        if(Q_f.val_leaf == EMPTY_LEAF || Q_f.val_leaf == FULL_LEAF){
            return this; // child is 0 or 1 and the quadtree is already computed
        }
        else if(Q_f.val_leaf == INTERNAL_NODE || Q_f.val_leaf == VALUE_NEED_CHILDREN){
            lqdag* lqdag_child = this->get_child_lqdag(child);
            if(!Q_f.children[child]){ // allocate memory for the child
                Q_f.children[child] = new quadtree_formula{};
            } // if it's not computed
            if(Q_f.children[child]->val_leaf == NO_VALUE_LEAF){
                quadtree_formula* qf_no_value = lqdag_child->compute_quadtree_formula_node(p);
                Q_f.children[child]->val_leaf = qf_no_value->val_leaf;
                Q_f.children[child]->children = qf_no_value->children;
            } // else: it's already computed
            return lqdag_child;
        }
        else{ // NO_VALUE --> compute the root and then, the child
            quadtree_formula* qf_no_value = this->compute_quadtree_formula_node(p);
            Q_f.val_leaf = qf_no_value->val_leaf;
            Q_f.children = qf_no_value->children;
            return this->lazy_child_completion(p, child, Q_f);
        }
    }

	/**
	 * Lazy completion with the selection operation in relational algebra
	 * @param p
	 * @param child
	 * @param Q_f
	 * @param qf_pred
	 * @return
	 */
	output_lqdag_pred* lazy_child_completion_with_pred(uint64_t p,
                                           uint64_t child,
                                           quadtree_formula& Q_f,
                                           quadtree_pred* qf_pred){
        if(Q_f.val_leaf == EMPTY_LEAF || Q_f.val_leaf == FULL_LEAF){ // no need to calculate its child
            return new output_lqdag_pred{this, qf_pred};
        }
        else if(Q_f.val_leaf == INTERNAL_NODE || Q_f.val_leaf == VALUE_NEED_CHILDREN){
            // calculate coordinates
            uint16_t diff_level = qf_pred->max_level-qf_pred->cur_level;
            uint16_t l = (uint16_t) log2(p);
            uint16_t* coordinatesTemp = new uint16_t[l];
            for(uint16_t j = 0; j < l; j++)
                coordinatesTemp[j] = qf_pred->coordinates[j];
            transformCoordinates(coordinatesTemp, l, diff_level, child);
            qf_pred->coordinates = coordinatesTemp;
            qf_pred->cur_level += 1;
            // value predicate
            uint64_t quadrant_side = qf_pred->grid_side/pow(qf_pred->k,qf_pred->cur_level);
            double val_eval_pred = eval_pred(qf_pred->pred, qf_pred->coordinates, quadrant_side, qf_pred->nAttr);

            if(val_eval_pred == 0.5  || val_eval_pred == 1){ // if it passes the predicate, we compute the child
                lqdag* lqdag_child = this->get_child_lqdag(child);
                if(!Q_f.children[child]){ // allocate memory for the child
                    Q_f.children[child] = new quadtree_formula{};
                } // if it's not computed
                if(Q_f.children[child]->val_leaf == NO_VALUE_LEAF){
                    quadtree_formula* qf_no_value = lqdag_child->compute_quadtree_formula_node(p);
                    Q_f.children[child]->val_leaf = qf_no_value->val_leaf;
                    Q_f.children[child]->children = qf_no_value->children;
                } // else: it's already computed
//				return lqdag_child;
                return new output_lqdag_pred{lqdag_child, qf_pred};
            } else{
                Q_f.children[child] = create_leaf(EMPTY_LEAF);
				return nullptr;
                return new output_lqdag_pred{nullptr, qf_pred};
            }


        }
        else{ // NO_VALUE --> compute the root and then, the child
            quadtree_formula* qf_no_value = this->compute_quadtree_formula_node(p);
            Q_f.val_leaf = qf_no_value->val_leaf;
            Q_f.children = qf_no_value->children;
//			return this->lazy_child_completion_with_pred(p, child, Q_f, qf_pred);
            return new output_lqdag_pred{this->lazy_child_completion(p, child, Q_f), qf_pred};
        }

    }

};

#endif