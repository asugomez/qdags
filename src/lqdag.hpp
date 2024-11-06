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
//	typedef vector<position> position_qdags;

	/**
	 * this is the materialization: it has the value of the node, the level, the coordinates, the number of children and the lqdag children.

	 */
	struct node_completion{
		double val_node = NO_VALUE_LEAF; // EMPTY_LEAF, FULL_LEAF, INTERNAL_NODE, NO_VALUE_LEAF
		uint16_t level;
		uint64_t n_children; // number of children
		uint256_t path;
		lqdag** lqdag_children = nullptr;
		// default constructor
		node_completion() : val_node(NO_VALUE_LEAF), level(0),  path(0),n_children(0), lqdag_children(nullptr) {}
		// constructor with p children
		node_completion(uint64_t p) : val_node(NO_VALUE_LEAF), level(0), path(0), n_children(p), lqdag_children(new lqdag*[p]) {
			for(uint64_t i = 0; i < p; i++)
				lqdag_children[i] = nullptr;
		}
		// constructor with value and p children
		node_completion(double val, uint64_t p) : val_node(val),  level(0), path(0), n_children(p), lqdag_children(new lqdag*[p]) {
			for(uint64_t i = 0; i < p; i++)
				lqdag_children[i] = nullptr;
		}
		// copy assignment operator (deep copy)
		node_completion& operator=(const node_completion& n){
			// TODO: see when we call this copy
			if(this != &n){
				// Clean up existing resources
				if(this->lqdag_children != nullptr)
					delete[] this->lqdag_children;

				this->val_node = n.val_node;
				this->level = n.level;
				this->n_children = n.n_children;
				this->path = n.path;

				// Allocate and copy coordinates
				this->lqdag_children = new lqdag*[n.n_children];
				for(uint64_t i = 0; i < n.n_children; i++)
//					this->lqdag_children[i] = n.lqdag_children[i];
					// TODO: create an equal operator!
					this->lqdag_children[i] = new lqdag(*n.lqdag_children[i]);
			}
			return *this;
		}
		// destructor
		~node_completion(){
			if(lqdag_children != nullptr){
				delete[] lqdag_children;
			}
		}
	};


private:
	vector<qdag> arr_qdags; // array of qdags
    formula_lqdag *form_lqdag; // expression, syntax tree, formula of the lqdag (the same for its children)
	vector<position> pos_qdags; // each qdag will have a space in this vector with its position (level and coordinates)
	node_completion* node_completion_lqdag; // children of the lqdag (computed)


public:
    lqdag() = default;

	/**
	 *
	 * @param arr_qdags array of qdags (the leaves)
	 * @param form_lqdag syntax tree
	 * @param p k^d number of children
	 */
	lqdag(vector<qdag> arr_qdags, formula_lqdag *form_lqdag) {
		this->arr_qdags = arr_qdags;
		this->form_lqdag = form_lqdag;
		this->pos_qdags =  vector<position>(arr_qdags.size());
		uint16_t p = std::pow(getK(), nAttr()); // TODO: test pos_qdags
		this->node_completion_lqdag = new node_completion(p);
	}

    lqdag(vector<qdag> arr_qdags, formula_lqdag *form_lqdag, vector<position> pos_qdags) {
		this->arr_qdags = arr_qdags;
		this->form_lqdag = form_lqdag;
		this->pos_qdags = pos_qdags;
		uint16_t p = std::pow(getK(), nAttr());
		this->node_completion_lqdag = new node_completion{NO_VALUE_LEAF, p};
	}

	lqdag(vector<qdag> arr_qdags, formula_lqdag *form_lqdag, vector<position> pos_qdags, node_completion* completion_children_lqdag) {
		this->arr_qdags = arr_qdags;
		this->form_lqdag = form_lqdag;
		this->pos_qdags = pos_qdags;
		this->node_completion_lqdag = completion_children_lqdag;
	}

	~lqdag(){
		// TODO: finish this
//		delete pos_qdags;
		delete node_completion_lqdag;
	}

	// ------------------ GETTERS ------------------ //

	vector<qdag> get_arr_qdags(){
		return this->arr_qdags;
	}

	uint64_t get_index(){
		return this->form_lqdag->get_index();
	}

	uint8_t get_functor(){
		return this->form_lqdag->get_functor();
	}

	formula_lqdag* get_formula(){
		return this->form_lqdag;
	}

	vector<position> get_position_qdags(){
		return this->pos_qdags;
	}

	node_completion* get_node_completion_lqdag(){
		return this->node_completion_lqdag;
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

	uint16_t getCurLevel() const { // TODO: check this, i think it's ok, because all qdags have the same grid side (domaine)
		return this->node_completion_lqdag->level;
	}

	uint64_t getGridSide() const {
		return this->form_lqdag->getGridSide();
	}

	uint64_t n_qdags(){
		return this->pos_qdags.size(); // number of qdags
	}

	uint64_t nChildren(){
		return this->node_completion_lqdag->n_children;
	}

	// ------------------ CHILD FUNCTION ------------------ //

	/**
	 * Compute the value of the i_qdag-th qdag to the current level and current node of the i-th qdag.
	 * And it returns true if the curr_node has a child in the i-th position.
	 * @param i_child the i-th child of the current node.
	 * @param pos_i_child the position of the i-th child in the next level.
	 * @param index_qdag the index of the current qdag
	 * Return value of the i-th child: FULL_LEAF, EMPTY_LEAF, NO_VALUE_LEAF
	 */
	double get_child_qdag(formula_lqdag* formula_pos, uint64_t index_qdag,uint64_t i_child, uint64_t& pos_i_child){
		if(index_qdag >= this->pos_qdags.size())
			cout << "error: get_child_qdag index out of range" << endl;
		uint16_t cur_level = (this->pos_qdags).at(index_qdag).level;
		uint64_t cur_node = (this->pos_qdags).at(index_qdag).cur_node;
		uint16_t max_level = formula_pos->getMaxLevel();
		bool node_exists = true;
		if(cur_level == max_level) {
			pos_i_child = cur_node + i_child;
			node_exists = this->arr_qdags.at(index_qdag).get_ith_bit(cur_level, pos_i_child);
			return node_exists ? FULL_LEAF : EMPTY_LEAF;
		} else{
			// ith_child: start position of the ith child of cur_node in the next level (cur_level + 1)
			pos_i_child = this->arr_qdags.at(index_qdag).get_child(cur_level, cur_node, i_child, node_exists);
			if(!node_exists)
				return EMPTY_LEAF;
			return is_a_qdag_leaf(formula_pos, index_qdag, cur_level+1, pos_i_child);
//			return node_exists ? NO_VALUE_LEAF : EMPTY_LEAF; // if node exists, we will have to compute its value after
		}
	}

	/**
	 * Update the coordiantes for a lqdag.position_qdags when we go down to a lqdag child
	 * @param ith_child
	 * @return
	 */
	position* get_position_data(uint64_t index_qdag, uint64_t ith_child, uint64_t curr_node){
		if(index_qdag >= this->pos_qdags.size())
			cout << "error: get_position_data index out of range" << endl;
		uint64_t nAttr_qdag = this->arr_qdags.at(index_qdag).nAttr();;
		uint16_t cur_level = (pos_qdags).at(index_qdag).level;
//		uint16_t diff_level = getMaxLevel() - cur_level;

		uint256_t cur_path = (pos_qdags).at(index_qdag).path;

		getNewMortonCodePath(cur_path, nAttr_qdag, cur_level, (uint256_t)ith_child);
		return new position{++cur_level, curr_node, cur_path}; // TODO: check curr_level ++
	}


	/**
	 * Return a new lqdag with the same array of qdags and the same formula_lqdag.
	 * We have to update the pos_qdags (level++, with a new node a new coordinates for each qdag)
	 * Assume we call this function when the current node is not an EMPTY_LEAF
	 * @param ith_child
	 * @return a new lqdag that corresponds to the ith-child of the current lqdag.
	 */
	lqdag* get_child_lqdag(uint64_t ith_child){
		// TODO: check empty and full, isnt it?
		if(this->node_completion_lqdag->val_node == EMPTY_LEAF || this->node_completion_lqdag->val_node == FULL_LEAF)
			return this;
		// check if the ith_child has been computed
		// TODO: check if get_child_lqdag is always called when != NO_VALUE_LEAF
		if (this->node_completion_lqdag->val_node != NO_VALUE_LEAF
				&& this->node_completion_lqdag->lqdag_children[0] != nullptr
				&& this->node_completion_lqdag->lqdag_children[ith_child] != nullptr) {
			return this->node_completion_lqdag->lqdag_children[ith_child];
		}

		// TODO: fix this! new_pos qdags is not stored, jere is the bug
		vector<position> new_pos_qdags;// = position_qdags[arr_qdags.size()];
		// copy the vector position of the current lqdag for the new lqdag
		if(arr_qdags.size()-1 >= this->pos_qdags.size())
			cout << "error: get_child_lqdag index out of range: arr_qdags.size" << endl;
		for(uint64_t i = 0; i < arr_qdags.size(); i++){
			try{
				new_pos_qdags.push_back(this->pos_qdags.at(i));
			} catch (const std::out_of_range& oor) {
				cout << "Out of Range error: " << i << " " << oor.what() << '\n';
			}
		}
//		new_pos_qdags[0] = pos_qdags[0]; // TODO see what happens

//		new_pos_qdags->at(0).cur_node = 40;

		double val_node = get_child_lqdag_aux(this->form_lqdag, ith_child, new_pos_qdags);

		lqdag* child_lqdag = new lqdag(this->arr_qdags, this->get_formula(), new_pos_qdags);
		// TODO: see if this->arr_qdags is a pointer (for vector) or not. It has to be a pointer!
		child_lqdag->node_completion_lqdag->val_node = val_node;
		child_lqdag->node_completion_lqdag->level = this->node_completion_lqdag->level +1;
		if(val_node != EMPTY_LEAF){
			uint256_t newPath = this->node_completion_lqdag->path;
			getNewMortonCodePath(newPath, nAttr(), this->getCurLevel(), (uint256_t)ith_child);
			child_lqdag->node_completion_lqdag->path = newPath;
		}

		this->node_completion_lqdag->lqdag_children[ith_child] = child_lqdag;
		return child_lqdag;

	}

    // algorithm 8 child(L,i)
    /**
     * We descend recursively by the syntax tree (formula_lqdag). Once we arrive to the leaves we compute the pos_qdags[j]
     * Can only be invoked when value is 1/2 or VALUE_NEED_CHILDREN.
     * The formula will remain the same (we use a pointer to the original formula)
     * @param i
     * @return the value of val_node.
     */
	double get_child_lqdag_aux(formula_lqdag* formula_pos, uint64_t ith_child, vector<position>& m_pos_qdag){
		// case base
		switch (formula_pos->get_functor()) {
			case FUNCTOR_QTREE: {
				uint64_t jth_child;
				double my_val  = get_child_qdag(formula_pos, formula_pos->get_index(), ith_child, jth_child);
				formula_pos->set_val_lqdag1(my_val);
				if(my_val!=EMPTY_LEAF)
					if(formula_pos->get_index() >= m_pos_qdag.size())
						cout << "error: get_child_lqdag_aux index out of range" << endl;
					m_pos_qdag.at(formula_pos->get_index()) = *get_position_data(formula_pos->get_index(), ith_child, jth_child);
				return my_val;
			}
			case FUNCTOR_NOT: {
				uint64_t jth_child;
				double my_val = get_child_qdag(formula_pos, formula_pos->get_index(), ith_child, jth_child);
				if(my_val!=NO_VALUE_LEAF) // apply NOT to value
					my_val = 1 - my_val;
				formula_pos->set_val_lqdag1(my_val);
				if(my_val!=EMPTY_LEAF)
					m_pos_qdag.at(formula_pos->get_index()) = *get_position_data(formula_pos->get_index(), ith_child, jth_child);
				return my_val;
			}
			case FUNCTOR_AND: {
				if(formula_pos->get_val_lqdag1() == EMPTY_LEAF || formula_pos->get_val_lqdag2() == EMPTY_LEAF)
					return EMPTY_LEAF;
				if(formula_pos->get_val_lqdag1() == FULL_LEAF){
					double val_l2 =  get_child_lqdag_aux(formula_pos->get_formula_lqdag2(), ith_child, m_pos_qdag);
					return val_l2;
				}
				if(formula_pos->get_val_lqdag2() == FULL_LEAF){
					double val_l1 = get_child_lqdag_aux(formula_pos->get_formula_lqdag1(), ith_child, m_pos_qdag);
					return val_l1;
				} else{
					double val_l1 = get_child_lqdag_aux(formula_pos->get_formula_lqdag1(), ith_child, m_pos_qdag);
					double val_l2 = get_child_lqdag_aux(formula_pos->get_formula_lqdag2(), ith_child, m_pos_qdag);
					if(!val_l1 || ! val_l2)
						return EMPTY_LEAF;
					else if(val_l1 == FULL_LEAF && val_l2 == FULL_LEAF)
						return FULL_LEAF;
					else
						return VALUE_NEED_CHILDREN;
				}
			}
			case FUNCTOR_OR: {
				if(formula_pos->get_val_lqdag1() == FULL_LEAF || formula_pos->get_val_lqdag2() == FULL_LEAF)
					return FULL_LEAF;
				if(formula_pos->get_val_lqdag1() == EMPTY_LEAF){
					double val_l2 =  get_child_lqdag_aux(formula_pos->get_formula_lqdag2(), ith_child, m_pos_qdag);
					return val_l2;
				}
				if(formula_pos->get_val_lqdag2() == EMPTY_LEAF){
					double val_l1 = get_child_lqdag_aux(formula_pos->get_formula_lqdag1(), ith_child, m_pos_qdag);
					return val_l1;
				} else{
					double val_l1 = get_child_lqdag_aux(formula_pos->get_formula_lqdag1(), ith_child, m_pos_qdag);
					double val_l2 = get_child_lqdag_aux(formula_pos->get_formula_lqdag2(), ith_child, m_pos_qdag);
					if(val_l1 == FULL_LEAF || val_l2 == FULL_LEAF)
						return FULL_LEAF;
					else if(!val_l1 && !val_l2)
						return EMPTY_LEAF;
					else
						return VALUE_NEED_CHILDREN;
				}
			}
			case FUNCTOR_EXTEND: { // Extend quadtree to attributes att_set
				// i --> i' (mapping)
				return get_child_lqdag_aux(formula_pos->get_formula_lqdag1(), formula_pos->getM(ith_child),m_pos_qdag);
			}
			case FUNCTOR_LQDAG:{
				return get_child_lqdag_aux(formula_pos->get_formula_lqdag1(), ith_child, m_pos_qdag);
			}
			default:
				throw "error: value_lqdag non valid functor";
		}
	}


	// ------------------ VALUE FUNCTION ------------------ //

    /**
     * algorithm 7 value(L).
     * In lqdags we introduce a new idea: full leaves, that denote subgrids full of 1s.
     * We save the value in val_lqdag_1, val_lqdag_2 or subQuatree->value. If we've already computed the value, we return it.
     * @return the value of the root of the lqdag. Can be EMPTY_LEAF, FULL_LEAF, INTERNAL_NODE or VALUE_NEED_CHILDREN.
     */
    double value_lqdag(formula_lqdag* formula_pos){
		if(this->node_completion_lqdag->val_node != NO_VALUE_LEAF)
			return this->node_completion_lqdag->val_node;

        switch(formula_pos->get_functor()){
            case FUNCTOR_QTREE: {
				if(formula_pos->get_val_lqdag1() == NO_VALUE_LEAF)
					formula_pos->set_val_lqdag1(this->value_qdag(formula_pos,formula_pos->get_index()));
                return formula_pos->get_val_lqdag1();
            }
            case FUNCTOR_NOT: {
				if(formula_pos->get_val_lqdag1() == NO_VALUE_LEAF)
					formula_pos->set_val_lqdag1(1-this->value_qdag(formula_pos, formula_pos->get_index()));
                return formula_pos->get_val_lqdag1();
            }
            case FUNCTOR_AND: {
				if(formula_pos->get_val_lqdag1() == EMPTY_LEAF || formula_pos->get_val_lqdag2() == EMPTY_LEAF)
					return EMPTY_LEAF;
				if(formula_pos->get_val_lqdag1() == FULL_LEAF)
					return this->value_lqdag(formula_pos->get_formula_lqdag2());
				if(formula_pos->get_val_lqdag2() == FULL_LEAF)
					return this->value_lqdag(formula_pos->get_formula_lqdag1());
                return VALUE_NEED_CHILDREN; // if both roots have Value 1⁄2, one cannot be sure of the Value of the resulting root until the AND between the children of Q1 and Q2 has been computed
            }
            case FUNCTOR_OR: {
				if(formula_pos->get_val_lqdag1() == FULL_LEAF || formula_pos->get_val_lqdag2() == FULL_LEAF)
					return FULL_LEAF;
				if(formula_pos->get_val_lqdag1() == EMPTY_LEAF)
					return this->value_lqdag(formula_pos->get_formula_lqdag2());
				if(formula_pos->get_val_lqdag2() == EMPTY_LEAF)
					return this->value_lqdag(formula_pos->get_formula_lqdag1());
				return VALUE_NEED_CHILDREN; // if both roots have Value 1⁄2, one cannot be sure of the Value of the resulting root until the AND between the children of Q1 and Q2 has been computed
			}
            case FUNCTOR_EXTEND: {
                if (formula_pos->get_val_lqdag1() == NO_VALUE_LEAF)
					formula_pos->set_val_lqdag1(this->value_lqdag(formula_pos->get_formula_lqdag1()));
                return formula_pos->get_val_lqdag1();
            }
			case FUNCTOR_LQDAG:{
				if (formula_pos->get_val_lqdag1() == NO_VALUE_LEAF)
					formula_pos->set_val_lqdag1(this->value_lqdag(formula_pos->get_formula_lqdag1()));
				return formula_pos->get_val_lqdag1();
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
     * We need to go down till we arrive to a FUNCTOR_TREE or FUNCTOR_NOT, because the value of the lqdag will be different after we evaluate the entire formula.
     * @return EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE.
     */
    double value_qdag(formula_lqdag* formula_pos, uint64_t index_qdag){
        if(formula_pos->get_functor() == FUNCTOR_QTREE || formula_pos->get_functor() == FUNCTOR_NOT) { // only is called
			if(index_qdag >= pos_qdags.size())
				cout << "error: value_qdag index out of range" << endl;
			uint16_t cur_level = (this->pos_qdags).at(index_qdag).level;
            uint64_t cur_node = (this->pos_qdags).at(index_qdag).cur_node;
            if(formula_pos->get_val_lqdag1() == NO_VALUE_LEAF)
                return is_a_qdag_leaf(formula_pos, index_qdag, cur_level, cur_node);
            return formula_pos->get_val_lqdag1();
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
     * @return whether is an EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE. For a Q_TREE (if the functor is NOT, we will have to invert the value).
     */
    double is_a_qdag_leaf(formula_lqdag* formula_pos, uint64_t index_qdag, uint16_t level, uint64_t node){
        uint16_t max_level = formula_pos->getMaxLevel();
        if(level > max_level){
			if(index_qdag >= pos_qdags.size())
				cout << "error: 1) is_a_qdag_leaf index out of range" << endl;
            return arr_qdags.at(index_qdag).get_ith_bit(max_level, node);
        }
		if(index_qdag >= pos_qdags.size())
			cout << "error: 2) is_a_qdag_leaf index out of range" << endl;
        uint8_t node_description = arr_qdags.at(index_qdag).get_node_last_level(level, node);
        if((node_description | 0) == 0)
            return EMPTY_LEAF;
        uint16_t val_full_ones;
        uint64_t k_d = arr_qdags.at(index_qdag).getKD();
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
            uint64_t siblings = arr_qdags.at(index_qdag).rank(level,node);
            for(uint64_t i = 0; i < k_d; i++){
                if(is_a_qdag_leaf(formula_pos,index_qdag, level+1, (siblings+i)*k_d) != 1)
                    return INTERNAL_NODE;
            }
            return FULL_LEAF;
        }
        return INTERNAL_NODE;
    }

	/**
	 * Algorithm 9 of the paper.
	 * We compute the completion of the lqdag and store the values in the node_completion_lqdag.
	 * @param p
	 * @param max_level
	 * @param cur_level the current level of the quadtree. If cur_level = max_level + 1, we are looking for the bit representing the point in the grid.
	 * @param UPPER_BOUND
	 * @param results
	 * @return
	 * PD: we could delete max_level and cur_level, but it's more simple to visualize it this way
	 */
	lqdag* completion_dfs(uint16_t max_level,
						  uint16_t cur_level,
						  uint64_t UPPER_BOUND,
						  uint64_t &results){
		if(results >= UPPER_BOUND){
			return this;
		}
		double val_lqdag = this->value_lqdag(this->form_lqdag);
		this->node_completion_lqdag->val_node = val_lqdag;
		if(val_lqdag == EMPTY_LEAF || val_lqdag == FULL_LEAF){
			if(val_lqdag == FULL_LEAF){
				results += (cur_level > max_level) ? 1 :  pow(nChildren(), (max_level - cur_level + 1));
			}
			return this;
		}

		double max_value = 0;
		double min_value = 1;
//		lqdag** test = new lqdag*[nChildren()];
		for(uint64_t i = 0; i < nChildren(); i++){
			// compute the child and store it in node_completion->lqdag_children[i]
			lqdag* child_lqdag = this->get_child_lqdag(i);
			// recursive call
//			test[i] = child_lqdag;
			child_lqdag->completion_dfs(max_level, cur_level+1, UPPER_BOUND, results);
			max_value = max(max_value, child_lqdag->node_completion_lqdag->val_node);
			min_value = min(min_value, child_lqdag->node_completion_lqdag->val_node);
		}
//		cout << "cur level: "<< cur_level << endl;
//		for(uint64_t i = 0; i < nChildren(); i++){
//			cout << test[i]->node_completion_lqdag->val_node << " ";
//		}
//		cout << endl;
		// pruning step --> return a leaf
		// TODO: i think we dont need to do this, because we do it computing the full leaf
		if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){
			delete[] this->node_completion_lqdag->lqdag_children;
			this->node_completion_lqdag->lqdag_children = nullptr;
			double val_leaf = (max_value == EMPTY_LEAF) ? EMPTY_LEAF : FULL_LEAF;
			this->node_completion_lqdag->val_node = val_leaf;
		}
		return this;

	}

	/**
	 * Lazy evaluation of the predicate
	 * @param ith_child the i-th child we want to check if it satisfies the predicate.
	 * @param pred
	 * @return The i-th child if it satisfies the predicates. Otherwise, return an empty leaf.
	 */
	lqdag* get_child_selection(uint64_t ith_child, predicate *pred){
		if(this->node_completion_lqdag->val_node == EMPTY_LEAF || (this->node_completion_lqdag->lqdag_children[ith_child] && this->node_completion_lqdag->lqdag_children[ith_child]->node_completion_lqdag->val_node == EMPTY_LEAF)){
			return this;
		}

//		uint16_t diff_level = this->getMaxLevel() - this->getCurLevel();
//		uint16_t* coordinatesTemp = new uint16_t[this->nAttr()];
//		for(uint16_t k = 0; k < this->nAttr(); k++)
//			coordinatesTemp[k] = this->node_completion_lqdag->coordinates[k];
		// update the current coordinates to be the path to the ith_child
//		transformCoordinates(coordinatesTemp, nAttr(), diff_level, ith_child);
		uint256_t newPath = this->node_completion_lqdag->path;
		getNewMortonCodePath(newPath, this->nAttr(), this->getCurLevel(), (uint256_t) ith_child);

		uint64_t quadrant_side = this->getGridSide() /
								 pow(this->getK(),this->getCurLevel()+1);

		// we evaluate the predicate first for the ith-child (with the coordinates and the quadrant)
		double val_eval_pred = eval_pred(pred, newPath, quadrant_side, this->nAttr());
		lqdag* child_lqdag;
		if(val_eval_pred == 0){ // return an empty child
			child_lqdag = new lqdag(this->arr_qdags, this->form_lqdag);
			child_lqdag->node_completion_lqdag->val_node = EMPTY_LEAF;
			child_lqdag->node_completion_lqdag->level = this->node_completion_lqdag->level +1;
			// TODO: check this
			if(this->node_completion_lqdag->lqdag_children && this->node_completion_lqdag->lqdag_children[ith_child]){
				delete[] this->node_completion_lqdag->lqdag_children[ith_child]->node_completion_lqdag->lqdag_children;
				this->node_completion_lqdag->lqdag_children[ith_child]->node_completion_lqdag->lqdag_children = nullptr;
			}
		} else if(val_eval_pred == 1){ // if it passes the predicate, we compute the child
			child_lqdag = this->get_child_lqdag(ith_child);
			return child_lqdag;
		} else { // value = 0.5
			child_lqdag = this->get_child_lqdag(ith_child);
			// update the value of the child (because the value is 0.5, if the original ith-child were 1, after the condition will be 0.5
			if(child_lqdag->node_completion_lqdag->val_node == FULL_LEAF){
				child_lqdag->node_completion_lqdag->val_node = VALUE_NEED_CHILDREN;
				this->node_completion_lqdag->val_node = VALUE_NEED_CHILDREN;
			}

		}
		this->node_completion_lqdag->lqdag_children[ith_child] = child_lqdag;
		return child_lqdag;
	}

	lqdag* completion_selection_dfs(predicate *pred,
									uint16_t max_level,
									uint16_t cur_level,
									uint64_t UPPER_BOUND,
									uint64_t &results){
		if(results >= UPPER_BOUND || this->node_completion_lqdag->val_node == EMPTY_LEAF){
			return this;
		}

		double max_value = 0;
		double min_value = 1;
		for(uint64_t i = 0; i < nChildren(); i++){
			lqdag* child_lqdag_pred = this->get_child_selection(i, pred);
			if(child_lqdag_pred->node_completion_lqdag->val_node == FULL_LEAF){
				// TODO: check get cur level is updated
				results += (cur_level > max_level) ? 1 :  pow(nChildren(), (max_level - cur_level + 1));
//				results += (this->getCurLevel() > this->getMaxLevel()) ? 1 :  pow(nChildren(), (this->getMaxLevel() - this->getCurLevel() + 1));
			}
			child_lqdag_pred->completion_selection_dfs(pred, max_level, cur_level+1, UPPER_BOUND, results);
			max_value = max(max_value, child_lqdag_pred->node_completion_lqdag->val_node);
			min_value = min(min_value, child_lqdag_pred->node_completion_lqdag->val_node);
		}
		if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){
			delete[] this->node_completion_lqdag->lqdag_children;
			this->node_completion_lqdag->lqdag_children = nullptr;
			double val_leaf = (max_value == EMPTY_LEAF) ? EMPTY_LEAF : FULL_LEAF;
			this->node_completion_lqdag->val_node = val_leaf;
		}
		return this;

	}

};

#endif
// TODO: see diff with a not!! we have to push dow  the operand