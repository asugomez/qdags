//
// Created by Asunción Gómez on 06-01-25.
//
#ifndef INCLUDED_PROJECTION
#define INCLUDED_PROJECTION

#include "../lqdag.hpp"
#include "formula_pi.hpp"
#include "utils_lqdags.hpp"

class projection {

public:

	typedef vector<uint64_t> att_set;
	

private:
	formula_pi* form_or; // expression, syntax tree, formula of the OR tree.
	att_set attribute_set_A; // A
	att_set attribute_set_A_prime; // A'
	lqdag* lqdag_root; // the lqdag where we apply the projection
	double val_pi = NO_VALUE_LEAF; // value of each OR
	projection** projection_children; // it has 2^|A| children (leaves of projection)
	double* value_children;

public:

	projection() = default;

	projection(lqdag* lqdag_root, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
//		this->materialization = new pi_materialization(this->number_projection_children);
		this->lqdag_root = lqdag_root;

		// TODO: see if memory is correct allocated!
		this->form_or = new formula_pi(attribute_set_A, attribute_set_A_prime, lqdag_root);

		this->value_children = new double[nChildren()];

		for(uint64_t i = 0; i < nChildren(); i++){
			this->value_children[i] = NO_VALUE_LEAF;
		}

		uint64_t number_projection_children = 2 << attribute_set_A.size();

		projection_children = new projection*[number_projection_children];

		for(uint64_t i = 0; i < number_projection_children; i++){ // number of leaves
			projection_children[i] = nullptr;
		}

		// TODO: OFT
		projection* test_child = get_child_pi(0);
	}

	projection(lqdag* lqdag_root, formula_pi* form_or, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
//		this->number_projection_children = 1 << attribute_set_A_prime.size();
//		this->materialization = new pi_materialization(this->number_projection_children);
		this->form_or = form_or;
		this->lqdag_root = lqdag_root;

//		this->val_pi = new double[number_projection_children]; // check val_pi is 0
//		for(uint64_t i = 0; i < number_projection_children; i++){
//			this->val_pi[i] = NO_VALUE_LEAF;
//		}
//		projection_children = new projection*[number_projection_children];
//		for(uint64_t i = 0; i < number_projection_children; i++){
//			projection_children[i] = nullptr;
//		}
	}

	uint64_t nAttrA(){
		return this->attribute_set_A.size();
	}

	uint64_t nAttrA_prime(){
		return this->attribute_set_A_prime.size();
	}

	uint64_t nChildren(){
		return (1 << nAttrA_prime()); // 2^|A'|
	}

	uint64_t number_projection_leaves(){
		return this->form_or->get_number_leaves();
	}

	uint8_t get_functor(){
		return this->form_or->get_functor();
	}

//	uint16_t getCurLevel(){
//		return this->materialization->level;
//	}

	/**
	 * The value of the root of the projection: 0, 1 or VALUE_NEED_CHILDREN
	 * @param val
	 */
//	void set_val_node_pi(double val){
//		this->materialization->val_node = val;
//	}

	/**
	 * Get the root of the i-th OR of the projection
	 * @param i
	 * @return
	 */
	formula_lqdag* get_OR_formula_child(uint64_t i){
		return this->form_or->get_OR_formula_child(i);
	}

	/**
	 * Get the correspondant index of the i-th child of the lqdag.
	 * @param i
	 * @return
	 */
	uint64_t get_index_children(uint64_t i){
		return this->form_or->get_index_children(i);
	}


	/**
	 * Get the value of the current projection (EMPTY_LEAF, FULL_LEAF or VALUE_NEED_CHILDREN)
	 * if the value of the original lqdag is 0 or 1, return that value. Otherwise, return VALUE_NEED_CHILDREN
	 * @return
	 */
	double value_pi(){
		double val_lqdag = this->lqdag_root->val_node();
		if(val_lqdag == FULL_LEAF || val_lqdag == EMPTY_LEAF){
			this->val_pi = val_lqdag;
		}
		return this->val_pi; // initially it will return NO_VALUE_LEAF, and then VALUE_NEED_CHILDREN
	}

	/**
	 * Compute a multi-OR between the lqdags chidlren of the i-th projection
	 * @param i the i-th child of the projection.
	 * @return EMPTY_LEAF, FULL_LEAF or VALUE_NEED_CHILDREN
	 */
	double get_value_multi_or_child_pi(uint64_t i){
		assert(i < this->nChildren() );

		if(this->value_children[i] != NO_VALUE_LEAF){
			return this->value_children[i];
		}

		double max_value = 0;
		double min_value = 1;
		uint64_t n_leaves_per_children = 1 << (nAttrA() - nAttrA_prime());
		uint64_t init_index = i * n_leaves_per_children;
		for(uint64_t j = init_index; j < init_index + n_leaves_per_children; j++) {
			lqdag *lqdag_j = this->lqdag_root->get_child_lqdag(this->get_index_children(j));
			double val_lqdag_j = lqdag_j->value_lqdag(lqdag_j->get_formula());
			max_value = max(max_value, val_lqdag_j);
			min_value = min(min_value, val_lqdag_j);
		}
		double val_or;
		if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){
			val_or = (max_value == EMPTY_LEAF) ? EMPTY_LEAF : FULL_LEAF;
		}
		else{
			val_or = VALUE_NEED_CHILDREN;
		}
		return val_or;

	}

	/**
	 *
	 * @return a projection with the output as a traditional quadtree (with the values of val_pi)
	 */
	projection* eval_projection(){
		double val_pi = this->value_pi(); // normally it has to be VALUE_NEED_CHILDREN
		double max_value = 0;
		double min_value = 1;
		if(val_pi == EMPTY_LEAF || val_pi == FULL_LEAF){
			return this;
		}

		for(uint64_t i = 0; i < nChildren(); i++){
			double val_child_pi = get_value_multi_or_child_pi(i); // evaluate an OR
			value_children[i] = val_child_pi;
			if(val_child_pi != EMPTY_LEAF && val_child_pi != FULL_LEAF){
				projection* child_pi = get_child_pi(i);
				child_pi->eval_projection();
			}
			else{
				// is a EMPTY_LEAF or a FULL_LEAF, no further computation is needed
				// TOOD: count the results?
			}
			max_value = max(max_value, val_child_pi);
			min_value = min(min_value, val_child_pi);
		}

		if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){
			delete [] projection_children;
			delete [] value_children;
			double val_node = (max_value == EMPTY_LEAF) ? EMPTY_LEAF : FULL_LEAF;
			this->val_pi = val_node;
		}
	}

	projection* get_ith_projection(uint64_t i){
		assert(i < 1 << this->attribute_set_A.size());
		lqdag* ith_child_lqdag = this->lqdag_root->get_child_lqdag(this->get_index_children(i));
		projection* child_pi = new projection(ith_child_lqdag, this->form_or, this->attribute_set_A, this->attribute_set_A_prime);
		return child_pi;

	}

	/**
	 * Compute the i-th child of the projection.
	 * Assuming the i-th child is not already computed (value 0 or 1).
	 * @param i
	 * @return
	 */
	projection* get_child_pi(uint64_t i){
		assert(i < nChildren() ); // TODO check assert

		if(this->value_pi() == EMPTY_LEAF || this->value_pi() == FULL_LEAF)
			return this;

		// check if the i-th child is already computed in projection_children, then return that
		if(this->value_pi() == VALUE_NEED_CHILDREN
			&& this->projection_children[0] != nullptr
			&& this->projection_children[i] != nullptr){
			return  this->projection_children[i];
		}

		// TODO: it has to create X projection that will be a multi - OR projection

		double val_child_pi = get_child_pi_aux(i, this->get_OR_formula_child(i));
		projection* child_pi = get_ith_projection(i);
		child_pi->val_pi = val_child_pi;

		if(val_child_pi != EMPTY_LEAF){
			// add coordinates to the i-th child
			uint256_t newPath = this->materialization->path;
			getNewMortonCodePath(newPath, nChildren, this->getCurLevel(), (uint256_t)i);
			// TODO: add this coordinates to the new projection
		}

		// TODO: what to return!
//		if(lqdag_child_pi != nullptr && this->val_pi[i] == VALUE_NEED_CHILDREN) // TODO: check this value is correct or we need to check 0.5 as well
//			this->projection_children[i] = new projection(lqdag_child_pi, this->form_or, this->attribute_set_A, this->attribute_set_A_prime);
	}

	/**
	 * Compute the i-th child of the projection
	 * @param i
	 * @param formula_OR
	 * @return nullptr if the value is FULL_LEAF or EMPTY_LEAF, otherwise the lqdag of the child
	 */
	double get_child_pi_aux(uint64_t i, formula_lqdag* formula_OR){
		switch (formula_OR->get_functor()) {
			case FUNCTOR_LQDAG: {// leaf
				uint64_t index_child_pi = formula_OR->get_index();
				// TODO: test which lqdag!
//				lqdag* test = this->lqdag_root->get_child_lqdag(index_child_pi);
				lqdag* child_lqdag_pi = formula_OR->get_lqdag()->get_child_lqdag(this->form_or->get_index_children(index_child_pi));
				double val_child_lqdag_pi = child_lqdag_pi->val_node(); // TODO: check this wont be 0.5
				if(val_child_lqdag_pi != FULL_LEAF && val_child_lqdag_pi != EMPTY_LEAF){// create a projection
					this->projection_children[index_child_pi] = new projection(child_lqdag_pi, this->form_or, this->attribute_set_A, this->attribute_set_A_prime);
					this->projection_children[index_child_pi]->set_val_node_pi(val_child_lqdag_pi);
				}
				return VALUE_NEED_CHILDREN;
			}
			case FUNCTOR_OR: {
				if (formula_OR->get_val_lqdag1() == FULL_LEAF || formula_OR->get_val_lqdag2() == FULL_LEAF) {
					return FULL_LEAF;
				}
				if (formula_OR->get_val_lqdag1() == EMPTY_LEAF) {
					// TODO: save the value only if is 0,1 otherwise put VALUE_NEED_CHILDREN
//					formula_OR->set_val_lqdag1(XX);
					double val_pi_2 = get_child_pi_aux(i, formula_OR->get_formula_lqdag2());
					return val_pi_2;
				}
				if (formula_OR->get_val_lqdag2() == EMPTY_LEAF) {
					double val_pi_1 = get_child_pi_aux(i, formula_OR->get_formula_lqdag1());
					return val_pi_1;
				} else {
					double val_pi_1 = get_child_pi_aux(i, formula_OR->get_formula_lqdag1());
					double val_pi_2 = get_child_pi_aux(i, formula_OR->get_formula_lqdag2());
					if(val_pi_1 == FULL_LEAF || val_pi_2 == FULL_LEAF)
						return FULL_LEAF;
					else if(!val_pi_1 && !val_pi_2)
						return EMPTY_LEAF;
					else
						return VALUE_NEED_CHILDREN;
				}
			}
			default:
				throw "error: get_child_pi_aux non valid functor";
		}
	}
};

#endif //INCLUDED_PROJECTION
