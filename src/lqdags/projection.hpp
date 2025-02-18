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
		this->lqdag_root = lqdag_root;

		// TODO: see if memory is correct allocated!
		this->form_or = new formula_pi(attribute_set_A, attribute_set_A_prime, lqdag_root);
		// TODO: maybe we only need the indexes children:
		// this->indexes_children = create_pi_index_children(attribute_set_A, attribute_set_A_prime);

		this->value_children = new double[nChildren()];
		for(uint64_t i = 0; i < nChildren(); i++){
			this->value_children[i] = NO_VALUE_LEAF;
		}

		uint64_t number_projection_children = 1 << attribute_set_A.size();

		projection_children = new projection*[number_projection_children];
		for(uint64_t i = 0; i < number_projection_children; i++){ // number of leaves
			projection_children[i] = nullptr;
		}

		// TODO: OFT
		projection* test_child = get_child_pi(0);
	}

	// TODO: maybe instead of the form, we need the indexes of the children
	projection(lqdag* lqdag_root, formula_pi* form_or, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
		this->form_or = form_or;
		this->lqdag_root = lqdag_root;

		this->value_children = new double[nChildren()];
		for(uint64_t i = 0; i < nChildren(); i++){
			this->value_children[i] = NO_VALUE_LEAF;
		}

		uint64_t number_projection_children = 1 << attribute_set_A.size();

		projection_children = new projection*[number_projection_children];
		for(uint64_t i = 0; i < number_projection_children; i++){ // number of leaves
			projection_children[i] = nullptr;
		}
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
		assert(i < this->nChildren());

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

	void compute_multi_projection_child_pi(uint64_t i){
		assert(i < this->nChildren());
		assert(this->value_children[i] == VALUE_NEED_CHILDREN);
		uint64_t n_leaves_per_children = 1 << (nAttrA() - nAttrA_prime());
		uint64_t init_index = i * n_leaves_per_children;
		for(uint64_t j = init_index; j < init_index + n_leaves_per_children; j++) {
			lqdag *lqdag_j = this->lqdag_root->get_child_lqdag(this->get_index_children(j));
			projection* child_pi = new projection(lqdag_j, this->form_or, this->attribute_set_A, this->attribute_set_A_prime);
			this->projection_children[j] = child_pi;
		}

	}

	projection* get_ith_projection(uint64_t i){
		assert(i < 1 << nAttrA());
		return this->projection_children[i];
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
				compute_multi_projection_child_pi(i);
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


};

#endif //INCLUDED_PROJECTION
