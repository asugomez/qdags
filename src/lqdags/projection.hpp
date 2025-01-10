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
//	position** pos_qdags; // each qdag will have a space in this vector with its position (level and coordinates)
//	node_completion* node_completion_lqdag; // children of the lqdag (computed)
	att_set attribute_set_A; // A
	att_set attribute_set_A_prime; // A'
//	lqdag** lqdags; // lqdag where we apply the projection
	lqdag* lqdag_root; // root of the lqdag
	double* val_pi; // value of each OR
	projection** projection_children = nullptr;


public:

	projection() = default;

	projection(lqdag* lqdag_root, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
		this->lqdag_root = lqdag_root;

		// TODO: see if memory is correct allocated!
		this->form_or = new formula_pi(attribute_set_A, attribute_set_A_prime, lqdag_root);

		this->val_pi = new double[1 << attribute_set_A_prime.size()];

		for(uint64_t i = 0; i < (1 << attribute_set_A.size()); i++){
			this->val_pi[i] = NO_VALUE_LEAF;
		}

		projection_children = new projection*[1 << attribute_set_A_prime.size()];

		for(uint64_t i = 0; i < (1 << attribute_set_A_prime.size()); i++){
			projection_children[i] = nullptr;
		}

		projection* test_child = get_child_pi(0);
	}

	projection(lqdag* lqdag_root, formula_pi* form_or, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
		this->lqdag_root = lqdag_root;
		this->form_or = form_or;

		this->val_pi = new double[1 << attribute_set_A_prime.size()]; // check val_pi is 0
		for(uint64_t i = 0; i < (1 << attribute_set_A.size()); i++){
			this->val_pi[i] = NO_VALUE_LEAF;
		}
		projection_children = new projection*[1 << attribute_set_A_prime.size()];
		for(uint64_t i = 0; i < (1 << attribute_set_A_prime.size()); i++){
			projection_children[i] = nullptr;
		}
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

	uint8_t get_functor_child(uint64_t i){
		return this->get_OR_formula_child(i)->get_functor();
	}


	double value_child_pi(uint64_t i){
		if(this->val_pi[i] == NO_VALUE_LEAF)
			this->get_child_pi(i);
		return this->val_pi[i];
	}

	/**
	 * Compute the i-th child of the projection. Assuming the i-th child is not already computed (value 0 or 1)
	 * @param i
	 * @return
	 */
	projection* get_child_pi(uint64_t i){
		assert(i < (1 << this->attribute_set_A_prime.size()) );
		// crear el formular_or
		// create new projection(...)
		lqdag* lqdag_child_pi = get_child_pi_aux(i, this->get_OR_formula_child(i));
		if(lqdag_child_pi != nullptr)
			this->val_pi[i] = lqdag_child_pi->val_node();

		if(this->val_pi[i] == VALUE_NEED_CHILDREN && lqdag_child_pi != nullptr) // TODO: check this value is correct or we need to check 0.5 as well
			this->projection_children[i] = new projection(lqdag_child_pi, this->form_or, this->attribute_set_A, this->attribute_set_A_prime);

		return this->projection_children[i]; // TODO: check what if this is null
	}

	/**
	 * Compute the i-th child of the projection
	 * @param i
	 * @param formula_OR
	 * @return
	 */
	lqdag* get_child_pi_aux(uint64_t i,formula_lqdag* formula_OR){
		switch (formula_OR->get_functor()) {
			case FUNCTOR_LQDAG: {// leaf
				uint64_t index_child_pi = formula_OR->get_index();
				// TODO: test which lqdag!
//				lqdag* test = this->lqdag_root->get_child_lqdag(index_child_pi);
				lqdag* child_lqdag_pi = formula_OR->get_lqdag()->get_child_lqdag(this->form_or->get_index_children(index_child_pi));
				return child_lqdag_pi; //child_lqdag_pi->value_lqdag(formula_OR)
			}
			case FUNCTOR_OR: {
				if (formula_OR->get_val_lqdag1() == FULL_LEAF || formula_OR->get_val_lqdag2() == FULL_LEAF) {
					this->val_pi[i] = FULL_LEAF;
					return nullptr;
				}
				if (formula_OR->get_val_lqdag1() == EMPTY_LEAF) {
					return get_child_pi_aux(i, formula_OR->get_formula_lqdag2());
				}
				if (formula_OR->get_val_lqdag2() == EMPTY_LEAF) {
					return get_child_pi_aux(i, formula_OR->get_formula_lqdag1());
				} else {
					lqdag l1 = get_child_pi_aux(i, formula_OR->get_formula_lqdag1());
					lqdag l2 = get_child_pi_aux(i, formula_OR->get_formula_lqdag2());
					if (l1.val_node() == FULL_LEAF || l2.val_node() == FULL_LEAF)
						this->val_pi[i] = FULL_LEAF;
					else if (!l1.val_node() && !l2.val_node())
						this->val_pi[i] = EMPTY_LEAF;
					else
						this->val_pi[i] = VALUE_NEED_CHILDREN;
					return nullptr;
				}
			}
			default:
				throw "error: get_child_pi_aux non valid functor";
		}
	}
};

#endif //INCLUDED_PROJECTION
