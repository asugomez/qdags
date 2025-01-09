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
	lqdag** lqdags; // lqdag where we apply the projection
	double* val_pi; // value of each OR
	projection** projection_children;


public:

	projection() = default;

	projection(lqdag* lqdag_root, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
		for(uint64_t i = 0; i < (1 << attribute_set_A.size()); i++){
			lqdags[i] = new lqdag(lqdag_root); // TODO: check this is a copy!
		}
		this->form_or = new formula_pi(attribute_set_A, attribute_set_A_prime, lqdag_root);
//		this->pos_qdags = new position*[1 << attribute_set_A.size()];
		this->val_pi = new double[1 << attribute_set_A_prime.size()]; // check val_pi is 0

		for(uint64_t i = 0; i < (1 << attribute_set_A.size()); i++){
			this->val_pi[i] = NO_VALUE_LEAF;
		}

		for(uint64_t i = 0; i < (1 << attribute_set_A_prime.size()); i++){
			projection_children[i] = nullptr;
		}
	}

	projection(formula_pi* form_or, att_set attribute_set_A, att_set attribute_set_A_prime, lqdag* lqdag_PROBLEM){
		// TODO: not implemented yet --> see lqdag
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
		for(uint64_t i = 0; i < (1 << attribute_set_A.size()); i++){
//			lqdags[i] = new lqdag(lqdag_root);
		}
		this->form_or = form_or;

		this->val_pi = new double[1 << attribute_set_A_prime.size()]; // check val_pi is 0
		for(uint64_t i = 0; i < (1 << attribute_set_A.size()); i++){
			this->val_pi[i] = NO_VALUE_LEAF;
		}
		for(uint64_t i = 0; i < (1 << attribute_set_A_prime.size()); i++){
			projection_children[i] = nullptr;
		}
	}

	uint8_t get_functor(){
		return this->form_or->get_functor();
	}

	formula_lqdag* get_OR_formula_child(uint64_t i){
		return this->form_or->get_OR_formula_child(i);
	}

	uint8_t get_functor_child(uint64_t i){
		return this->get_OR_formula_child(i)->get_functor();
	}


	double value_child_pi(uint64_t i){
		return this->val_pi[i];
	}

	projection* get_child_pi(uint64_t i){
		assert(i < (1 << this->attribute_set_A_prime.size()) );
		// crear el formular_or
		// create new projection(...)
		this->val_pi[i] = get_child_pi_aux(i, this->get_OR_formula_child(i));
		this->projection_children[i] = new projection(this->form_or, this->attribute_set_A, this->attribute_set_A_prime, this->lqdags[i]);

	}

	double get_child_pi_aux(uint64_t i,formula_lqdag* formula_OR){ // TODO: see if we need i
		switch (formula_OR->get_functor()) {
			case FUNCTOR_LQDAG: {// leaf
				uint64_t index_child_pi = formula_OR->get_index();
				lqdag *child_lqdag_pi = formula_OR->get_lqdag()->get_child_lqdag(this->form_or->get_index_children(index_child_pi));
				return child_lqdag_pi->value_lqdag(formula_OR);
			}
			case FUNCTOR_OR: {
				if (formula_OR->get_val_lqdag1() == FULL_LEAF || formula_OR->get_val_lqdag2() == FULL_LEAF)
					return FULL_LEAF;
				if (formula_OR->get_val_lqdag1() == EMPTY_LEAF) {
					double val_l2 = get_child_pi_aux(i, formula_OR->get_formula_lqdag2());
					return val_l2;
				}
				if (formula_OR->get_val_lqdag2() == EMPTY_LEAF) {
					double val_l1 = get_child_pi_aux(i, formula_OR->get_formula_lqdag1());
					return val_l1;
				} else {
					double val_l1 = get_child_pi_aux(i, formula_OR->get_formula_lqdag1());
					double val_l2 = get_child_pi_aux(i, formula_OR->get_formula_lqdag2());
					if (val_l1 == FULL_LEAF || val_l2 == FULL_LEAF)
						return FULL_LEAF;
					else if (!val_l1 && !val_l2)
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
