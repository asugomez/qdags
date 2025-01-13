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

	struct pi_materialization{
		double val_node = NO_VALUE_LEAF;
		uint16_t level;
		uint64_t n_children; // number of children
		uint256_t path;
		projection** projection_children; // TODO: aun no estoy segura si usariamos esto

		pi_materialization(): val_node(NO_VALUE_LEAF), level(0), n_children(0), path(0), projection_children(nullptr) {}

		pi_materialization(uint64_t n_children): val_node(NO_VALUE_LEAF), level(0), n_children(n_children), path(0), projection_children(new projection * [n_children]) {
			for(uint64_t i = 0; i < n_children; i++){
				projection_children[i] = nullptr;
			}
		}
	};

private:
	formula_pi* form_or; // expression, syntax tree, formula of the OR tree.
	att_set attribute_set_A; // A
	att_set attribute_set_A_prime; // A'
	lqdag* lqdag_root; // the lqdag where we apply the projection
	double* val_pi; // value of each OR
	uint64_t number_projection_children; // 2^(A)
	projection** projection_children = nullptr;
	pi_materialization* materialization;


public:

	projection() = default;

	projection(lqdag* lqdag_root, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
		this->number_projection_children = 1 << attribute_set_A_prime.size();
		this->materialization = new pi_materialization(this->number_projection_children);
		this->lqdag_root = lqdag_root;

		// TODO: see if memory is correct allocated!
		this->form_or = new formula_pi(attribute_set_A, attribute_set_A_prime, lqdag_root);

		this->val_pi = new double[number_projection_children];

		for(uint64_t i = 0; i < number_projection_children; i++){
			this->val_pi[i] = NO_VALUE_LEAF;
		}

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
		this->number_projection_children = 1 << attribute_set_A_prime.size();
		this->materialization = new pi_materialization(this->number_projection_children);
		this->form_or = form_or;
		this->lqdag_root = lqdag_root;

		this->val_pi = new double[number_projection_children]; // check val_pi is 0
		for(uint64_t i = 0; i < number_projection_children; i++){
			this->val_pi[i] = NO_VALUE_LEAF;
		}
		projection_children = new projection*[number_projection_children];
		for(uint64_t i = 0; i < number_projection_children; i++){
			projection_children[i] = nullptr;
		}
	}

	uint8_t get_functor(){
		return this->form_or->get_functor();
	}

	/**
	 * The value of the root of the projection: 0, 1 or VALUE_NEED_CHILDREN
	 * @param val
	 */
	void set_val_node_pi(double val){
		this->materialization->val_node = val;
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

	/**
	 * Get the value of the current projection
	 * @return
	 */
	double get_value_pi(){
		double val_lqdag = this->lqdag_root->val_node();
		if(val_lqdag == FULL_LEAF || val_lqdag == EMPTY_LEAF)
			return val_lqdag;
		else
			return VALUE_NEED_CHILDREN;
	}

	/**
	 * Get the value of the i-th child of the projection
	 * @param i
	 * @return
	 */
	double get_value_child_pi(uint64_t i){
		return this->val_pi[i];
	}

	/**
	 * Compute the i-th child of the projection.
	 * Assuming the i-th child is not already computed (value 0 or 1).
	 * @param i
	 * @return
	 */
	projection* get_child_pi(uint64_t i){
		assert(i < (1 << this->attribute_set_A_prime.size()) );
		if(this->get_value_pi() == EMPTY_LEAF || this->get_value_pi() == FULL_LEAF)
			return this;

		double val_child_pi = get_child_pi_aux(i, this->get_OR_formula_child(i));

		this->val_pi[i] = val_child_pi;

		// create new projection

		if(val_child_pi != EMPTY_LEAF){
			// add coordinates to the i-th child

		}
		if(lqdag_child_pi != nullptr && this->val_pi[i] == VALUE_NEED_CHILDREN) // TODO: check this value is correct or we need to check 0.5 as well
			this->projection_children[i] = new projection(lqdag_child_pi, this->form_or, this->attribute_set_A, this->attribute_set_A_prime);

	}

	/**
	 * Compute the i-th child of the projection
	 * @param i
	 * @param formula_OR
	 * @return nullptr if the value is FULL_LEAF or EMPTY_LEAF, otherwise the lqdag of the child
	 */
	double get_child_pi_aux(uint64_t i,formula_lqdag* formula_OR){
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
