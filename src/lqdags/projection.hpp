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
	position** pos_qdags; // each qdag will have a space in this vector with its position (level and coordinates)
//	node_completion* node_completion_lqdag; // children of the lqdag (computed)
	att_set attribute_set_A; // A
	att_set attribute_set_A_prime; // A'
	lqdag* lqdag_root; // lqdag where we apply the projection

public:

	projection() = default;

	projection(lqdag* lqdag_root, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
		this->lqdag_root = lqdag_root;
		lqdag** lqdag_children = new lqdag*[1 << attribute_set_A.size()];
		for(uint64_t i = 0; i < 1 << attribute_set_A.size(); i++){
			lqdag_children[i] = new lqdag(lqdag_root);
		}
		this->form_or = new formula_pi(attribute_set_A.size(), attribute_set_A_prime.size(), lqdag_children);
		this->pos_qdags = new position*[1 << attribute_set_A.size()];
	}

};

#endif //INCLUDED_PROJECTION
