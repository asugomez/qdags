//
// Created by Asunción Gómez on 08-01-25.
//

#ifndef FORMULA_PI_HPP
#define FORMULA_PI_HPP

#include<bits/stdc++.h>
#include "formula_lqdag.hpp"
#include "utils_lqdags.hpp"

//const uint8_t FUNCTOR_OR = 3;
const uint8_t FUNCTOR_PI = 5; // root of the projection formula

class formula_pi{ // tree of ORs for the projection operation

private:
	uint8_t functor; // PI
	formula_lqdag** formula_pi_children; // formula of ORs
	uint64_t number_children_pi; // number of OR in the first level
	uint16_t max_level_formula; // height - 1 of the formula (tree of ORs)
	uint64_t number_leaves; // number of lqdag_children at the last level (2^d)
	uint16_t nAttr_A;
	uint16_t nAttr_A_prime;
	indexes indexes_children; // indexes of the j for each i

public:

	formula_pi() = default;

	formula_pi(att_set attribute_set_A, att_set attribute_set_A_prime, lqdag* lqdag){
		this->functor = FUNCTOR_PI;
		this->nAttr_A = attribute_set_A.size();
		this->nAttr_A_prime = attribute_set_A_prime.size();

		this->number_children_pi = 1 << nAttr_A_prime;
		this->max_level_formula = nAttr_A - nAttr_A_prime;
		this->number_leaves = 1 << nAttr_A;
		this->indexes_children = create_pi_index_children(attribute_set_A, attribute_set_A_prime);

		uint64_t i;
		uint16_t current_level = max_level_formula;
		uint64_t current_number_OR = number_leaves >> 1;

		formula_lqdag* lqdags_leaves_formula[number_leaves];
		// TODO: check level_OR_Formula!
		formula_lqdag* level_OR_formula[max_level_formula][current_number_OR];
		for(i = 0; i < number_leaves; i++){ // create the last level of the tree
			lqdags_leaves_formula[i] = new formula_lqdag(FUNCTOR_LQDAG, lqdag, i); // TODO: check index lqdag
			if(i%2 == 1){
				level_OR_formula[current_level-1][i/2] = new formula_lqdag(FUNCTOR_OR, lqdags_leaves_formula[i-1], lqdags_leaves_formula[i]);
			}
		}
		current_level --;
		current_number_OR >>= 1;

		while(current_level > 0){
			for(i = 0; i < current_number_OR + 1; i++) {
				level_OR_formula[current_level-1][i] = new formula_lqdag(FUNCTOR_OR, level_OR_formula[current_level][2*i], level_OR_formula[current_level][2*i + 1]);
			}
			current_level --;
			current_number_OR >>= 1;
		}

		this->formula_pi_children = level_OR_formula[0];



	}

	uint8_t get_functor() const {
		return this->functor;
	}

	formula_lqdag* get_OR_formula_child(uint64_t i){
		return formula_pi_children[i];
	}

	uint64_t get_index_children(uint64_t i){
		return indexes_children[i];
	}

};

#endif //FORMULA_PI_HPP
