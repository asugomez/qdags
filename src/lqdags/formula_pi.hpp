//
// Created by Asunción Gómez on 08-01-25.
//

#ifndef FORMULA_PI_HPP
#define FORMULA_PI_HPP

#include<bits/stdc++.h>
#include "formula_lqdag.hpp"

//const uint8_t FUNCTOR_OR = 3;
const uint8_t FUNCTOR_PI = 5; // root of the projection formula

class formula_pi{ // tree of ORs for the projection operation

private:
	uint8_t functor; // PI
	formula_lqdag** formula_pi_children; // formula of ORs
	uint64_t number_OR_first_level; // number of OR in the first level
	uint16_t max_level_formula; // height - 1 of the formula (tree of ORs)
	lqdag** lqdag_children;
	uint64_t number_leaves; // number of lqdag_children at the last level (2^d)
	uint16_t nAttr_A;
	uint16_t nAttr_A_prime;

public:

	formula_pi() = default;

	formula_pi(uint16_t nAttr_A, uint16_t nAttr_A_prime, lqdag** lqdag_children){
		this->functor = FUNCTOR_PI;
		this->number_OR_first_level = 1 << nAttr_A_prime;
		this->max_level_formula = nAttr_A - nAttr_A_prime;
		this->lqdag_children = lqdag_children;
		this->number_leaves = 1 << nAttr_A;
		this->nAttr_A = nAttr_A;
		this->nAttr_A_prime = nAttr_A_prime;

		uint64_t i;
		uint16_t current_level = max_level_formula - 1;
		uint64_t current_number_OR = number_leaves >> 1; // TODO: check

		formula_lqdag* lqdags_leaves_formula[1 << nAttr_A];
		formula_lqdag* level_OR_formula[max_level_formula][current_number_OR];
		for(i = 0; i < number_leaves + 1; i++){
			lqdags_leaves_formula[i] = new formula_lqdag(FUNCTOR_LQDAG, lqdag_children[i]);
			if(i%2 == 1){ // check i%2
				level_OR_formula[current_level][i/2] = new formula_lqdag(FUNCTOR_OR, lqdags_leaves_formula[i-1], lqdags_leaves_formula[i]);
			}
		}
		current_level --;
		current_number_OR >>= 1; // TODO: check

		while(current_level >= 0){
			for(i = 0; i < current_number_OR + 1; i++) {
				level_OR_formula[current_level][i] = new formula_lqdag(FUNCTOR_OR, level_OR_formula[current_level][2*i], level_OR_formula[current_level][2*i + 1]);
			}
			current_level --;
		}

		this->formula_pi_children = level_OR_formula[0];


	}


};

#endif //FORMULA_PI_HPP
