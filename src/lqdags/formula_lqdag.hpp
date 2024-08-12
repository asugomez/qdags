//
// Created by Asunción Gómez on 24-07-24.
//

#ifndef QDAGS_FORMULA_HPP
#define QDAGS_FORMULA_HPP

#include<bits/stdc++.h>
//#include "../qdags.hpp"
#include "../lqdag.hpp"

const uint8_t FUNCTOR_QTREE = 0; // leaf
const uint8_t FUNCTOR_LQDAG = 5; // leaf
const uint8_t FUNCTOR_NOT = 1; // leaf
const uint8_t FUNCTOR_AND = 2; // internal quadtree_formula
const uint8_t FUNCTOR_OR = 3; // internal quadtree_formula
const uint8_t FUNCTOR_EXTEND = 4; // internal quadtree_formula


class formula_lqdag{

public:
	typedef std::vector<uint64_t> att_set;
	typedef uint8_t type_mapping_M;

private:
	// lqdag L=(f,o), where f is a functor.
	uint64_t index; // position for the QTREE, or LQDAG
	uint8_t functor; // QTREE, NOT, AND, OR, EXTEND, LQDAG
	formula_lqdag *formula_lqdag1; // we save value of the left lqdag and also of the leaves (QTREE, NOT and LQDAG)
	formula_lqdag *formula_lqdag2;
	qdag* formula_leaf_qdag;
//	lqdags::lqdag* formula_leaf_lqdag;
	lqdag* formula_leaf_lqdag;
//	std::variant<qdag*, lqdag*> formula_leaf;
	double val_lqdag1 = NO_VALUE_LEAF;
	double val_lqdag2 = NO_VALUE_LEAF;
	uint8_t k;
	uint16_t nAtt;
	uint16_t max_level;
	uint64_t grid_side;
	type_mapping_M *M; // mapping for extend functor. Use it when we need to extend the subquadtree, Otherwise M[i] = i (case not extended).

public:

	formula_lqdag() = default;

	// F = (QTREE, Q_r)
	formula_lqdag(qdag *qdag) {
		this->functor = FUNCTOR_QTREE;
		this->formula_lqdag1 = nullptr;
		this->formula_lqdag2 = nullptr;
		this->formula_leaf_qdag = qdag;
		this->k = qdag->getK();
		this->nAtt = qdag->nAttr();
		this->max_level = qdag->getHeight()-1;
		this->grid_side = qdag->getGridSide();
	}

	// F = (QTREE, Q_r) o F = (NOT, Q_r)
	formula_lqdag(uint8_t functor, qdag *qdag) {
		assert(functor == FUNCTOR_QTREE || functor == FUNCTOR_NOT);
		this->functor = functor;
		this->formula_lqdag1 = nullptr;
		this->formula_lqdag2 = nullptr;
		this->formula_leaf_qdag = qdag;
		this->k = qdag->getK();
		this->nAtt = qdag->nAttr();
		this->max_level = qdag->getHeight()-1;
		this->grid_side = qdag->getGridSide();
	}

	// F = (LQDAG, L_q)
	formula_lqdag(lqdag *lqdag) {
		this->functor = FUNCTOR_LQDAG;
		this->formula_lqdag1 = nullptr;
		this->formula_lqdag2 = nullptr;
		this->formula_leaf_lqdag = lqdag;

		this->k = formula_leaf_lqdag->get_formula().getK();
		this->nAtt = formula_leaf_lqdag->get_formula().nAttr();
		this->max_level = formula_leaf_lqdag->get_formula().getHeight() - 1;
		this->grid_side = formula_leaf_lqdag->get_formula().getGridSide();
	}

	// F = (LQDAG, L_q)
	formula_lqdag(uint8_t functor, lqdag *lqdag) {
		assert(functor == FUNCTOR_LQDAG);
		this->functor = functor;
		this->formula_lqdag1 = nullptr;
		this->formula_lqdag2 = nullptr;
		this->formula_leaf_lqdag = lqdag;
		this->k = formula_leaf_lqdag->get_formula().getK();
		this->nAtt = formula_leaf_lqdag->get_formula().nAttr();
		this->max_level = formula_leaf_lqdag->get_formula().getHeight() - 1;
		this->grid_side = formula_leaf_lqdag->get_formula().getGridSide();
	}

	// F = (AND, F1, F2) o F = (OR, F1, F2)
	formula_lqdag(uint8_t functor, formula_lqdag *l1, formula_lqdag *l2, uint8_t k = 0, uint16_t nAtt = 0, uint16_t max_level = 0, uint64_t grid_side = 0) {
		assert(functor == FUNCTOR_AND || functor == FUNCTOR_OR);
		this->functor = functor;
		this->formula_lqdag1 = l1;
		this->formula_lqdag2 = l2;
		this->k = k;
		this->nAtt = nAtt;
		this->max_level = max_level;
		this->grid_side = grid_side;
	}


	// L = (EXTEND, L1, A)
	formula_lqdag(uint8_t functor, formula_lqdag* l, uint8_t k, uint16_t nAtt, uint16_t max_level=0, uint64_t grid_side=0){
		assert(functor == FUNCTOR_EXTEND);
		this->functor = functor;
		this->formula_lqdag1 = l;
		this->formula_lqdag2 = nullptr;
		this->k = k;
		this->nAtt = nAtt; // d
		this->max_level = max_level;
		this->grid_side = grid_side;

		uint64_t p = std::pow(k, nAtt);

		if(this->formula_lqdag1->is_qdag()) {
			if(M == nullptr) {
//				uint64_t msize = std::pow(this->formula_lqdag1->get_qdag()->getK(),
//										  this->formula_lqdag1->get_qdag()->getD());

//				uint16_t dim = this->attribute_set_A.size(); // d
				uint16_t dim_prime = this->formula_lqdag1->get_qdag()->nAttr(); // d'

				M = new type_mapping_M[p];

				uint64_t mask;
				uint64_t i, i_prime;

				for (i = 0; i < p; ++i) {
					// todos los bits están en cero excepto el bit en la posición dim_prime - 1.
					mask = 1 << (dim_prime - 1); // equivalent to 2^(dim_prime-1)
					i_prime = 0;

					for (uint16_t j = 0; j < dim_prime; ++j) {
						if (i & (1 << (nAtt - this->formula_lqdag1->get_qdag()->getAttr(j) - 1)))
							i_prime |= mask;

						mask >>= 1;
					}

					M[i] = i_prime; // = M[i_prime];
				}
			}
			this->M = M;

			if (nAtt == 3)
				this->formula_lqdag1->get_qdag()->createTableExtend3();
			else if (nAtt == 4)
				this->formula_lqdag1->get_qdag()->createTableExtend4();
			else if (nAtt == 5)
				this->formula_lqdag1->get_qdag()->createTableExtend5();
			else {
				cout << "EXTEND only works for queries of up to 5 attributes..." << endl;
				exit(1);
			}
		}
		else{
			// TODO: ver el extend para LQDAG como hoja
		}
	}

	uint64_t get_index(){
		return this->index;
	}

	uint8_t get_functor() const {
		return this->functor;
	}

	formula_lqdag* get_formula_lqdag1() const {
		return this->formula_lqdag1;
	}

	formula_lqdag* get_formula_lqdag2() const {
		return this->formula_lqdag2;
	}

	qdag* get_qdag() const {
		if(is_qdag())
			return formula_leaf_qdag;
	}

	lqdag* get_lqdag() const {
		if(is_lqdag())
			return this->formula_leaf_lqdag;
	}

	uint8_t getK() const {
		return this->k;
	}

	uint16_t nAttr() const {
		return this->nAtt;
	}

	uint16_t getMaxLevel() const {
		return this->max_level;
	}

	uint64_t getGridSide() const {
		return this->grid_side;
	}

	double get_val_lqdag1() const {
		return this->val_lqdag1;
	}


	double get_val_lqdag2() const {
		return this->val_lqdag2;
	}

	void set_val_lqdag1(double newVal) {
		this->val_lqdag1 = newVal;
	}

	void set_val_lqdag2(double newVal){
		this->val_lqdag2 = newVal;
	}

	bool is_qdag() const { // Or functor is FUNCTOR_QTREE
		return functor == FUNCTOR_QTREE;
	}

	bool is_lqdag() const {
		return functor == FUNCTOR_LQDAG;
//		return std::holds_alternative<lqdag*>(this->formula_leaf);
	}


	// TODO change it to void
	/**
	 *  Return the indexes of the corresponding child for each branch of the syntax tree.
	 * @param i
	 * @param new_i vector with the new indexes from left to right
	 * @return
	 */
	uint16_t get_index_i(uint16_t i, vector<uint16_t>& new_i){
		switch (this->get_functor()) {
			case FUNCTOR_QTREE: case FUNCTOR_NOT: {
				new_i.push_back(i);
			}
			case FUNCTOR_LQDAG:{
				new_i.push_back(this->formula_leaf_lqdag->get_formula()->get_index_i(i, new_i));
			}
			case FUNCTOR_AND: case FUNCTOR_OR: {
				new_i.push_back(this->formula_lqdag1->get_index_i(i, new_i ));
				new_i.push_back(this->formula_lqdag2->get_index_i(i, new_i));
			}
			case FUNCTOR_EXTEND: { // Extend quadtree to attributes att_set
				// i --> i' (the mapping)
				new_i.push_back(this->formula_lqdag1->get_index_i(this->M[i], new_i));
			}
			default:
				throw "error: value_lqdag non valid functor";

		}
	}

	// algorithm 8 child(L,i)
	/**
	 * Can only be invoked when value is 1/2 or VALUE_NEED_CHILDREN.
	 * The formula will remain the same (we use a pointer to the original formula)
	 * @param i
	 * @return
	 */
	formula_lqdag* get_formula_child_lqdag(uint64_t i){
		// case base
		switch (this->get_functor()) {
			case FUNCTOR_QTREE: case FUNCTOR_NOT: {
				uint16_t cur_level = this->subQuadtree->level;
				uint64_t cur_node = this->subQuadtree->node;
				uint16_t max_level = this->subQuadtree->qdag->getHeight() - 1;
				bool node_exists = true;
				uint64_t ith_child;
				double value;
				if(cur_level == max_level) {
					ith_child = cur_node + i;
					node_exists = this->subQuadtree->qdag->get_ith_bit(cur_level, cur_node + i);
					value = node_exists ? FULL_LEAF : EMPTY_LEAF;
				} else{
					// ith_child: start position of the ith child of cur_node in the next level (cur_level + 1)
					ith_child = this->subQuadtree->qdag->get_child(cur_level, cur_node, i, node_exists);
					value = node_exists ? NO_VALUE_LEAF : EMPTY_LEAF; // if node exists, we will have to compute its value after
				}
				cur_level++;
				subQuadtreeChild *subQuadtree = new subQuadtreeChild{this->subQuadtree->qdag, cur_level, ith_child, value};
				lqdag *l = new lqdag(this->form_lqdag, subQuadtree);
				return l;
			}
			case FUNCTOR_AND: {
				double val_l1 = this->val_lqdag1 != NO_VALUE_LEAF ? this->val_lqdag1 : this->lqdag1->value_lqdag();
				if (val_l1== 1)
					return lqdag2->get_child_lqdag(i);
				double val_l2 = this->val_lqdag2 != NO_VALUE_LEAF ? this->val_lqdag2 : this->lqdag2->value_lqdag();
				if (val_l2 == 1)
					return lqdag1->get_child_lqdag(i);
				lqdag *l = new lqdag(this->form_lqdag, this->lqdag1->get_child_lqdag(i), this->lqdag2->get_child_lqdag(i));
				return l;
			}
			case FUNCTOR_OR: {
				double val_l1 = this->val_lqdag1 != NO_VALUE_LEAF ? this->val_lqdag1 : this->lqdag1->value_lqdag();
				if (val_l1 == 0)
					return lqdag2->get_child_lqdag(i);
				double val_l2 = this->val_lqdag2 != NO_VALUE_LEAF ? this->val_lqdag2 : this->lqdag2->value_lqdag();
				if (val_l2 == 0)
					return lqdag1->get_child_lqdag(i);
				lqdag *l = new lqdag(this->form_lqdag, this->lqdag1->get_child_lqdag(i), this->lqdag2->get_child_lqdag(i));
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



};

#endif //QDAGS_FORMULA_HPP
