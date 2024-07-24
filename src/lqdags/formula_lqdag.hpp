//
// Created by Asunción Gómez on 24-07-24.
//

#ifndef QDAGS_FORMULA_HPP
#define QDAGS_FORMULA_HPP

#include<bits/stdc++.h>
#include "../qdags.hpp"
#include "../lqdag.hpp"

const uint8_t FUNCTOR_QTREE = 0; // leaf
const uint8_t FUNCTOR_NOT = 1; // leaf
const uint8_t FUNCTOR_AND = 2; // internal quadtree_formula
const uint8_t FUNCTOR_OR = 3; // internal quadtree_formula
const uint8_t FUNCTOR_EXTEND = 4; // internal quadtree_formula

const double NO_VALUE_LEAF = 3;
const double EMPTY_LEAF = 0;
const double FULL_LEAF = 1;
const double INTERNAL_NODE = 0.5;

class formula_lqdag{
public:
	typedef std::vector<uint64_t> att_set;
	typedef uint8_t type_mapping_M;

private:
	// lqdag L=(f,o), where f is a functor.
	uint8_t functor; // QTREE, NOT, AND, OR, EXTEND
	formula_lqdag *formula_lqdag1;
	formula_lqdag *formula_lqdag2;
	std::variant<qdag*, lqdag*> formula_leaf;
	double val_lqdag1 = NO_VALUE_LEAF;
	double val_lqdag2 = NO_VALUE_LEAF;
	att_set attribute_set_A; // for extend functor
	type_mapping_M *M; // mapping for extend functor. Use it when we need to extend the subquadtree, Otherwise M[i] = i (case not extended).

public:

	formula_lqdag() = default;

	// F = (QTREE, Q_r)
	formula_lqdag(qdag *qdag) {
		this->functor = FUNCTOR_QTREE;
		this->formula_lqdag1 = nullptr;
		this->formula_lqdag2 = nullptr;
		this->formula_leaf = qdag;
	}

	// F = (LQDAG, L_q)
	formula_lqdag(lqdag *lqdag) {
		this->functor = FUNCTOR_QTREE;
		this->formula_lqdag1 = nullptr;
		this->formula_lqdag2 = nullptr;
		this->formula_leaf = lqdag;
	}

	// F = (QTREE, Q_r) o F = (NOT, Q_r)
	formula_lqdag(uint8_t functor, qdag *qdag) {
		assert(functor == FUNCTOR_QTREE || functor == FUNCTOR_NOT);
		this->functor = functor;
		this->formula_lqdag1 = nullptr;
		this->formula_lqdag2 = nullptr;
		this->formula_leaf = qdag;
	}

	// F = (LQDAG, L_q) o F = (NOT, L_q)
	formula_lqdag(uint8_t functor, lqdag *lqdag) {
		assert(functor == FUNCTOR_QTREE || functor == FUNCTOR_NOT);
		this->functor = functor;
		this->formula_lqdag1 = nullptr;
		this->formula_lqdag2 = nullptr;
		this->formula_leaf = lqdag;
	}

	// F = (AND, F1, F2) o F = (OR, F1, F2)
	formula_lqdag(uint8_t functor, formula_lqdag *l1, formula_lqdag *l2) {
		assert(functor == FUNCTOR_AND || functor == FUNCTOR_OR);
		this->functor = functor;
		this->formula_lqdag1 = l1;
		this->formula_lqdag2 = l2;
		this->attribute_set_A = {};
//		this->formula_leaf = nullptr;
	}

	// L = (EXTEND, L1, A)
	formula_lqdag(uint8_t functor, formula_lqdag* l, att_set &attribute_set_A, type_mapping_M *M = nullptr) {
		assert(functor == FUNCTOR_EXTEND);
		this->functor = functor;
		this->formula_lqdag1 = l;
		this->formula_lqdag2 = nullptr;
		this->attribute_set_A = attribute_set_A;
//		this->formula_leaf = nullptr;

		if(M == nullptr) {
			uint64_t msize = std::pow(this->lqdag1->subQuadtree->qdag->getK(),
									  this->lqdag1->subQuadtree->qdag->getD());

			uint16_t dim = this->attribute_set_A.size(); // d
			uint16_t dim_prime = this->lqdag1->subQuadtree->qdag->nAttr(); // d'
			uint64_t p = std::pow(this->lqdag1->subQuadtree->qdag->getK(), dim);

			M = new type_mapping_M[p];

			uint64_t mask;
			uint64_t i, i_prime;

			for (i = 0; i < p; ++i) {
				// todos los bits están en cero excepto el bit en la posición dim_prime - 1.
				mask = 1 << (dim_prime - 1); // equivalent to 2^(dim_prime-1)
				i_prime = 0;

				for (uint16_t j = 0; j < dim_prime; ++j) {
					if (i & (1 << (dim - this->lqdag1->subQuadtree->qdag->getAttr(j) - 1)))
						i_prime |= mask;

					mask >>= 1;
				}

				M[i] = i_prime; // = M[i_prime];
			}
		}
		this->M = M;

		if (attribute_set_A.size() == 3)
			this->lqdag1->subQuadtree->qdag->createTableExtend3();
		else if (attribute_set_A.size() == 4)
			this->lqdag1->subQuadtree->qdag->createTableExtend4();
		else if (attribute_set_A.size() == 5)
			this->lqdag1->subQuadtree->qdag->createTableExtend5();
		else {
			cout << "EXTEND only works for queries of up to 5 attributes..." << endl;
			exit(1);
		}
	}

	uint8_t get_functor() const {
		return this->functor;
	}

};

#endif //QDAGS_FORMULA_HPP
