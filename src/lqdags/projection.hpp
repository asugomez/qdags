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

	struct TreeNode {
		double value = NO_VALUE_LEAF;
		TreeNode **children; // Dynamic array of child pointers
		uint64_t nChildren;      // Number of children

		// Constructor
		TreeNode(uint64_t n) : value(NO_VALUE_LEAF), children(nullptr), nChildren(n) {}

		TreeNode(double val, uint64_t n) : value(val), children(nullptr), nChildren(n) {}

		double get_val_node(){
			return this->value;
		}

		void set_val_node(double val){
			this->value = val;
		}

		// Add a child or change the value of a child
		void addChild(uint64_t i, double val_children_i){
			if(children == nullptr){
				children = new TreeNode*[nChildren];
				for(uint64_t j = 0; j < nChildren; j++){
					children[j] = nullptr;
				}
			}
			if(children[i] == nullptr){
				children[i] = new TreeNode(val_children_i, nChildren);
			}
			if(children[i] != nullptr && children[i]->value == NO_VALUE_LEAF){
				children[i]->value = val_children_i;
			}
		}



		// Destructor to free memory
		~TreeNode() {
			for (uint64_t i = 0; i < nChildren; ++i) {
				delete children[i];
			}
			delete[] children;
		}
	};


private:
//	formula_pi* form_or; // expression, syntax tree, formula of the OR tree.
	att_set attribute_set_A; // A
	att_set attribute_set_A_prime; // A'
	lqdag* lqdag_root; // the lqdag where we apply the projection
	projection** projection_children; // it has 2^|A| children (leaves of projection)
	indexes indexes_children; // indexes of the j for each i
	TreeNode* quadtree_materialization;

public:

	projection() = default;

	projection(lqdag* lqdag_root, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
		this->lqdag_root = lqdag_root;
		this->quadtree_materialization = new TreeNode(nChildren());

		// TODO: see if memory is correct allocated!
//		this->form_or = new formula_pi(attribute_set_A, attribute_set_A_prime, lqdag_root);
		this->indexes_children = create_pi_index_children(attribute_set_A, attribute_set_A_prime);
		// TODO: maybe we only need the indexes children:
		// this->indexes_children = create_pi_index_children(attribute_set_A, attribute_set_A_prime);

		uint64_t number_projection_children = 1 << attribute_set_A.size();

		projection_children = new projection*[number_projection_children];
		for(uint64_t i = 0; i < number_projection_children; i++){ // number of leaves
			projection_children[i] = nullptr;
		}

		// TODO: OFT
		this->eval_projection();
	}

	// TODO: maybe instead of the form, we need the indexes of the children
	projection(lqdag* lqdag_root, indexes &ind, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
//		this->form_or = form_or;
// TODO: check memory leaks for vector indexes
		this->indexes_children = ind;
		this->lqdag_root = lqdag_root;

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


	projection* get_ith_projection(uint64_t i){
		assert(i < 1 << nAttrA());
		return this->projection_children[i];
	}

	/**
	 * Get the correspondant index of the i-th child of the lqdag.
	 * @param i
	 * @return
	 */
	uint64_t get_index_children(uint64_t i){
		assert(i < (1<<nAttrA()));
		return indexes_children[i];
//		return this->form_or->get_index_children(i);
	}

	void set_val_node(double val){
		this->quadtree_materialization->set_val_node(val);
	}

	double get_val_node(){
		return this->quadtree_materialization->get_val_node();
	}

	/**
	 * Get the value of the current projection (EMPTY_LEAF, FULL_LEAF or VALUE_NEED_CHILDREN)
	 * if the value of the original lqdag is 0 or 1, return that value. Otherwise, return VALUE_NEED_CHILDREN
	 * @return
	 */
	double value_pi(){
		double val_lqdag = this->lqdag_root->val_node();
		if(val_lqdag == FULL_LEAF || val_lqdag == EMPTY_LEAF){
			this->set_val_node(val_lqdag);
		}
		return this->get_val_node(); // initially it will return NO_VALUE_LEAF, and then VALUE_NEED_CHILDREN
	}

	double get_value_children(uint64_t i){
		assert(i < this->nChildren());
		if(this->quadtree_materialization->children != nullptr
			&& this->quadtree_materialization->children[i] != nullptr) {
			return this->quadtree_materialization->children[i]->value;
		}
		return NO_VALUE_LEAF; // TODO: check its return this!
	}

	void set_value_children(uint64_t i, double val){
		assert(i < this->nChildren());
		this->quadtree_materialization->addChild(i, val);
	}



	/**
	 * Compute a multi-OR between the lqdags chidlren of the i-th projection
	 * @param i the i-th child of the projection.
	 * @return EMPTY_LEAF, FULL_LEAF or VALUE_NEED_CHILDREN
	 */
	double get_value_multi_or_child_pi(uint64_t i){
		assert(i < this->nChildren());

		if(this->get_value_children(i) != NO_VALUE_LEAF){
			return this->get_value_children(i);
		}

		double max_value = 0;
		uint64_t n_leaves_per_children = 1 << (nAttrA() - nAttrA_prime());
		uint64_t init_index = i * n_leaves_per_children;
		for(uint64_t j = init_index; j < init_index + n_leaves_per_children; j++) {
			lqdag *lqdag_j = this->lqdag_root->get_child_lqdag(this->get_index_children(j));
			double val_lqdag_j = lqdag_j->value_lqdag(lqdag_j->get_formula());
			if(val_lqdag_j == FULL_LEAF){
				this->set_value_children(i, FULL_LEAF);
				return FULL_LEAF;
			}
			max_value = max(max_value, val_lqdag_j);
		}
		double val_or;
		if(max_value == EMPTY_LEAF){
			val_or = EMPTY_LEAF;
		}
		else{
			val_or = VALUE_NEED_CHILDREN;
		}
		this->set_value_children(i, val_or);
		return val_or;

	}

	void compute_multi_projection_child_pi(uint64_t i){
		assert(i < this->nChildren());
		assert(this->value_children[i] == VALUE_NEED_CHILDREN);

		uint64_t n_leaves_per_children = 1 << (nAttrA() - nAttrA_prime());
		uint64_t init_index = i * n_leaves_per_children;
		for(uint64_t j = init_index; j < init_index + n_leaves_per_children; j++) {
			lqdag *lqdag_j = this->lqdag_root->get_child_lqdag(this->get_index_children(j));
			double val_lqdag_j = lqdag_j->value_lqdag(lqdag_j->get_formula());
			if(val_lqdag_j == EMPTY_LEAF){
				this->projection_children[j] = nullptr;
			} else {
				projection* child_pi = new projection(lqdag_j, this->indexes_children, this->attribute_set_A, this->attribute_set_A_prime);
				this->projection_children[j] = child_pi;

			}
		}

	}



	/**
	 *
	 * @return a projection with the output as a traditional quadtree (with the values in value_children)
	 */
	projection* eval_projection(uint16_t cur_level = 0){
		this->set_val_node(this->value_pi());
		double max_value = 0;
		double min_value = 1;

		for(uint64_t i = 0; i < nChildren(); i++){
			double val_child_pi = get_value_multi_or_child_pi(i); // evaluate an OR
			this->set_value_children(i, val_child_pi);
			if(val_child_pi == VALUE_NEED_CHILDREN){
				compute_multi_projection_child_pi(i);
				uint64_t n_leaves_per_children = 1 << (nAttrA() - nAttrA_prime());
				uint64_t init_index = i * n_leaves_per_children;
				for(uint64_t j = init_index; j < init_index + n_leaves_per_children; j++) {
					if(this->projection_children[j] != nullptr) {
						projection *proj_child_j = this->projection_children[j]->eval_projection(++cur_level);
//						if(proj_child_j->val_pi == FULL_LEAF){
//							val_child_pi = FULL_LEAF;
//							value_children[i] = FULL_LEAF;
//							break; // TODO: check it gets out of the for
//						}
					}

				}
			}
			else{
				// is a EMPTY_LEAF or a FULL_LEAF, no further computation is needed
				// TODO: count the results?
			}
			max_value = max(max_value, val_child_pi);
			min_value = min(min_value, val_child_pi);
		}

		if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){
			delete [] projection_children;
			delete [] quadtree_materialization->children;
			double val_node = (max_value == EMPTY_LEAF) ? EMPTY_LEAF : FULL_LEAF;
			this->set_val_node(val_node);
		}
		return this;
	}


};

#endif //INCLUDED_PROJECTION
