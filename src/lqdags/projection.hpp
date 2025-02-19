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

	// TODOO: maybe replace it, idk
	struct NodeMaterialization {
		double value = NO_VALUE_LEAF; // the evaluation
		uint64_t n_children;      // Number of children
		uint64_t n_projection;      // Number of children
		NodeMaterialization **children; // Dynamic array of child pointers
		// the projection to evaluate and compute value and the values of the n_children
		projection** projection_children = nullptr; // it has 2^|A| children (leaves of projection)

		// Constructor
		NodeMaterialization(uint64_t n, uint64_t p) : value(NO_VALUE_LEAF), n_children(n), n_projection(p), children(nullptr), projection_children(nullptr) {}

		NodeMaterialization(double val, uint64_t n, uint64_t p) : value(val), n_children(n), n_projection(p), children(nullptr), projection_children(nullptr) {}

		void init_children(){
			children = new NodeMaterialization*[n_children];
			for(uint64_t j = 0; j < n_children; j++){
				children[j] = nullptr;
			}
		}

		void init_projection_children(){
			projection_children = new projection*[n_projection];
			for(uint64_t j = 0; j < n_projection; j++){
				projection_children[j] = nullptr;
			}
		}

		double get_val_node(){
			return this->value;
		}

		void set_val_node(double val){
			this->value = val;
		}

		double get_value_child_node(uint64_t i){
			assert(i < n_children);
			if(children != nullptr && children[i] != nullptr)
				return children[i]->value;
			return NO_VALUE_LEAF;
		}

		projection* get_ith_projection(uint64_t i){
			assert(i < n_projection);
			if(projection_children != nullptr && projection_children[i] != nullptr)
				return projection_children[i];
			return nullptr;
		}

		void set_ith_projection(uint64_t i, projection* proj){
			assert(i < n_projection);
			projection_children[i] = proj;
		}

		// Add a child or change the value of a child
		void addChild(uint64_t i, double val_children_i){
			if(children == nullptr){
				init_children();
				init_projection_children();
			}
			if(children[i] == nullptr){
				children[i] = new NodeMaterialization(val_children_i,n_children, n_projection);
			}
			children[i]->value = val_children_i;
		}


		// Destructor to free memory
		~NodeMaterialization() {
			if(children != nullptr){
				for (uint64_t i = 0; i < n_children; ++i) {
					if(children[i] != nullptr)
						delete children[i];
				}
				delete[] children;
			}
			if(projection_children != nullptr){
				for (uint64_t i = 0; i < n_projection; ++i) {
					if(projection_children[i] != nullptr)
						delete projection_children[i];
				}
				delete[] projection_children;
			}
		}
	};


private:
//	formula_pi* form_or; // expression, syntax tree, formula of the OR tree.
	att_set attribute_set_A; // A
	att_set attribute_set_A_prime; // A'
	lqdag* lqdag_root; // the lqdag where we apply the projection
//	projection** projection_children; // it has 2^|A| children (leaves of projection)
	indexes indexes_children; // indexes of the j for each i
	NodeMaterialization* node_materialization; // the materialization of the projection

public:

	projection() = default;

	projection(lqdag* lqdag_root, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
		this->lqdag_root = lqdag_root;
		this->node_materialization = new NodeMaterialization(nChildren(), 1<<attribute_set_A.size());

		// TODO: see if memory is correct allocated!
//		this->form_or = new formula_pi(attribute_set_A, attribute_set_A_prime, lqdag_root);
		this->indexes_children = create_pi_index_children(attribute_set_A, attribute_set_A_prime);

		// this->indexes_children = create_pi_index_children(attribute_set_A, attribute_set_A_prime);

		uint64_t number_projection_children = 1 << attribute_set_A.size();

//		projection_children = new projection*[number_projection_children];
//		for(uint64_t i = 0; i < number_projection_children; i++){ // number of leaves
//			projection_children[i] = nullptr;
//		}

		// TODO: OFT
		projection* test = this->eval_projection(0, this->node_materialization);
		// TODO: see materialization
		test->get_val_node();

	}

	projection(lqdag* lqdag_root, indexes &ind, att_set attribute_set_A, att_set attribute_set_A_prime){
		this->attribute_set_A = attribute_set_A;
		this->attribute_set_A_prime = attribute_set_A_prime;
		this->indexes_children = ind;
		this->lqdag_root = lqdag_root;
		this->node_materialization = new NodeMaterialization(nChildren(), 1 << attribute_set_A.size());

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
		return this->node_materialization->get_ith_projection(i);
	}

	void set_ith_projection(uint64_t i, projection* proj){
		assert(i < 1 << nAttrA());
		this->node_materialization->set_ith_projection(i, proj);
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


	double get_val_node(){
		return this->node_materialization->get_val_node();
	}

	void set_val_node(double val){
		this->node_materialization->set_val_node(val);
	}

	/**
	 * Get the value of the i-th child of the materialization quadtree
	 * @param i the i-th child.
	 * @return
	 */
	double get_value_child_node(uint64_t i){
		assert(i < this->nChildren());
		if(this->node_materialization->children != nullptr
		   && this->node_materialization->children[i] != nullptr) {
			return this->node_materialization->children[i]->value;
		}
		return NO_VALUE_LEAF; // TODO: check its return this!
	}

	/**
	 * Set the value for the i-th child of the materialization quadtree
	 * @param i the i-th child
	 * @param val the value of the i-th child
	 */
	void set_value_child_node(uint64_t i, double val){
		assert(i < this->nChildren());
		this->node_materialization->addChild(i, val);
	}

	/**
	 * Get the value of the current projection (EMPTY_LEAF, FULL_LEAF or VALUE_NEED_CHILDREN)
	 * if the value of the original lqdag is 0 or 1, return that value. Otherwise, return VALUE_NEED_CHILDREN
	 * @return
	 */
	double value_pi(){
		double val_lqdag = this->lqdag_root->value_lqdag(lqdag_root->get_formula());
		if(val_lqdag == FULL_LEAF || val_lqdag == EMPTY_LEAF){
			return val_lqdag;
		}
		return VALUE_NEED_CHILDREN; // initially it will return NO_VALUE_LEAF, and then VALUE_NEED_CHILDREN
	}


	/**
	 * Compute a multi-OR between the lqdags chidlren of the i-th projection
	 * @param i the i-th child of the projection.
	 * @return EMPTY_LEAF, FULL_LEAF or VALUE_NEED_CHILDREN
	 */
	double get_value_multi_or_child_pi(uint64_t i){
		assert(i < this->nChildren());

		if(this->get_value_child_node(i) != NO_VALUE_LEAF){
			return this->get_value_child_node(i);
		}

		double max_value = 0;
		uint64_t n_leaves_per_children = 1 << (nAttrA() - nAttrA_prime());
		uint64_t init_index = i * n_leaves_per_children;
		for(uint64_t j = init_index; j < init_index + n_leaves_per_children; j++) {
			lqdag *lqdag_j = this->lqdag_root->get_child_lqdag(this->get_index_children(j));
			double val_lqdag_j = lqdag_j->value_lqdag(lqdag_j->get_formula());
			if(val_lqdag_j == FULL_LEAF){
//				this->set_value_child_node(i, FULL_LEAF);
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
//		this->set_value_child_node(i, val_or);
		return val_or;

	}

	void compute_multi_or_projection_child_pi(uint64_t i){
		assert(i < this->nChildren());
		assert(this->get_value_child_node(i) == VALUE_NEED_CHILDREN);

		uint64_t n_leaves_per_children = 1 << (nAttrA() - nAttrA_prime());
		uint64_t init_index = i * n_leaves_per_children;
		for(uint64_t j = init_index; j < init_index + n_leaves_per_children; j++) {
			lqdag *lqdag_j = this->lqdag_root->get_child_lqdag(this->get_index_children(j));
			double val_lqdag_j = lqdag_j->value_lqdag(lqdag_j->get_formula());
			if(val_lqdag_j == EMPTY_LEAF){
				this->set_ith_projection(j, nullptr);
			} else {
				projection* child_pi = new projection(lqdag_j, this->indexes_children, this->attribute_set_A, this->attribute_set_A_prime);
				child_pi->set_val_node(val_lqdag_j);
				this->set_ith_projection(j,child_pi);
			}
		}

	}


	/**
	 *
	 * @return a projection with the output as a traditional quadtree (with the values in value_children)
	 */
	projection* eval_projection(uint16_t cur_level, NodeMaterialization* materialization, projection* last_proj = nullptr){
		this->set_val_node(this->value_pi());

		if(this->value_pi() == FULL_LEAF || this->value_pi() == EMPTY_LEAF){
			return this;
		}
		double max_value = 0;
		double min_value = 1;

		// we evaluate the previous projection and we keep the best value (because it's an OR)
		for(uint64_t i = 0; i < nChildren(); i++){
			double parent_child_i = last_proj ? last_proj->get_value_child_node(i) : 3;
			double val_child_pi = get_value_multi_or_child_pi(i); // evaluate an OR
			this->set_value_child_node(i, val_child_pi);

			if(val_child_pi == FULL_LEAF){
				if(parent_child_i != FULL_LEAF){
					materialization->addChild(i, val_child_pi);
				}
			}
			else if(val_child_pi == VALUE_NEED_CHILDREN){
				if(parent_child_i == EMPTY_LEAF || parent_child_i == NO_VALUE_LEAF){
					materialization->addChild(i, val_child_pi);
				}
				compute_multi_or_projection_child_pi(i); // compute the projection leaves for the i-th child
				uint64_t n_leaves_per_children = 1 << (nAttrA() - nAttrA_prime());
				uint64_t init_index = i * n_leaves_per_children;
				projection* last = nullptr;
				for(uint64_t j = init_index; j < init_index + n_leaves_per_children; j++) {
					if(this->get_ith_projection(j) != nullptr) {
						projection *proj_child_j = this->get_ith_projection(j);
						proj_child_j->eval_projection(++cur_level, materialization->children[i],last);
						last = proj_child_j;
						if(proj_child_j->get_val_node() == FULL_LEAF){
							// no further computation is needed for the other ORs for THAT child
							val_child_pi = FULL_LEAF;
							this->set_value_child_node(i, val_child_pi);
							materialization->addChild(i, val_child_pi);
							break;
						}
					}
				}
			}
			else if(val_child_pi == EMPTY_LEAF){
				if(parent_child_i == NO_VALUE_LEAF){
					materialization->addChild(i, val_child_pi);
				}
			}

			else{
				// is a EMPTY_LEAF or a FULL_LEAF, no further computation is needed
				// TODO: count the results?
			}
			max_value = max(max_value, val_child_pi);
			min_value = min(min_value, val_child_pi);
		}

		// if all the children are empty or full, then is a leaf
		if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){
			double val_node = (max_value == EMPTY_LEAF) ? EMPTY_LEAF : FULL_LEAF;
			this->set_val_node(val_node);
			delete[] node_materialization->children;
			node_materialization->children = nullptr;
			delete[] node_materialization->projection_children;
			node_materialization->projection_children = nullptr;
			delete[] materialization->children;
			materialization->children = nullptr;
			delete[] materialization->projection_children;
			materialization->projection_children = nullptr;

		}
		return this;
	}


};

#endif //INCLUDED_PROJECTION
