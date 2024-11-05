//
// Created by Asunción Gómez on 19-06-24.
//

#ifndef QDAGS_UTILS_LQDAGS_HPP
#define QDAGS_UTILS_LQDAGS_HPP

#include<bits/stdc++.h>
#include "../lqdag.hpp"
#include "../qdags.hpp"

//const double NO_VALUE_LEAF_UTILS = 3;
const double NO_VALUE_COMPUTED = 3;
const double NO_VALUE_LEAF = 3;
const double EMPTY_LEAF = 0;
const double FULL_LEAF = 1;
const double INTERNAL_NODE = 0.5;
const double VALUE_NEED_CHILDREN = 2; // internal node of the syntax tree

#define OP_EQUAL 0
#define OP_UNEQUAL 1
#define OP_LESS_EQUAL 2
#define OP_GREATER_EQUAL 3
#define OP_LESS 4
#define OP_GREATER 5

const uint8_t TYPE_ATT1_ATT2 = 0;
const uint8_t TYPE_ATT1_CONST = 1;
const uint8_t TYPE_ATT2_CONST = 2;

/**
 * Position of the pointer to the node of a qdag
 */
struct position{
	uint16_t level;
	uint64_t cur_node;
//	uint64_t n_children;
	uint256_t path;
	// default constructor
	position() : level(0), cur_node(0), path(0) {}
	position(uint16_t l, uint64_t n, uint256_t mPath): level(l), cur_node(n), path(mPath) {}

	// Copy constructor (deep copy)
	position(const position& p)
			: level(p.level), cur_node(p.cur_node), path(p.path) {}


	// Copy assignment operator (deep copy)
	position& operator=(const position& p) {
		if (this != &p) {  // Avoid self-assignment
			level = p.level;
			cur_node = p.cur_node;
			path = p.path;
		}
		return *this;
	}


	// Destructor
//	~position() {
//	}

};

//// children computed of the lqdag
//struct node_completion{
//	double val_node = NO_VALUE_LEAF; // EMPTY_LEAF, FULL_LEAF, INTERNAL_NODE, NO_VALUE_LEAF
//	uint64_t n_children; // number of children
//	lqdag** lqdag_children = nullptr;
//
//	// default constructor
//	node_completion() : val_node(NO_VALUE_LEAF), n_children(0), lqdag_children(nullptr) {}
//	// constructor with value and p children
//	node_completion(double val, uint64_t p) : val_node(val), n_children(p), lqdag_children(new lqdag*[p]) {}
//	// destructor
//	~node_completion(){
//		if(lqdag_children != nullptr){
//			delete[] lqdag_children;
//		}
//	}
//};

//struct quadtree_formula { // represents the output of the formula we are evaluating
//    double val_node = NO_VALUE_LEAF_UTILS; // EMPTY_LEAF, FULL_LEAF, INTERNAL_NODE, NO_VALUE_LEAF
//    quadtree_formula** children = nullptr;
//};

// represents a subtree of the quadtree


/**
 * represents a predicate in the query
 * can be:
 * - att_1 op att_2 (0) TYPE_ATT1_ATT2
 * - att_1 op val_const (1) TYPE_ATT1_CONST
 * - att_2 op val_const (2) TYPE_ATT2_CONST
 */
struct predicate {
    uint64_t att_1;
    uint64_t att_2;
    double val_const;
    uint8_t op; // can be OP_EQUAL, OP_UNEQUAL, OP_LESS_EQUAL, OP_GREATER_EQUAL, OP_LESS, OP_GREATER
    uint8_t type_pred; // can be 0, 1, 2
};

/**
 * Value of the predicate
 * @param pred
 * @param coordinates the initial position x,y,z,... of the quadrant
 * @param quadrant_side
 * @return 0 if the quadrant do not satisfy the predicate, 1 if it satisfies the predicate, 0.5 if it is partially inside the limits
 */
static double eval_pred(predicate* pred, uint16_t* coordinates, uint64_t quadrant_side, uint16_t nAtt) {
    vector<uint64_t> min_att;
    vector<uint64_t> max_att;
    for(uint16_t i = 0; i < nAtt; i++){
        min_att.push_back(coordinates[i]);
        max_att.push_back(coordinates[i] + quadrant_side - 1);
    }

    uint64_t min_att_1, max_att_1, min_att_2, max_att_2;
    if(pred->type_pred == TYPE_ATT1_ATT2){ // att1 op att2
        min_att_1 = min_att[pred->att_1];
        max_att_1 = max_att[pred->att_1];
        min_att_2 = min_att[pred->att_2];
        max_att_2 = max_att[pred->att_2];
        switch (pred->op){
            case OP_EQUAL:
                if(quadrant_side == 1){
                    return min_att_1 == min_att_2;
                }
                else if(min_att_1 == min_att_2 || max_att_1 == max_att_2) // || min_att_1 == max_att_2 || max_att_1 == min_att_2)
                    return 0.5;
                else
                    return 0;
            case OP_UNEQUAL:
                if(quadrant_side == 1){
                    return min_att_1 != min_att_2;
                }
                else if(min_att_1 == min_att_2 || min_att_1 == max_att_2 || max_att_1 == min_att_2 || max_att_1 == max_att_2)
                    return 0.5;
                else
                    return 1;
            case OP_LESS_EQUAL:
                if(quadrant_side == 1){
                    return min_att_1 <= min_att_2;
                }
                if(max_att_1 <= min_att_2 && min_att_1 <= max_att_2)
                    return 1;
                else if(min_att_1 <= max_att_2 && max_att_1 >= min_att_2) // a part of the quadrant is inside the limits
                    return 0.5;
                else
                    return 0;
            case OP_LESS:
                if(quadrant_side == 1){
                    return min_att_1 < min_att_2;
                }
                if(max_att_1 < min_att_2 && min_att_1 < max_att_2)
                    return 1;
                else if(min_att_1 < max_att_2 && max_att_1 > min_att_2) // a part of the quadrant is inside the limits
                    return 0.5;
                else
                    return 0;
            case OP_GREATER_EQUAL:
                if(quadrant_side == 1){
                    return min_att_1 >= min_att_2;
                }
                if(min_att_1 >= max_att_2 && max_att_1 >= min_att_2)
                    return 1;
                else if(max_att_1 >= min_att_2 &&  min_att_1 <= max_att_2) // a part of the quadrant is inside the limits
                    return 0.5;
                else
                    return 0;
            case OP_GREATER:
                if(quadrant_side == 1){
                    return min_att_1 > min_att_2;
                }
                if(min_att_1 > max_att_2 && max_att_1 > min_att_2)
                    return 1;
                else if(max_att_1 > min_att_2 &&  min_att_1 < max_att_2) // a part of the quadrant is inside the limits
                    return 0.5;
                else
                    return 0;
            default:
                throw "error: eval_pred default att1 op att2";
        }
    }
    else if(pred->type_pred == TYPE_ATT1_CONST){ // att1 op val_const
        min_att_1 = min_att[pred->att_1];
        max_att_1 = max_att[pred->att_1];
        switch (pred->op) {
            case OP_EQUAL:
                if(quadrant_side == 1){
                    return min_att_1 == pred->val_const;
                } else{
                    if(min_att_1<= pred->val_const && pred->val_const <= max_att_1)
                        return 0.5;
                    else
                        return 0;
                }
            case OP_UNEQUAL:
                if(quadrant_side == 1){
                    return min_att_1 != pred->val_const;
                } else{
                    if(min_att_1<= pred->val_const && pred->val_const <= max_att_1)
                        return 0.5;
                    else
                        return 1;
                }
            case OP_LESS_EQUAL: // att1 <= val_const
                if(quadrant_side == 1)
                    return min_att_1 <= pred->val_const;
                if(max_att_1 <= pred->val_const)
                    return 1;
                else if(min_att_1 <= pred->val_const)
                    return 0.5;
                else if(min_att_1 > pred->val_const)
                    return 0;
                else
                    throw "error: eval_pred OP_LESS_EQUAL"; // no debiese llegar aqui
            case OP_LESS: // att1 < val_const
                if(quadrant_side == 1)
                    return min_att_1 < pred->val_const;
                if(max_att_1 < pred->val_const)
                    return 1;
                else if(min_att_1 < pred->val_const)
                    return 0.5;
                else if(min_att_1 >= pred->val_const)
                    return 0;
                else
                    throw "error: eval_pred OP_LESS";
            case OP_GREATER_EQUAL: // att1 >= val_const
                if(quadrant_side == 1)
                    return min_att_1 >= pred->val_const;
                else if(min_att_1 >= pred->val_const)
                    return 1;
                else if(min_att_1 < pred->val_const && pred->val_const <= max_att_1)
                    return 0.5;
                else if(max_att_1 < pred->val_const)
                    return 0;
                else
                    throw "error: eval_pred OP_GREATER_EQUAL";
            case OP_GREATER: // att1 > val_const
                if(quadrant_side == 1){
                    return min_att_1 > pred->val_const;
                }
                if(min_att_1 > pred->val_const)
                    return 1;
                else if(min_att_1 <= pred->val_const && pred->val_const <= max_att_1) // TPDO: ver estas
                    return 0.5;
                else if(max_att_1 < pred->val_const)
                    return 0;
                else
                    throw "error: eval_pred OP_GREATER";
            default:
                throw "error: eval_pred default att1 op val_const";
        }
    }
    else{ // att2 op val_const TYPE_ATT2_CONST
		min_att_2 = min_att[pred->att_2];
		max_att_2 = max_att[pred->att_2];
        switch (pred->op) {
            case OP_EQUAL:
                if(quadrant_side == 1){
                    return min_att_2 == pred->val_const;
                } else{
                    if(min_att_2<= pred->val_const && pred->val_const <= max_att_2)
                        return 0.5;
                    else
                        return 0;
                }
            case OP_UNEQUAL:
                if(quadrant_side == 1){
                    return min_att_2 != pred->val_const;
                } else{
                    if(min_att_2 <= pred->val_const && pred->val_const <= max_att_2)
                        return 0.5;
                    else
                        return 1;
                }
            case OP_LESS_EQUAL: // att2 <= val_const
                if(quadrant_side == 1)
                    return min_att_2 <= pred->val_const;
                if(max_att_2 <= pred->val_const)
                    return 1;
                else if(min_att_2 <= pred->val_const)
                    return 0.5;
                else if(min_att_2 > pred->val_const)
                    return 0;
                else
                    throw "error: eval_pred OP_LESS_EQUAL"; // no debiese llegar aqui
            case OP_LESS: // att2 < val_const
                if(quadrant_side == 1)
                    return min_att_2 < pred->val_const;
                if(max_att_2 < pred->val_const)
                    return 1;
                else if(min_att_2 < pred->val_const)
                    return 0.5;
                else if(min_att_2 >= pred->val_const)
                    return 0;
                else
                    throw "error: eval_pred OP_LESS";
            case OP_GREATER_EQUAL: // att2 >= val_const
                if(quadrant_side == 1){
                    return min_att_2 >= pred->val_const;
                }
                if(min_att_2 >= pred->val_const)
                    return 1;
                else if(min_att_2 < pred->val_const && pred->val_const <= max_att_2)
                    return 0.5;
                else if(max_att_2 < pred->val_const)
                    return 0;
                else
                    throw "error: eval_pred OP_GREATER_EQUAL";
            case OP_GREATER: // att2 > val_const
                if(quadrant_side == 1){
                    return min_att_2 > pred->val_const;
                }
                if(min_att_2 > pred->val_const)
                    return 1;
                else if(min_att_2 <= pred->val_const && pred->val_const <= max_att_2)
                    return 0.5;
                else if(max_att_2 < pred->val_const)
                    return 0;
                else
                    throw "error: eval_pred OP_GREATER";
            default:
                throw "error: eval_pred default att2 op val_const";
        }
    }

}


//struct quadtree_pred{
//    ::predicate *pred;
//    uint16_t *coordinates;
//    uint16_t cur_level;
//    uint16_t max_level;
//    uint64_t grid_side;
//    uint8_t k;
//    uint16_t nAttr;
//};



#endif //QDAGS_UTILS_LQDAGS_HPP
