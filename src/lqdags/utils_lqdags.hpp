//
// Created by Asunción Gómez on 19-06-24.
//

#ifndef QDAGS_UTILS_LQDAGS_HPP
#define QDAGS_UTILS_LQDAGS_HPP

#include<bits/stdc++.h>
#include "../qdags.hpp"
#include "../lqdag.hpp"

const double NO_VALUE_LEAF_UTILS = 3;

#define OP_EQUAL 0
#define OP_UNEQUAL 1
#define OP_LESS_EQUAL 2
#define OP_GREATER_EQUAL 3
#define OP_LESS 4
#define OP_GREATER 5

const uint8_t TYPE_ATT1_ATT2 = 0;
const uint8_t TYPE_ATT1_CONST = 1;
const uint8_t TYPE_ATT2_CONST = 2;

struct quadtree_formula { // represents the output of the formula we are evaluating
    double val_leaf = NO_VALUE_LEAF_UTILS; // EMPTY_LEAF, FULL_LEAF, INTERNAL_NODE, NO_VALUE_LEAF
    quadtree_formula** children = nullptr;
};

// represents a subtree of the quadtree
struct subQuadtreeChild {
    ::qdag *qdag;
    uint16_t level; // the level of the quadtree_formula. 0 to start.
    uint64_t node; // absolute position in the bv[level] of the quadtree
    double value;// = NO_VALUE_LEAF;
    // Constructor sin inicializar node_description
    subQuadtreeChild(::qdag* qdagPtr, uint16_t lvl, uint64_t nd, double val = NO_VALUE_LEAF_UTILS)
            : qdag(qdagPtr), level(lvl), node(nd), value(val) {}
};

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
static double eval_pred(predicate* pred, uint16_t* coordinates, uint64_t quadrant_side, uint64_t nAttr) {
    vector<uint64_t> min_att;
    vector<uint64_t> max_att;
    for(uint16_t i = 0; i < nAttr; i++){
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
                if(quadrant_side == 1){
                    return min_att_1 <= pred->val_const;
                }
                if(max_att_1 <= pred->val_const)
                    return 1;
                else if(min_att_1 <= pred->val_const && pred->val_const < max_att_1)
                    return 0.5;
                else if(min_att_1 > pred->val_const)
                    return 0;
                else
                    throw "error: eval_pred OP_LESS_EQUAL"; // no debiese llegar aqui
            case OP_LESS: // att1 < val_const
                if(quadrant_side == 1){
                    return min_att_1 < pred->val_const;
                }
                if(max_att_1 < pred->val_const)
                    return 1;
                else if(min_att_1 <= pred->val_const && pred->val_const <= max_att_1)
                    return 0.5;
                else if(min_att_1 > pred->val_const)
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
        switch (pred->op) {
            min_att_2 = min_att[pred->att_2];
            max_att_2 = max_att[pred->att_2];
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
                if(quadrant_side == 1){
                    return min_att_2 <= pred->val_const;
                }
                if(max_att_2 <= pred->val_const)
                    return 1;
                else if(min_att_2 <= pred->val_const && pred->val_const < max_att_2)
                    return 0.5;
                else if(min_att_2 > pred->val_const)
                    return 0;
                else
                    throw "error: eval_pred OP_LESS_EQUAL"; // no debiese llegar aqui
            case OP_LESS: // att2 < val_const
                if(quadrant_side == 1){
                    return min_att_2 < pred->val_const;
                }
                if(max_att_2 < pred->val_const)
                    return 1;
                else if(min_att_2 <= pred->val_const && pred->val_const <= max_att_2)
                    return 0.5;
                else if(min_att_2 > pred->val_const)
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


struct quadtree_pred{
    ::predicate *pred;
    uint16_t *coordinates;
    uint16_t cur_level;
    uint16_t max_level;
    uint64_t grid_side;
    uint8_t k;
    uint64_t nAttr;
};



#endif //QDAGS_UTILS_LQDAGS_HPP
