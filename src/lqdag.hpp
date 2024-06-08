//
// We will not use qdags, only quadtrees
//
#ifndef INCLUDED_LQDAGS
#define INCLUDED_LQDAGS

#include<bits/stdc++.h>
#include "./qdags.hpp"
#include "./utils.hpp"

const uint8_t FUNCTOR_QTREE = 0; // leaf
const uint8_t FUNCTOR_NOT = 1; // leaf
const uint8_t FUNCTOR_AND = 2; // internal quadtree_formula
const uint8_t FUNCTOR_OR = 3; // internal quadtree_formula
const uint8_t FUNCTOR_EXTEND = 4; // internal quadtree_formula

const double NO_VALUE_LEAF = 3;
const double EMPTY_LEAF = 0;
const double FULL_LEAF = 1;
const double INTERNAL_NODE = 0.5;
// This indicates that one cannot determine the value of the quadtree_formula without computing the values of its children.
const double VALUE_NEED_CHILDREN = 2; // it indicates we need to compute the values of its children

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
    double val_leaf = NO_VALUE_LEAF; // EMPTY_LEAF, FULL_LEAF, INTERNAL_NODE, NO_VALUE_LEAF
    quadtree_formula** children = nullptr;
};

// represents a subtree of the quadtree
struct subQuadtreeChild {
    ::qdag *qdag;
    uint16_t level; // the level of the quadtree_formula. 0 to start.
    uint64_t node; // absolute position in the bv[level] of the quadtree
    double value;// = NO_VALUE_LEAF;
    // Constructor sin inicializar node_description
    subQuadtreeChild(::qdag* qdagPtr, uint16_t lvl, uint64_t nd, double val = NO_VALUE_LEAF)
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

// lqdag as a syntax tree
class lqdag {

public:
    typedef vector<uint64_t> att_set;

    /**
     *
     * @param pred
     * @param coordinates the initial position x,y,z,... of the quadrant
     * @param quadrant_side
     * @return
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

private:
    // lqdag L=(f,o), where f is a functor.
    uint8_t functor; // QTREE, NOT, AND, OR, EXTEND
    subQuadtreeChild *subQuadtree;
    lqdag *lqdag1;
    lqdag *lqdag2;
    double val_lqdag1 = NO_VALUE_LEAF;
    double val_lqdag2 = NO_VALUE_LEAF;
    att_set attribute_set_A; // for extend functor
    type_mapping_M *M; // mapping for extend functor. Use it when we need to extend the subquadtree, Otherwise M[i] = i (case not extended).

public:
    lqdag() = default;

    // L = (QTREE, Q_r)
    lqdag(subQuadtreeChild* subQuadtree) {
        this->functor = FUNCTOR_QTREE;
        this->subQuadtree = subQuadtree;
        this->lqdag1 = nullptr;
        this->lqdag2 = nullptr;
//        for(uint64_t i = 0; i < subQuadtree->qdag->nAttr(); i++) // TODO: pass att set as param
//            this->attribute_set_A.push_back(subQuadtree->qdag->getAttr(i));
    }

    // L = (QTREE, Q_r) o L = (NOT, Q_r)
    lqdag(uint8_t functor, subQuadtreeChild* subQuadtree) {
        assert(functor == FUNCTOR_QTREE || functor == FUNCTOR_NOT);
        this->functor = functor;
        this->subQuadtree = subQuadtree;
        this->lqdag1 = nullptr;
        this->lqdag2 = nullptr;
//        for(uint64_t i = 0; i < subQuadtree->qdag->nAttr(); i++) // TODO: pass att set as param
//            this->attribute_set_A.push_back(subQuadtree->qdag->getAttr(i));
    }

    // L = (AND, L1, L2) o L = (OR, L1, L2)
    lqdag(uint8_t functor, lqdag *l1, lqdag *l2) {
        assert(functor == FUNCTOR_AND || functor == FUNCTOR_OR);
        this->functor = functor;
        this->lqdag1 = l1;
        this->lqdag2 = l2;
        this->subQuadtree = nullptr;
        this->attribute_set_A = {};
    }

    // L = (EXTEND, L1, A)
    // TODO: puedo asumir q dps delextend viene un subquadtreE?
    lqdag(uint8_t functor, lqdag* l, att_set &attribute_set_A, type_mapping_M *M = nullptr) {
        assert(functor == FUNCTOR_EXTEND);
        this->functor = functor;
        this->lqdag1 = l; // has
        this->attribute_set_A = attribute_set_A;
        this->subQuadtree = nullptr;
        this->lqdag2 = nullptr;

        if(M == nullptr) {
            uint64_t msize = std::pow(this->lqdag1->subQuadtree->qdag->getK(),
                                        this->lqdag1->subQuadtree->qdag->getD());
            this->M = new type_mapping_M[msize];

            uint16_t dim = this->attribute_set_A.size(); // d
            uint16_t dim_prime = this->lqdag1->attribute_set_A.size(); // d'
            uint64_t p = std::pow(this->lqdag1->subQuadtree->qdag->getK(), dim);

            M = new type_mapping_M[p];

            uint64_t mask;
            uint64_t i, i_prime;

            for (i = 0; i < p; ++i) {
                // todos los bits están en cero excepto el bit en la posición dim_prime - 1.
                mask = 1 << (dim_prime - 1); // equivalent to 2^(dim_prime-1)
                i_prime = 0;

                for (uint16_t j = 0; j < dim_prime; ++j) {
                    if (i & (1 << (dim - this->lqdag1->attribute_set_A[j] - 1)))
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

    uint8_t getFunctor(){
        return this->functor;
    }

    att_set getAttributeSet(){
        return this->attribute_set_A;
    }

    type_mapping_M* getM(){
        return this->M;
    }

    subQuadtreeChild* getSubQuadtree(){
        return this->subQuadtree;
    }

    // algorithm 8 child(L,i)
    /**
     * Can only be invoked when value is 1/2 or VALUE_NEED_CHILDREN.
     * @param i
     * @return
     */
    lqdag* get_child_lqdag(uint64_t i){
        // case base
        switch (this->functor) {
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
                lqdag *l = new lqdag(this->functor, subQuadtree);
                return l;
            }
            case FUNCTOR_AND: {
                double val_l1 = this->val_lqdag1 != NO_VALUE_LEAF ? this->val_lqdag1 : this->lqdag1->value_lqdag();
                if (val_l1== 1)
                    return lqdag2->get_child_lqdag(i);
                double val_l2 = this->val_lqdag2 != NO_VALUE_LEAF ? this->val_lqdag2 : this->lqdag2->value_lqdag();
                if (val_l2 == 1)
                    return lqdag1->get_child_lqdag(i);
                lqdag *l = new lqdag(FUNCTOR_AND, this->lqdag1->get_child_lqdag(i), this->lqdag2->get_child_lqdag(i));
                return l;
            }
            case FUNCTOR_OR: {
                double val_l1 = this->val_lqdag1 != NO_VALUE_LEAF ? this->val_lqdag1 : this->lqdag1->value_lqdag();
                if (val_l1 == 0)
                    return lqdag2->get_child_lqdag(i);
                double val_l2 = this->val_lqdag2 != NO_VALUE_LEAF ? this->val_lqdag2 : this->lqdag2->value_lqdag();
                if (val_l2 == 0)
                    return lqdag1->get_child_lqdag(i);
                lqdag *l = new lqdag(FUNCTOR_OR, this->lqdag1->get_child_lqdag(i), this->lqdag2->get_child_lqdag(i));
                return l;
            }
            case FUNCTOR_EXTEND: { // Extend quadtree to attributes att_set
                // i --> i' (hago el mapeo)
                lqdag *l = new lqdag(FUNCTOR_EXTEND, this->lqdag1->get_child_lqdag(this->M[i]), this->attribute_set_A, this->M);
                return l;
            }
            default:
                throw "error: value_lqdag non valid functor";

        }
    }

    /**
     * algorithm 7 value(L).
     * In lqdags we introduce a new idea: full leaves, that denote subgrids full of 1s.
     * We save the value in val_lqdag_1, val_lqdag_2 or subQuatree->value. If we've already computed the value, we return it.
     * @return the value of the root of the lqdag. Can be EMPTY_LEAF, FULL_LEAF, INTERNAL_NODE or VALUE_NEED_CHILDREN.
     */
    double value_lqdag(){
        switch(this->functor){
            case FUNCTOR_QTREE: {
                return this->value_quadtree();
            }
            case FUNCTOR_NOT: {
                return 1 - this->value_quadtree();
            }
            case FUNCTOR_AND: {
                if(this->val_lqdag1 == NO_VALUE_LEAF)
                    this->val_lqdag1 = this->lqdag1->value_lqdag();
                if(this->val_lqdag1 == 0)
                    return EMPTY_LEAF;

                if(this->val_lqdag2 == NO_VALUE_LEAF)
                    this->val_lqdag2 = this->lqdag2->value_lqdag();
                if(this->val_lqdag2 == 0)
                    return EMPTY_LEAF;

                if (this->val_lqdag1 == 1)
                    return this->val_lqdag2;
                if (this->val_lqdag2 == 1)
                    return this->val_lqdag1;

                return VALUE_NEED_CHILDREN; // if both roots have Value 1⁄2, one cannot be sure of the Value of the resulting root until the AND between the children of Q1 and Q2 has been computed
            }
            case FUNCTOR_OR: {
                if(this->val_lqdag1 == NO_VALUE_LEAF)
                    this->val_lqdag1 = this->lqdag1->value_lqdag();
                if(this->val_lqdag1 == 1)
                    return FULL_LEAF;

                if(this->val_lqdag2 == NO_VALUE_LEAF)
                    this->val_lqdag2 = this->lqdag2->value_lqdag();
                if(this->val_lqdag2 == 1)
                    return FULL_LEAF;

                if (this->val_lqdag1 == 0)
                    return this->val_lqdag2;
                if (this->val_lqdag2 == 0)
                    return this->val_lqdag1;

                return VALUE_NEED_CHILDREN; // if both roots have Value 1⁄2, one cannot be sure of the Value of the resulting root until the OR between the children of Q1 and Q2 has been computed
            }
            case FUNCTOR_EXTEND: {
                if (this->val_lqdag1 == NO_VALUE_LEAF)
                    this->val_lqdag1 = this->lqdag1->value_lqdag();
                return this->val_lqdag1;
            }
            default:
                throw "error: value_lqdag non valid functor";
        }
    }

    /**
     * Algorithm 1 Value(Q), adapted for lqdag
     * Value of a qdag
     * In lqdags we introduce a new idea: full leaves, that denote subgrids full of 1s.
     * Leaves can be in a higher level than the last one, when the subgrids are all 0s or all 1s.
     * FULL_LEAF if the qdag represents a full single cell, EMPTY_LEAF if it is an empty leaf, INTERNAL_NODE if is an internal quadtree_formula.
     * @return EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE.
     */
    double value_quadtree(){
        if(this->functor == FUNCTOR_QTREE) { // only is called
            uint16_t cur_level = this->subQuadtree->level;
            uint64_t cur_node = this->subQuadtree->node;
            if(this->subQuadtree->value == NO_VALUE_LEAF)
                this->subQuadtree->value = is_leaf(cur_level, cur_node);
            return this->subQuadtree->value;

        } else {
            //cout << "error: value_quadtree called on non-quadtree" << endl;
            throw "error: value_quadtree called on non-quadtree";
        }
    }

    /**
     *
     * @param level
     * @param node
     * @param k
     * @param k_d
     * If level = 0 and node = 0, we are looking the root of the quadtree.
     * If level > max_level, we are looking for the bit representing a point in the grid.
     * @return whether is an EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE
     */
    double is_leaf(uint16_t level, uint64_t node){
        uint16_t max_level = this->subQuadtree->qdag->getHeight() - 1;
        if(level > max_level){
            return this->subQuadtree->qdag->get_ith_bit(max_level, node);
        }
        uint8_t node_description = this->subQuadtree->qdag->get_node_last_level(level, node);
        if((node_description | 0) == 0)
            return EMPTY_LEAF;
        uint16_t val_full_ones;
        uint64_t k_d = this->subQuadtree->qdag->getKD();
        switch (k_d) {
            case 2:
                val_full_ones = 3;
                break;
            case 4:
                val_full_ones = 15;
                break;
            default:
                throw "error: invalid k_d";
        }
        if((node_description & val_full_ones) == val_full_ones){ // check the children
            if(level == max_level)
                return FULL_LEAF;
            uint64_t siblings = this->subQuadtree->qdag->rank(level,node);
            for(uint64_t i = 0; i < k_d; i++){
                if(is_leaf(level+1, (siblings+i)*k_d) != 1)
                    return INTERNAL_NODE;
            }
            return FULL_LEAF;
        }
        return INTERNAL_NODE;
    }



    /**
     * Create an empty or full leaf.
     * @param val can be EMPTY_LEAF or FULL_LEAF.
     * @return quadtree_formula with the value of the leaf.
     */
    quadtree_formula* create_leaf(double val){
        quadtree_formula* newNode = new quadtree_formula;
        newNode->val_leaf = val;
        newNode->children = nullptr;
        return newNode;
    }

    /**
     * Algorithm 9 of paper.
     * The completion Q_f is the quadtree representing the output of the formula F, represented as an lqdag.
     *
     * @param p
     * @param max_level
     * @param cur_level the current level of the quadtree. If cur_level = max_level + 1, we are looking for the bit representing the point in the grid.
     * @param UPPER_BOUND
     * @param results
     * @return
     */
    quadtree_formula* completion(uint64_t p,
                                 uint16_t max_level,
                                 uint16_t cur_level,
                                 uint64_t UPPER_BOUND,
                                 uint64_t &results) {

        if(results >= UPPER_BOUND){ // top k results only
            quadtree_formula* newNode = create_leaf(NO_VALUE_LEAF);
            return newNode;
        }
        double val_lqdag = this->value_lqdag();
        if(val_lqdag == EMPTY_LEAF || val_lqdag == FULL_LEAF){
            quadtree_formula* newNode = create_leaf(val_lqdag);
            if(val_lqdag == FULL_LEAF){
                results += (cur_level > max_level) ? 1 :  pow(p, (max_level - cur_level + 1));
            }
            return newNode;
        }

        double max_value = 0;
        double min_value = 1;
        quadtree_formula **Q_f_children = new quadtree_formula*[p];
        // if all quadtres are empty, return a 0 leaf.
        // C_i ← Completion(Child(L_F,i))
        for(uint64_t i = 0; i < p; i++){
            Q_f_children[i] = this->get_child_lqdag(i)->completion(p, max_level, cur_level+1, UPPER_BOUND, results);
            max_value = max(max_value, Q_f_children[i]->val_leaf); // val_leaf can be EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE
            min_value = min(min_value, Q_f_children[i]->val_leaf);
        }
        // return a leaf
        if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){ // return a leaf
            for(uint64_t i = 0; i < p; i++)
                delete Q_f_children[i];
            delete[] Q_f_children;
            double val_leaf = (max_value == EMPTY_LEAF) ? EMPTY_LEAF : FULL_LEAF;
            quadtree_formula* newNode = create_leaf(val_leaf);
            return newNode;
        }
        else {
            // return an internal node with children
            quadtree_formula* quadtree = new quadtree_formula();
            quadtree->val_leaf = val_lqdag; //INTERNAL_NODE;
            quadtree->children = Q_f_children;
            return quadtree;
        }
    }


    /**
     *
     * @param p
     * @param max_level
     * @param cur_level
     * @param UPPER_BOUND
     * @param results
     * @param pred
     * @param coordinates an array of nAtt coordinates of the node representing the quadrant.
     * @param checkPred once we have evaluated the predicate and if this evaluation returns 1, we do not have to evaluate it again for its children.
     * @return
     */
    quadtree_formula *
    completion_with_pred(uint64_t p, uint16_t max_level, uint16_t cur_level, uint64_t grid_side, uint8_t k, uint64_t UPPER_BOUND, uint64_t &results,
                         predicate *pred, uint16_t *coordinates, uint64_t nAttr) {

        if(results >= UPPER_BOUND){ // top k results only
            quadtree_formula* newNode = create_leaf(EMPTY_LEAF);
            return newNode;
        }
        // if results < UPPER_BOUND
        uint64_t quadrant_side = grid_side/pow(k,cur_level);
        double val_eval_pred = eval_pred(pred, coordinates, quadrant_side, nAttr);
        if(val_eval_pred == 0){
            quadtree_formula* newNode = create_leaf(EMPTY_LEAF);
            return newNode;
        } else if(val_eval_pred == 1){
            return this->completion(p, max_level, cur_level, UPPER_BOUND, results);
        } else {
            double max_value = 0;
            double min_value = 1;
            quadtree_formula **Q_f_children = new quadtree_formula*[p];
            // if all quadtres are empty, return a 0 leaf.
            // C_i ← Completion(Child(L_F,i))
            for(uint64_t i = 0; i < p; i++){
                uint16_t diff_level = max_level-cur_level;
                uint16_t l = (uint16_t) log2(p);
                uint16_t* coordinatesTemp = new uint16_t[l];
                for(uint16_t j = 0; j < l; j++)
                    coordinatesTemp[j] = coordinates[j];
                transformCoordinates(coordinatesTemp, l, diff_level, i);
                quadtree_formula* newNode = this->get_child_lqdag(i)->completion_with_pred(p, max_level, cur_level + 1,
                                                                                           grid_side, k,
                                                                                           UPPER_BOUND, results, pred,
                                                                                           coordinatesTemp, nAttr);
                Q_f_children[i] = newNode;
                max_value = max(max_value, newNode->val_leaf); // val_leaf can be EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE
                min_value = min(min_value, newNode->val_leaf);
            }
            // return a leaf
            if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){ // return a leaf
                for(uint64_t i = 0; i < p; i++)
                    delete Q_f_children[i];
                delete[] Q_f_children;
                double val_leaf = (max_value == EMPTY_LEAF)? EMPTY_LEAF : FULL_LEAF;
                quadtree_formula* newNode = create_leaf(val_leaf);
                return newNode;
            }
            else {
                // return an internal node with children
                quadtree_formula* quadtree = new quadtree_formula;
                quadtree->val_leaf = INTERNAL_NODE;
                quadtree->children = Q_f_children;
                return quadtree;
            }

        }

    }


    /**
     * It builds the quadtree_formula Q_f of the lqdag L_F, as we ask for a particular node.
     * @param p
     * @param max_level
     * @param cur_level
     * @param UPPER_BOUND
     * @param results
     * @param level
     * @param child the i-th child of the current node
     * @param Q_f the quadtree_formula built so far.
     * @return a pointer to the child of the lqdag
     */
    lqdag* lazy_child_completion(uint64_t p,
                                   uint64_t child,
                                   quadtree_formula& Q_f){
        if(Q_f.val_leaf == EMPTY_LEAF || Q_f.val_leaf == FULL_LEAF){
            return this; // child is 0 or 1 and the quadtree is already computed
        }
        else if(Q_f.val_leaf == INTERNAL_NODE || Q_f.val_leaf == VALUE_NEED_CHILDREN){
            lqdag* lqdag_child = this->get_child_lqdag(child);
            if(!Q_f.children[child]){ // allocate memory for the child
                Q_f.children[child] = new quadtree_formula{};
            } // if it's not computed
            if(Q_f.children[child]->val_leaf == NO_VALUE_LEAF){
                quadtree_formula* qf_no_value = lqdag_child->compute_quadtree_formula_node(p);
                Q_f.children[child]->val_leaf = qf_no_value->val_leaf;
                Q_f.children[child]->children = qf_no_value->children;
            } // else: it's already computed
            return lqdag_child;
        }
        else{ // NO_VALUE --> compute the root and then, the child
            quadtree_formula* qf_no_value = this->compute_quadtree_formula_node(p);
            Q_f.val_leaf = qf_no_value->val_leaf;
            Q_f.children = qf_no_value->children;
            return this->lazy_child_completion(p, child, Q_f);
        }
    }

    /**
     * Compute the value of a node of the quadtree_formula Q_f of the lqdag.
     * @return
     */
    quadtree_formula* compute_quadtree_formula_node(uint64_t p){
        double val_lqdag = this->value_lqdag();
        if(val_lqdag == EMPTY_LEAF || val_lqdag == FULL_LEAF){
            return create_leaf(val_lqdag);
        } else{
            quadtree_formula* newNode = new quadtree_formula;
            newNode->val_leaf = val_lqdag;
            newNode->children = new quadtree_formula*[p];
            return newNode;
        }
    }
};

#endif