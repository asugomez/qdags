//
// We will not use qdags, only quadtrees
//
#ifndef INCLUDED_LQDAGS
#define INCLUDED_LQDAGS

#include<bits/stdc++.h>
#include "./qdags.hpp"

const uint8_t FUNCTOR_QTREE = 0; // leaf
const uint8_t FUNCTOR_NOT = 1; // leaf
const uint8_t FUNCTOR_AND = 2; // internal quadtree_formula
const uint8_t FUNCTOR_OR = 3; // internal quadtree_formula
const uint8_t FUNCTOR_EXTEND = 4; // internal quadtree_formula

const uint8_t NOT_A_LEAF = 2;

const double NO_VALUE_LEAF = 3;
const double EMPTY_LEAF = 0;
const double FULL_LEAF = 1;
const double INTERNAL_NODE = 0.5;
// This indicates that one cannot determine the value of the quadtree_formula without computing the values of its children.
const double VALUE_NEED_CHILDREN = 2; // it indicates we need to compute the values of its children

struct quadtree_formula { // represents the output of the formula we are evaluating
    double val_leaf; // EMPTY_LEAF, FULL_LEAF, INTERNAL_NODE
    vector<quadtree_formula*> children;
    uint8_t k;
    uint8_t d; // number of attributes
};

// represents a subtree of the quadtree
struct subQuadtreeChild {
    ::qdag *qdag;
    uint16_t level; // the level of the quadtree_formula. 0 to start.
    uint64_t node; // absolute position in the bv[level] of the quadtree
    double value;// = NO_VALUE_LEAF;
    // can only be a number from the original quadtree (not extended)
};

// lqdag as a syntax tree
class lqdag {

public:
    typedef vector<uint64_t> att_set;

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
        this->subQuadtree->value = NO_VALUE_LEAF;
        this->lqdag1 = nullptr;
        this->lqdag2 = nullptr;
        for(uint64_t i = 0; i < subQuadtree->qdag->nAttr(); i++) // TODO: pass att set as param
            this->attribute_set_A.push_back(subQuadtree->qdag->getAttr(i));
    }

    // L = (QTREE, Q_r) o L = (NOT, Q_r)
    lqdag(uint8_t functor, subQuadtreeChild* subQuadtree) {
        assert(functor == FUNCTOR_QTREE || functor == FUNCTOR_NOT);
        this->functor = functor;
        this->subQuadtree = subQuadtree;
        this->subQuadtree->value = NO_VALUE_LEAF;
        this->lqdag1 = nullptr;
        this->lqdag2 = nullptr;
        for(uint64_t i = 0; i < subQuadtree->qdag->nAttr(); i++) // TODO: pass att set as param
            this->attribute_set_A.push_back(subQuadtree->qdag->getAttr(i));
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
    lqdag(uint8_t functor, lqdag* l, att_set &attribute_set_A, type_mapping_M *M = nullptr) {
        assert(functor == FUNCTOR_EXTEND);
        this->functor = functor;
        this->lqdag1 = l; // has
        this->attribute_set_A = attribute_set_A;
        this->subQuadtree = nullptr;
        this->lqdag2 = nullptr;

        if(M == nullptr) {
            uint16_t msize = std::pow(this->lqdag1->subQuadtree->qdag->getK(),
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
        if(this->functor == FUNCTOR_QTREE) {
            uint64_t nAtt = this->attribute_set_A.size();
            uint16_t cur_level = this->subQuadtree->level;
            uint64_t cur_node = this->subQuadtree->node;

            this->subQuadtree->value = is_leaf(cur_level, cur_node, this->subQuadtree->qdag->getK());
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
     * @return whether is an EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE
     */
    double is_leaf(uint16_t level, uint64_t node, uint8_t k){
        uint16_t max_level = this->subQuadtree->qdag->getHeight() - 1;
//        this->subQuadtree->qdag->get_ith_bit(level, node);
        if(level > max_level){ // level + 1 == max_level
            return this->subQuadtree->qdag->get_ith_bit(level-1, node);//(node_description & 255) == 255;
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
            case 8:
                val_full_ones = 255;
                break;
            case 16:
                val_full_ones = 65535;
            default:
                throw "error: invalid k_d";
        }
        if((node_description & val_full_ones) == val_full_ones){ // check the children
            if(level == max_level)
                return FULL_LEAF;
            uint64_t siblings = this->subQuadtree->qdag->rank(level,node);
            for(uint8_t i = 0; i < k_d; i++){
                if(is_leaf(level+1, siblings*k_d + i, k) != 1)
                    return INTERNAL_NODE;
            }
            return FULL_LEAF;
        }
        return INTERNAL_NODE;
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
                if(this->subQuadtree->value == NO_VALUE_LEAF)
                    this->subQuadtree->value = this->value_quadtree();
                return this->subQuadtree->value;
            }
            case FUNCTOR_NOT: {
                if(this->subQuadtree->value == NO_VALUE_LEAF)
                    this->subQuadtree->value= this->value_quadtree();
                return 1 - this->subQuadtree->value;
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
                uint16_t cur_node = this->subQuadtree->node;
                bool node_exists = true;
                uint64_t ith_child = this->subQuadtree->qdag->get_child(cur_level, cur_node, i, node_exists);
                cur_level++;
                double value = node_exists ? NO_VALUE_LEAF : EMPTY_LEAF;
                subQuadtreeChild *subQuadtree = new subQuadtreeChild{this->subQuadtree->qdag, cur_level, ith_child, value};
//                this->get_child_quadtree(i,subQuadtree);
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
                lqdag *l = new lqdag(FUNCTOR_EXTEND, this->lqdag1->get_child_lqdag(this->M[i]), this->attribute_set_A, this->M);
                return l;
            }
            default:
                throw "error: value_lqdag non valid functor";

        }
    }

    /**
     * Create an empty or full leaf.
     * @param val can be EMPTY_LEAF or FULL_LEAF.
     * @return quadtree_formula with the value of the leaf.
     */
    quadtree_formula* create_leaf(double val){
        quadtree_formula* newNode = new quadtree_formula();
        newNode->val_leaf = val;
        newNode->children = {};
        newNode->k = 0;
        newNode->d = 0;
        return newNode;
    }

    /**
     * Algorithm 9 of paper.
     * The completion Q_f is the quadtree representing the output of the formula F, represented as an lqdag.
     * @param k
     *
     */
    quadtree_formula* completion(uint64_t p,
                                 uint16_t max_level,
                                 uint16_t level,
                                 uint64_t UPPER_BOUND,
                                 uint64_t &results) {

        if(results >= UPPER_BOUND){ // top k results only
            quadtree_formula* newNode = create_leaf(EMPTY_LEAF);
            return newNode;
        }
        // if results < UPPER_BOUND
        double val_lqdag = this->value_lqdag();
        // return a leaf
        if(val_lqdag == EMPTY_LEAF || val_lqdag == FULL_LEAF){
            quadtree_formula* newNode = create_leaf(val_lqdag);
            if(val_lqdag == FULL_LEAF){
                if(level == max_level)
                    results++;
                else if(level < max_level)
                    results += p * (max_level - level);
            }

            return newNode;
        }

        double max_value = 0;
        double min_value = 1;
        vector<quadtree_formula*> Q_f_children;
        // if all quadtres are empty, return a 0 leaf.
        // C_i ← Completion(Child(L_F,i))
        for(uint64_t i = 0; i < p; i++){
            quadtree_formula* newNode = this->get_child_lqdag(i)->completion(p, max_level, level + 1, UPPER_BOUND, results);
            Q_f_children.push_back(newNode);
            max_value = max(max_value, newNode->val_leaf); // val_leaf can be EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE
            min_value = min(min_value, newNode->val_leaf);
        }
        // return a leaf
        if(max_value == EMPTY_LEAF || min_value == FULL_LEAF){ // return a leaf
            Q_f_children.clear(); // memory
            double val_leaf = (max_value == EMPTY_LEAF)? EMPTY_LEAF : FULL_LEAF;
            quadtree_formula* newNode = create_leaf(val_leaf);
            if(val_lqdag == FULL_LEAF){
                if(level == max_level)
                    results++;
                else if(level < max_level)
                    results += p * (max_level - level);
            }
            return newNode;
        }
        else {
            // return an internal node with children
            quadtree_formula *quadtree = new quadtree_formula();
            quadtree->val_leaf = val_lqdag; //INTERNAL_NODE;
            quadtree->children = Q_f_children;
            quadtree->k = 0;
            quadtree->d = 0;
            return quadtree;
        }
    }


    void print(int depth = 0) {
        switch (this->functor) {
            case FUNCTOR_QTREE: {
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                cout << " QTREE" << endl;
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                cout << " " << this->subQuadtree->level << " " << this->subQuadtree->node << endl;
                break;
            }
            case FUNCTOR_NOT: {
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                cout << " NOT" << endl;
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                this->lqdag1->print(depth + 1);
                break;
            }
            case FUNCTOR_AND: {
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                cout << " AND" << endl;
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                this->lqdag1->print(depth + 1);
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                this->lqdag2->print(depth + 1);
                break;
            }
            case FUNCTOR_OR: {
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                cout << " OR" << endl;
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                this->lqdag1->print(depth + 1);
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                this->lqdag2->print(depth + 1);
                break;
            }
            case FUNCTOR_EXTEND: {
                for (int i = 0; i < depth; ++i) {
                    std::cout << ">";
                }
                cout << " EXTEND";
                cout << " att_a: " << this->attribute_set_A << endl;
                this->lqdag1->print(depth + 1);
                break;
            }
            default:
                throw "error: value_lqdag non valid functor";

        }
    }
};

#endif