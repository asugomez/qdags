//
// We will not use qdags, only quadtrees
//
#include<bits/stdc++.h>
#include "./qdags.hpp"

const uint8_t FUNCTOR_QTREE = 0; // leaf
const uint8_t FUNCTOR_NOT = 1; // leaf
const uint8_t FUNCTOR_AND = 2; // internal quadtree_formula
const uint8_t FUNCTOR_OR = 3; // internal quadtree_formula
const uint8_t FUNCTOR_EXTEND = 4; // internal quadtree_formula
const uint8_t NOT_A_LEAF = 2;

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
    qdag *qdag;
    uint16_t level; // the level of the quadtree_formula. 0 to start.
    uint64_t node;
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
    att_set attribute_set_A; // for extend functor

public:
    lqdag() = default;

    // L = (QTREE, Q_r)
    lqdag(subQuadtreeChild* subQuadtree) {
        this->functor = FUNCTOR_QTREE;
        this->subQuadtree = subQuadtree;
        this->lqdag1 = nullptr;
        this->lqdag2 = nullptr;
        for(uint64_t i = 0; i < subQuadtree->qdag->nAttr(); i++)
            this->attribute_set_A.push_back(subQuadtree->qdag->getAttr(i));
    }

    // L = (QTREE, Q_r) o L = (NOT, Q_r)
    lqdag(uint8_t functor, subQuadtreeChild* subQuadtree) {
        assert(functor == FUNCTOR_QTREE || functor == FUNCTOR_NOT);
        this->functor = functor;
        this->subQuadtree = subQuadtree;
        this->lqdag1 = nullptr;
        this->lqdag2 = nullptr;
        for(uint64_t i = 0; i < subQuadtree->qdag->nAttr(); i++)
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
    lqdag(uint8_t functor, lqdag* l, att_set &attribute_set_A) {
        assert(functor == FUNCTOR_EXTEND);
        this->functor = functor;
        this->lqdag1 = l;
        this->lqdag1->attribute_set_A = attribute_set_A;
        this->attribute_set_A = attribute_set_A;

        this->subQuadtree = nullptr;
        this->lqdag2 = nullptr;
    }

    // TODO: see if we are going to need this or not.
    uint64_t get_grid_side(){
        if(this->functor == FUNCTOR_QTREE) {
            uint16_t curr_level = this->subQuadtree->level;
            uint64_t k= this->subQuadtree->qdag->getK();
            // TODO: debug see if ceil is working!
            return ceil(this->subQuadtree->qdag->getGridSide() / (pow(k, curr_level)));
        } else{
            cout << "error: get_grid_side called on non-quadtree" << endl;
            return 0;
        }
    }

    /**
     * Algorithm 1 Value(Q), adapted for lqdag
     * Value of a qdag
     * In lqdags we introduce a new idea: full leaves, that denote subgrids full of 1s.
     * Leaves can be in a higher level than the last one, when the subgrids are all 0s or all 1s.
     * @return 1 if the qdag represents a full single cell, 0 if it is empty, 1/2 if is an internal quadtree_formula.
     */
    double value_quadtree(){
        if(this->functor == FUNCTOR_QTREE) {
            uint16_t max_level = this->subQuadtree->qdag->getHeight() - 1;
            uint16_t cur_level = this->subQuadtree->level;
            uint64_t cur_node = this->subQuadtree->node;
            // grid side is 1 or leaf at last level
            // this->get_grid_side() == 1 ?? // TODO: esta bien cambiar asi la grid side?, no es lo mismo ambas condiciones?
            if(cur_level == max_level)
                return this->subQuadtree->qdag->get_ith_bit(cur_level, cur_node); // TODO: see if problems bcause get_ith_bit returns bool.
            // leaf higher level (subgrid full of 1s or 0s)
            // TODO: see when it's a full leaf (all points are 1, not only subgrids)
            uint8_t val_leaf = is_leaf(cur_level, cur_node, this->subQuadtree->qdag->getK(), this->subQuadtree->qdag->getKD());

            if(val_leaf == 0 || val_leaf == 1)
                return val_leaf;
            else
                return 0.5;
        } else{
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
     * @return 0 if is a leaf full of 0s, 1 if is a leaf full of 1s and NOT_A_LEAF otherwise.
     */
    uint8_t is_leaf(uint16_t level, uint64_t node, uint8_t k, uint64_t k_d){
        uint8_t node_description = this->subQuadtree->qdag->get_node_last_level(level, node);
        uint16_t max_level = this->subQuadtree->qdag->getHeight() - 1;
        if(level == max_level){
            return (node_description & 255) == 255;
        }
        if((node_description & 255) == 0)
            return 0;
        if((node_description & 255) == 255){ // check the children
            for(uint8_t i = 0; i < k; i++){
                if(is_leaf(level+1, node*k_d + i, k, k_d) != 1)
                    return NOT_A_LEAF;
            }
            return 1;
        }
        return NOT_A_LEAF;
    }

    /**
     * Algorithm 2 child(Q,i)
     * @param i
     * @param subQuad a subQuadtreeChild which represents the i-th child of a quadtree_formula of the quadtree.
     * @return a subQuadtreeChild which represents the i-th child of a quadtree_formula of the quadtree.
     */
    void get_child_quadtree(uint64_t i, subQuadtreeChild* subQuad){
        // assert functor == QUADTREE
        if(this->functor != FUNCTOR_QTREE){
            throw "error: get_child_quadtree called on non-quadtree";
        }
        uint16_t cur_level = this->subQuadtree->level;
        uint64_t cur_node = this->subQuadtree->node;
        uint64_t siblings = this->subQuadtree->qdag->rank(cur_level, cur_node); // number of nodes of the left of the quadtree_formula in that level in the qdag
        uint16_t new_level = cur_level + 1;
        uint64_t new_node = siblings * this->subQuadtree->qdag->getKD() + i;
        subQuad->qdag = this->subQuadtree->qdag;
        subQuad->level = new_level;
        subQuad->node = new_node;
        //subQuad = subQuadtreeChild{this->subQuadtree->qdag, new_level, new_node};
        //return subQuad; // TODO: test que entrega bien el puntero
    }

    /**
     * algorithm 7 value(L)
     * In lqdags we introduce a new idea: full leaves, that denote subgrids full of 1s
     * @return the value of the root of the lqdag.
     */
    double value_lqdag(){
        switch(this->functor){
            case FUNCTOR_QTREE:
                return this->value_quadtree();
            case FUNCTOR_NOT:
                return 1 - this->value_quadtree();
            case FUNCTOR_AND: {
                double val_lqdag1 = this->lqdag1->value_lqdag();
                double val_lqdag2 = this->lqdag2->value_lqdag();
                if (val_lqdag1 == 0 || val_lqdag2 == 0)
                    return EMPTY_LEAF;
                if (val_lqdag1 == 1)
                    return val_lqdag2;
                if (val_lqdag2 == 1)
                    return val_lqdag1;
                return VALUE_NEED_CHILDREN; // if both roots have Value 1⁄2, one cannot be sure of the Value of the resulting root until the AND between the children of Q1 and Q2 has been computed
            }
            case FUNCTOR_OR: {
                double val_lqdag1 = this->lqdag1->value_lqdag();
                double val_lqdag2 = this->lqdag2->value_lqdag();
                if (val_lqdag1 == 1 || val_lqdag2 == 1)
                    return FULL_LEAF;
                if (val_lqdag1 == 0)
                    return val_lqdag2;
                if (val_lqdag2 == 0)
                    return val_lqdag1;
                return VALUE_NEED_CHILDREN; // if both roots have Value 1⁄2, one cannot be sure of the Value of the resulting root until the OR between the children of Q1 and Q2 has been computed
            }
            case FUNCTOR_EXTEND:
                return this->lqdag1->value_lqdag();
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
                subQuadtreeChild *subQuadtree = new subQuadtreeChild();
                this->get_child_quadtree(i,subQuadtree);
                lqdag *l = new lqdag(this->functor, subQuadtree);
                return l;
            }
            case FUNCTOR_AND: {
                if (this->lqdag1->value_lqdag() == 1)
                    return lqdag2->get_child_lqdag(i);
                if (this->lqdag2->value_lqdag() == 1)
                    return lqdag1->get_child_lqdag(i);
                lqdag *l = new lqdag(FUNCTOR_AND, this->lqdag1->get_child_lqdag(i), this->lqdag2->get_child_lqdag(i));
                return l;
            }
            case FUNCTOR_OR: {
                if (this->lqdag1->value_lqdag() == 0)
                    return lqdag2->get_child_lqdag(i);
                if (this->lqdag2->value_lqdag() == 0)
                    return lqdag1->get_child_lqdag(i);
                lqdag *l = new lqdag(FUNCTOR_OR, this->lqdag1->get_child_lqdag(i), this->lqdag2->get_child_lqdag(i));
                return l;
            }
            case FUNCTOR_EXTEND: { // Extend quadtree to attributes att_set
                // same mapping as in the extend function but only for this i. See extend in qdags.hpp
                uint16_t dim = this->attribute_set_A.size();
                uint16_t dim_prime = this->lqdag1->attribute_set_A.size();

                type_mapping_M M = 0;

                uint64_t mask;
                uint64_t i_prime;

                mask = 1 << (dim_prime - 1); // equivalent to 2^(dim_prime-1)
                i_prime = 0;

                for (uint16_t j = 0; j < dim_prime; ++j) {
                    if (i & (1 << (dim - this->attribute_set_A[j]- 1)))
                        i_prime |= mask;
                    mask >>= 1;
                }
                lqdag *l = new lqdag(FUNCTOR_EXTEND, this->lqdag1->get_child_lqdag(i_prime), this->attribute_set_A);
                return l;
            }
            default:
                throw "error: value_lqdag non valid functor";

        }
    }


    /**
     * Algorithm 9 of paper.
     * The completion Q_f is the quadtree representing the output of the formula F, represented as an lqdag
     *
     */
    quadtree_formula* completion(uint8_t dim){ // TODO: see how to get dimension
        double val_lqdag = this->value_lqdag();
        // return a leaf
        if(val_lqdag == EMPTY_LEAF || val_lqdag == FULL_LEAF){
            quadtree_formula* newNode = new quadtree_formula();
            newNode->val_leaf = val_lqdag;
            newNode->children = {};
            // TODO: check if we need to write k and d.
            newNode->k = 0;
            newNode->d = 0;
            return newNode;
        }
        double max_value = 0;   // 00000000
        double min_value = 255; // 11111111
        vector<quadtree_formula*> Q_f_children;
        // if all quadtres are empty, return a 0 leaf.
        // C_i ← Completion(Child(L_F,i))
        uint64_t p = std::pow(2, dim); // TODO: check if base is always 2
        for(uint64_t i = 0; i < p; i++){
            quadtree_formula* newNode = this->get_child_lqdag(i)->completion(dim);
            Q_f_children.push_back(newNode);
            max_value = max(max_value, newNode->val_leaf); // val_leaf can be EMPTY_LEAF, FULL_LEAF or INTERNAL_NODE
            min_value = min(min_value, newNode->val_leaf);
        }
        // return a leaf
        if(max_value == 0 || min_value == 1){ // return a leaf
            Q_f_children.clear(); // memory
            quadtree_formula* newNode = new quadtree_formula();
            newNode->val_leaf = (max_value == EMPTY_LEAF)? EMPTY_LEAF : FULL_LEAF;
            newNode->children = {};
            newNode->k = 0;
            newNode->d = 0;
            return newNode;
        }
        // return an internal node with children
        quadtree_formula* quadtree = new quadtree_formula();
        quadtree->val_leaf = INTERNAL_NODE;
        quadtree->children = Q_f_children;
        quadtree->k = 0;
        quadtree->d = 0;
        return quadtree;
    }


    // ---------- FULL RELATIONAL ALGEBRA ---------- //
    // join = (and(and(extend(Q_TREE, (A,B,C)), (extend(Q_TREE, (A,B,C)))), extend(Q_TREE, (A,B,C)))
    //

    void pred(){
        // TODO: no tengo idea como hacerla
    }

    void selection(){

    }

    //(JOIN, L1(A1), L2(A2)) = (AND, (EXTEND, L1, A1 ∪ A2), (EXTEND, L2, A1 ∪ A2))
    void join(){

    }

    // (DIFF, L1(A), L2(A)) = (AND, L1, (NOT, L2))
    void diff(){

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