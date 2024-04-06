//
// We will not use qdags, only quadtrees
//
#include<bits/stdc++.h>
#include "./qdags.hpp"

const uint8_t FUNCTOR_QTREE = 0; // leaf
const uint8_t FUNCTOR_NOT = 1; // leaf
const uint8_t FUNCTOR_AND = 2; // internal node
const uint8_t FUNCTOR_OR = 3; // internal node
const uint8_t FUNCTOR_EXTEND = 4; // internal node
const double VALUE_NEED_CHILDREN = 2; // it indicates we need to compute the values of its children

struct node {
    double val_leaf; // 0, 1, 1/2
    vector<node*> children;
    //att_set attribute_set;
    uint16_t height;
    uint8_t k;
    uint8_t d;
};

// represents a subtree of the quadtree
struct subQuadtreeChild {
    qdag *qdag;
    int16_t level; // the level of the node. -1 for the root
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
        this->attribute_set_A = {};
    }

    // L = (QTREE, Q_r) o L = (NOT, Q_r)
    lqdag(uint8_t functor, subQuadtreeChild* subQuadtree) {
        assert(functor == FUNCTOR_QTREE || functor == FUNCTOR_NOT);
        this->functor = functor;
        this->subQuadtree = subQuadtree;
        this->lqdag1 = nullptr;
        this->lqdag2 = nullptr;
        this->attribute_set_A = {};
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
        this->attribute_set_A = attribute_set_A;

        this->subQuadtree = nullptr;
        this->lqdag2 = nullptr;
    }


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
     * Now leaves can be in the last lever or higher.
     * @return 1 if the qdag represents a full single cell, 0 if it is empty, 1/2 if is an internal node.
     */
    double value_quadtree(){
        if(this->functor == FUNCTOR_QTREE) {
            uint16_t max_level = this->subQuadtree->qdag->getHeight() - 1;
            uint16_t cur_level = this->subQuadtree->level;
            uint64_t cur_node = this->subQuadtree->node;
            // grid side is 1 or leaf at last level
            if(this->get_grid_side() == 1 || cur_level == max_level) // TODO: esta bien cambiar asi la grid side?, no es lo mismo ambas condiciones?
                return this->subQuadtree->qdag->get_ith_bit(cur_level, cur_node); // TODO: see if problems bcause get_ith_bit returns bool.
            // leaf higher level (subgrid full of 1s or 0s)
            uint8_t is_leaf = this->subQuadtree->qdag->get_node_last_level(cur_level, cur_node);
            if(is_leaf == 0 || (is_leaf & 255) == 1)
                return is_leaf;
            else
                return 0.5;
        } else{
            //cout << "error: value_quadtree called on non-quadtree" << endl;
            throw "error: value_quadtree called on non-quadtree";
        }
    }

    /**
     * Algorithm 2 child(Q,i)
     * @param i
     * @param subQuad a subQuadtreeChild which represents the i-th child of a node of the quadtree.
     * @return a subQuadtreeChild which represents the i-th child of a node of the quadtree.
     */
    void get_child_quadtree(uint64_t i, subQuadtreeChild* subQuad){
        // assert functor == QUADTREE
        if(this->functor != FUNCTOR_QTREE){
            throw "error: get_child_quadtree called on non-quadtree";
        }
        uint16_t curr_level = this->subQuadtree->level;
        uint64_t curr_node = this->subQuadtree->node;
        uint64_t siblings = this->subQuadtree->qdag->rank(curr_level, curr_node); // number of nodes of the left of the node in that level in the qdag
        uint16_t new_level = curr_level+1;
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
        if(this->functor == FUNCTOR_QTREE){
            return this->value_quadtree();
        }
        else if(this->functor == FUNCTOR_NOT){
            return 1 - this->value_quadtree();
        }
        else if(this->functor == FUNCTOR_AND){
            if(this->lqdag1->value_lqdag() == 0 || this->lqdag2->value_lqdag() == 0)
                return 0;
            if(this->lqdag1->value_lqdag() == 1)
                return this->lqdag2->value_lqdag();
            if(this->lqdag2->value_lqdag() == 1)
                return this->lqdag1->value_lqdag();
            return VALUE_NEED_CHILDREN;
        }
        else if(this->functor == FUNCTOR_OR){
            if(this->lqdag1->value_lqdag() == 1 || this->lqdag2->value_lqdag() == 1)
                return 1;
            if(this->lqdag1->value_lqdag() == 0)
                return this->lqdag2->value_lqdag();
            if(this->lqdag2->value_lqdag() == 0)
                return this->lqdag1->value_lqdag();
            return VALUE_NEED_CHILDREN;
        }
        else if(this->functor == FUNCTOR_EXTEND){
            return this->lqdag1->value_lqdag();
        }
        throw "error: value_lqdag non valid functor";
    }

    // algorithm 8 child(L,i)
    /**
     * Can only be invoked when value is 1/2 or 2.
     * @param i
     * @return
     */
    lqdag* get_child_lqdag(uint64_t i){
        // case base
        if(this->functor == FUNCTOR_QTREE || this->functor == FUNCTOR_NOT){
            subQuadtreeChild* subQuadtree = new subQuadtreeChild();
            this->get_child_quadtree(i, subQuadtree); // TODO: check if it's necessary to create another subquadtreechild
            lqdag *l = new lqdag(this->functor, subQuadtree);
            return l; // TODO: see memory leaks
        }
        else if(this->functor == FUNCTOR_AND){
            if(this->lqdag1->value_lqdag() == 1)
                return lqdag2->get_child_lqdag(i);
            if(this->lqdag2->value_lqdag() == 1)
                return lqdag1->get_child_lqdag(i);
            lqdag *l = new lqdag(FUNCTOR_AND, this->lqdag1->get_child_lqdag(i), this->lqdag2->get_child_lqdag(i));
            return l;
        }
        else if(this->functor == FUNCTOR_OR){
            if(this->lqdag1->value_lqdag() == 0)
                return lqdag2->get_child_lqdag(i);
            if(this->lqdag2->value_lqdag() == 0)
                return lqdag1->get_child_lqdag(i);
            lqdag *l = new lqdag(FUNCTOR_OR, this->lqdag1->get_child_lqdag(i), this->lqdag2->get_child_lqdag(i));
            return l;
        }
        else if(this->functor == FUNCTOR_EXTEND){ // Extend quadtree to attributes att_set
            // TODO: completar esto con el paper
            // TODO: pregunta no entiendo lo del paper
            uint16_t dim = this->attribute_set_A.size();
            uint16_t dim_prime = this->lqdag1->attribute_set_A.size();

            type_mapping_M M = 0;

            uint64_t mask;
            uint64_t i_prime;

            mask = 1 << (dim_prime - 1); // equivalent to 2^(dim_prime-1)
            i_prime = 0;

            for (uint16_t j = 0; j < dim_prime; ++j) {
                if (i & (1 << (dim - this->subQuadtree->qdag->getAttr(j) - 1)))
                    i_prime |= mask;
                mask >>= 1;
            }
            lqdag *l = new lqdag(FUNCTOR_EXTEND, this->lqdag1->get_child_lqdag(i_prime), this->attribute_set_A);
            return l;
        }

        throw "error: value_lqdag non valid functor";
    }


    /**
     * Algorithm 9 of paper.
     * The completion Q_f is the quadtree representing the output of the formula F, represented as an lqdag
     *
     */
    node* completion(){
        if(this->value_lqdag() == 0 || this->value_lqdag() == 1){
            node* newNode = new node();
            newNode->val_leaf = this->value_lqdag();
            newNode->children = {};
            newNode->height = 1;
            newNode->k = 0;
            newNode->d = 0;
            return newNode;
            /* if we return a qdag* :
            vector<uint64_t> bv[1];
            //bv.push_back(this->value_lqdag());
            if(this->value_lqdag()) // if it's 1, put a 1 as a leaf
                bv[0] = vector<uint64_t>(0);
            // TODO: see que pasa cuando value = 0, bv empty... memory leakS?
            qdag* Q_f = new qdag(bv, this->attribute_set_A, 1, this->subQuadtree->qdag->getK(), this->subQuadtree->qdag->getD());
            return Q_f;*/
        }
        double max_value = 0;
        double min_value = 255;
        vector<node*> Q_f_children;
        // if all quadtres are empty, return a 0 leaf.
        // TODO: this is wrong the D dimension
        for(uint64_t i = 0; i < this->subQuadtree->qdag->getD(); i++){
            node* newNode = this->get_child_lqdag(i)->completion();
            Q_f_children.push_back(newNode);
            max_value = max(max_value, newNode->val_leaf);
            min_value = min(min_value, newNode->val_leaf);
        }
        if(max_value == 0 || min_value == 1){
            Q_f_children.clear(); // memory
            node* newNode = new node();
            newNode->val_leaf = max_value == 0 ? 0 : 1;
            newNode->children = {};
            newNode->height = 1;
            newNode->k = 0;
            newNode->d = 0;
            return newNode;
        }

        node* quadtree = new node();
        quadtree->val_leaf = 1/2;
        quadtree->children = Q_f_children;
        quadtree->height = 2;
        quadtree->k = this->subQuadtree->qdag->getK();
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
};