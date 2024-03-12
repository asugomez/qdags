//
// We will not use qdags, only quadtrees
//
#include<bits/stdc++.h>
#include "./qdags.hpp"

const uint8_t FUNCTOR_QTREE = 0;
const uint8_t FUNCTOR_NOT = 1;
const uint8_t FUNCTOR_AND = 2;
const uint8_t FUNCTOR_OR = 3;
const uint8_t FUNCTOR_EXTEND = 4;
const double VALUE_NEED_CHILDREN = 2;


// represents a subtree of the quadtree
struct subQuadtreeChild {
    qdag *Qdag;
    uint16_t level; // the level of the node. -1 for the root
    uint64_t node;
};

class lqdag {

public:
    typedef vector<uint64_t> att_set;

private:
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
    }

    // L = (QTREE, Q_r) o L = (NOT, Q_r)
    lqdag(uint8_t functor, subQuadtreeChild* subQuadtree) {
        // TODO: assert functor is QTREE, NOT,
        this->functor = functor;
        this->subQuadtree = subQuadtree;
    }

    // L = (AND, L1, L2) o L = (OR, L1, L2)
    lqdag(uint8_t functor, lqdag *l1, lqdag *l2) {
        // TODO: assert functor is AND, OR,
        this->functor = functor;
        this->lqdag1 = l1;
        this->lqdag2 = l2;
    }

    // L = (EXTEND, L1, A)
    lqdag(uint8_t functor, lqdag* l, att_set &attribute_set_A) {
        // TODO: assert functor is EXTEND
        this->functor = functor;
        this->lqdag1 = l;
        this->attribute_set_A = attribute_set_A;
    }

    /**
     * Algorithm 1 Value(Q)
     * Value of a qdag
     * @return 1 if the grid is a single point, 0 if it is empty, 1/2 otherwise.
     */
    double value_quadtree(){
        if(this->functor == FUNCTOR_QTREE) {
            // see the level and the node
            uint16_t height = this->subQuadtree->Qdag->getHeight();
            uint16_t curr_level = this->subQuadtree->level;
            uint64_t curr_node = this->subQuadtree->node;
            // if the node is a leaf (grid side = 1)
            // TODO: pregunta difference between l= 1 and a leaf
            // grid side = 1
            if(curr_level == height - 1)
                return 0;
            // leaf
            if(curr_level == height) // TODO: MALO ESTO!!
                return this->subQuadtree->Qdag->get_ith_bit(curr_level, curr_node); // TODO: see if problems bcause get_ith_bit returns bool.
            else
                return 0.5;
        } else{
            cout << "error: value_quadtree called on non-quadtree" << endl;
            return 2;
        }
    }

    /**
     * Algorithm 2 child(Q,i)
     * @param i
     * @return a subQuadtreeChild which represents the i-th child of a node of the quadtree.
     */
    void get_child_quadtree(uint64_t i, subQuadtreeChild* subQuad){
        // assert functor == QUADTREE
        if(this->functor != FUNCTOR_QTREE){
            throw "error: get_child_quadtree called on non-quadtree";
        }
        uint16_t curr_level = this->subQuadtree->level;
        uint64_t curr_node = this->subQuadtree->node;
        uint64_t siblings = this->subQuadtree->Qdag->rank(curr_level,curr_node); // number of nodes of the left of the node in that level in the qdag
        uint16_t new_level = curr_level+1;
        uint64_t new_node = siblings * this->subQuadtree->Qdag->getKD() + i;
        subQuad->Qdag = this->subQuadtree->Qdag;
        subQuad->level = new_level;
        subQuad->node = new_node;
        //subQuad = subQuadtreeChild{this->subQuadtree->Qdag, new_level, new_node};
        //return subQuad; // TODO: test que entrega bien el puntero
    }

    /**
     * algorithm 7 value(L)
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
            else if(this->lqdag1->value_lqdag() == 1)
                return this->lqdag2->value_lqdag();
            else if(this->lqdag2->value_lqdag() == 1)
                return this->lqdag1->value_lqdag();
            return VALUE_NEED_CHILDREN;
        }
        else if(this->functor == FUNCTOR_OR){
            if(this->lqdag1->value_lqdag() == 1 || this->lqdag2->value_lqdag() == 1)
                return 1;
            else if(this->lqdag1->value_lqdag() == 0)
                return this->lqdag2->value_lqdag();
            else if(this->lqdag2->value_lqdag() == 0)
                return this->lqdag1->value_lqdag();
            return VALUE_NEED_CHILDREN;
        }
        else if(this->functor == FUNCTOR_EXTEND){
            return this->lqdag1->value_lqdag();
        }
        throw "error: value_lqdag non valid functor";
    }

    // algorithm 8 child(L,i)
    lqdag* get_sub_lqdag(uint64_t i){
        if(this->functor == FUNCTOR_QTREE || this->functor == FUNCTOR_NOT){
            lqdag *l = new lqdag();
            l->functor = this->functor;
            subQuadtreeChild* subQuadtree = new subQuadtreeChild();
            this->get_child_quadtree(i, subQuadtree);
            l->subQuadtree = subQuadtree;
            return l; // TODO: see memory leaks
        }
        else if(this->functor == FUNCTOR_AND){
            if(this->lqdag1->value_lqdag() == 1)
                return lqdag2->get_sub_lqdag(i);
            else if(this->lqdag2->value_lqdag() == 1)
                return lqdag1->get_sub_lqdag(i);
            lqdag *l = new lqdag();
            l->functor = FUNCTOR_AND;
            l->lqdag1 = this->lqdag1->get_sub_lqdag(i);
            l->lqdag2 = this->lqdag2->get_sub_lqdag(i);
            return l;
        }
        else if(this->functor == FUNCTOR_OR){
            if(this->lqdag1->value_lqdag() == 0)
                return lqdag2->get_sub_lqdag(i);
            else if(this->lqdag2->value_lqdag() == 0)
                return lqdag1->get_sub_lqdag(i);
            lqdag *l = new lqdag();
            l->functor = FUNCTOR_OR;
            l->lqdag1 = this->lqdag1->get_sub_lqdag(i);
            l->lqdag2 = this->lqdag2->get_sub_lqdag(i);
            return l;
        }
        else if(this->functor == FUNCTOR_EXTEND){
            // TODO: completar esto con el paper
            // TODO: pregunta no entiendo lo del paper
            uint8_t d = this->attribute_set_A.size();
            uint8_t d_prime = this->lqdag1->attribute_set_A.size();
            // first d bits of i
            uint64_t m_d = i >> (64-d); // TODO: ver si esto esta bien
            // the projection of m_d to the positions in which the attributes of A′ appear in A
            uint64_t m_d_prime = d_prime & m_d; // TODO: aqui puse cualqueir cosa
            uint16_t i_prime = 0; // TODO: esto no esta bien
            lqdag *l = new lqdag();
            l->functor = FUNCTOR_EXTEND;
            l->lqdag1 = this->lqdag1->get_sub_lqdag(i_prime);
            l->attribute_set_A = this->attribute_set_A;
            return l;
        }

        throw "error: value_lqdag non valid functor";
    }


};