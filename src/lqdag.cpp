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

class lqdag {

public:
    typedef vector<uint64_t> att_set;

private:
    uint8_t functor; // QTREE, NOT, AND, OR, EXTEND
    qdag *Q;
    lqdag *lqdag1;
    lqdag *lqdag2;
    att_set attribute_set_A; // for extend functor

public:
    lqdag() = default;

    // L = (QTREE, Q_r)
    lqdag(qdag* Q) {
        this->functor = FUNCTOR_QTREE;
        this->Q = Q;
    }

    // L = (QTREE, Q_r) o L = (NOT, Q_r)
    lqdag(uint8_t functor, qdag* Q) {
        // TODO: assert functor is QTREE, NOT,
        this->functor = functor;
        this->Q = Q;
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

    double value(){
        if(this->functor == FUNCTOR_QTREE){
            return this->Q->value();
        }
        else if(this->functor == FUNCTOR_NOT){
            return 1 - this->Q->value();
        }
        else if(this->functor == FUNCTOR_AND){
            if(this->lqdag1->value() == 0 || this->lqdag2->value() == 0){
                return 0;
            }
            else if(this->lqdag1->value() == 1){
                return this->lqdag2->value();
            }
            else if(this->lqdag2->value() == 1){
                return this->lqdag1->value();
            }
            else{
                return VALUE_NEED_CHILDREN;
            }
        }
        else if(this->functor == FUNCTOR_OR){
            if(this->lqdag1->value() == 1 || this->lqdag2->value() == 1){
                return 1;
            }
            else if(this->lqdag1->value() == 0){
                return this->lqdag2->value();
            }
            else if(this->lqdag2->value() == 0){
                return this->lqdag1->value();
            }
            else{
                return VALUE_NEED_CHILDREN;
            }
        }
        else if(this->functor == FUNCTOR_EXTEND){
            return this->lqdag1->value();
        }
    }

    // algorithm 8 child(L,i)
    // TODO: pregunta get_child_qdag obtengo toooodo el qdag? o el primer nivel?
    lqdag* get_sub_lqdag(uint16_t i){
        if(this->functor == FUNCTOR_QTREE || this->functor == FUNCTOR_NOT){
            lqdag *l = new lqdag();
            l->functor = this->functor;
            l->Q = this->Q->get_child_qdag(i);
            return l;
        }
        else if(this->functor == FUNCTOR_AND){
            if(this->lqdag1->value() == 1){
                return lqdag2->get_sub_lqdag(i);
            }
            else if(this->lqdag2->value() == 1){
                return lqdag1->get_sub_lqdag(i);
            }
            lqdag *l = new lqdag();
            l->functor = FUNCTOR_AND;
            l->lqdag1 = this->lqdag1->get_sub_lqdag(i);
            l->lqdag2 = this->lqdag2->get_sub_lqdag(i);
            return l;
        }
        else if(this->functor == FUNCTOR_OR){
            if(this->lqdag1->value() == 0){
                return lqdag2->get_sub_lqdag(i);
            }
            else if(this->lqdag2->value() == 0){
                return lqdag1->get_sub_lqdag(i);
            }
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
        else{
            // error functor
            cout << "functor unrecognized" << endl;
        }
    }


};