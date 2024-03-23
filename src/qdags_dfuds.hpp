//
// Created by Asunción Gómez on 14-03-24.
//

#ifndef QDAGS_QDAGS_DFUDS_HPP
#define QDAGS_QDAGS_DFUDS_HPP

#include <sdsl/bit_vectors.hpp>
#include<bits/stdc++.h>
#include<ratio>
#include<chrono>
#include<ctime>

#include "se_quadtree_dfuds.hpp"

using namespace std::chrono;

extern high_resolution_clock::time_point start_rank, stop_rank;
extern double total_time_rank;
extern duration<double> time_span_rank;

typedef uint8_t type_mapping_M;

bool compare_pairs_dfuds(const pair<uint64_t, uint64_t> &i, const pair<uint64_t, uint64_t> &j) {
    return i.second < j.second;
}

class qdag_dfuds {
public:

    typedef vector<uint64_t> att_set;

private:
    se_quadtree_dfuds *Q;

    type_mapping_M *M;    // mapping

    att_set attribute_set;

    uint64_t grid_side;

    uint16_t Msize;  // number of children of every qdag node, k^d

    bool is_extended_qdag;

    vector<vector<type_mapping_M> *> M_prime;

    int32_t tab_extend_5[16];   // queries of 5 attributes, i.e., dimension 2^5=32
    int32_t tab_extend_4[16];
    int32_t tab_extend_3[16];

public:

    qdag_dfuds() = default;

    uint64_t size() {
        uint64_t s = Q->size() + Msize * sizeof(uint16_t) + attribute_set.size() * sizeof(uint64_t)
                     + M_prime.size() * sizeof(vector<type_mapping_M> *) + sizeof(uint64_t);
        return s;
    };

    void setAtts(uint64_t att1, uint64_t att2) {
        attribute_set[0] = att1;
        attribute_set[1] = att2;
    }

    qdag_dfuds(const qdag_dfuds &_Q){
        this->Q = _Q.Q;
        _Q.Q->inc_ref_count();
        this->M = _Q.M;
        for (uint64_t i = 0; i < _Q.attribute_set.size(); i++)
            this->attribute_set.push_back(_Q.attribute_set[i]);

        this->grid_side = _Q.grid_side;
        this->Msize = _Q.Msize;
        this->is_extended_qdag = _Q.is_extended_qdag;
    }

    qdag_dfuds(std::vector<std::vector<uint64_t>> &points,
        att_set &_attribute_set,
        const uint64_t _grid_side,
        uint8_t k, uint8_t d){

        Msize = std::pow(k, d);

        M = new type_mapping_M[Msize];

        uint64_t i, j;

        for (i = 0; i < Msize; i++)
            M[i] = i;  // identity mapping

        attribute_set = _attribute_set;

        vector<uint64_t> tuple_aux(d);
        vector<pair<uint64_t, uint64_t>> map_sort_att(d);

        for (i = 0; i < d; i++)
            map_sort_att[i] = make_pair(i, attribute_set[i]);

        std::sort(map_sort_att.begin(), map_sort_att.end(), compare_pairs_dfuds);

        for (i = 0; i < points.size(); i++) {
            //if (i%1000000==0) cout << i << endl;
            for (j = 0; j < d; j++)
                tuple_aux[j] = points[i][map_sort_att[j].first];

            for (j = 0; j < d; j++)
                points[i][j] = tuple_aux[j];

        }

        std::sort(attribute_set.begin(), attribute_set.end());

        Q = new se_quadtree_dfuds(points, _grid_side, k, d);

        grid_side = _grid_side;
        is_extended_qdag = false;
    }

    ~qdag_dfuds(){
        if (Q && !is_extended_qdag) {
            delete Q;
        }
        if (is_extended_qdag) delete M;
    }

    qdag_dfuds *extend(att_set &attribute_set_A){
        uint16_t dim = attribute_set_A.size(); // d
        uint16_t dim_prime = attribute_set.size(); // d'
        uint64_t p = std::pow(Q->getK(), dim);       // tamaño del mapeo, número de hijos por nodo en el quadtree original

        type_mapping_M *_M = new type_mapping_M[p]; // mapeo

        uint64_t mask;
        uint64_t i, i_prime;

        // TODO: see if we have to change something (order DFS, idk)
        for (i = 0; i < p; ++i) {
            // todos los bits están en cero excepto el bit en la posición dim_prime - 1.
            mask = 1 << (dim_prime - 1); // equivalent to 2^(dim_prime-1)
            i_prime = 0;

            for (uint16_t j = 0; j < dim_prime; ++j) {
                if (i & (1 << (dim - attribute_set[j] - 1)))
                    i_prime |= mask;

                mask >>= 1;
            }

            _M[i] = M[i_prime];
        }

        qdag_dfuds *q = new qdag_dfuds();

        q->Q = this->Q;
        q->M = _M;
        q->attribute_set = attribute_set_A;
        q->grid_side = this->grid_side;
        q->is_extended_qdag = true;
        q->Msize = p; // this.Msize;

        return q; //el nuevo qdag extendido
    }

    uint64_t nAttr() {
        return attribute_set.size();
    }


    uint64_t getAttr(uint64_t i) {
        return attribute_set[i];
    }


    uint64_t getGridSide() {
        return grid_side;
    }


    uint64_t getHeight() {
        return Q->getHeight();
    }


    uint8_t getK() {
        return Q->getK();
    }

    uint8_t getD() {
        return Q->getD();
    }


    uint16_t nChildren() {
        return Msize;
    }

    uint64_t getKD() {
        return Q->getKD();
    }


    /** the mapping between the children
     * for example:
     * original quadtree: 1011
     * qdag: 1011 1011
     * getM(5) = 1
     **/
    uint16_t getM(uint16_t i) {
        return M[i];
    }

    // This is for a binary relation, i.e., a k^2-tree with 4 children per node
    // son tablas precomputadas para materializar los nodos
    void createTableExtend5()  // para join con 5 atributos, se constuyen al momento del join
    {
        uint64_t i, j;
        uint32_t x, B;

        if (Q->getKD() == 2)
            B = 4;
        else
            B = 16;

        for (i = 0; i < B; i++) {
            x = 0;
            for (j = 0; j < 32; j++)
                if (i & (1 << M[j]))
                    x = (x << 1) | 1;
                else
                    x = (x << 1);

            tab_extend_5[i] = x;
            //cout << x << endl;
        }
    }

    // This is for a binary relation, i.e., a k^2-tree with 4 children per node
    void createTableExtend4() // para join con 4 atributos
    {
        uint64_t i, j;
        uint32_t x, B;

        if (Q->getKD() == 2)
            B = 4;
        else
            B = 16;

        for (i = 0; i < B; i++) {
            x = 0;
            for (j = 0; j < 16; j++)
                if (i & (1 << M[j]))
                    x = (x << 1) | 1;
                else
                    x = (x << 1);

            tab_extend_4[i] = x << 16;
            //cout << std::hex << x << endl;
        }
    }

    // This is for a binary relation, i.e., a k^2-tree with 4 children per node
    void createTableExtend3() // para join con 3 atributos
    {
        uint64_t i, j;
        uint32_t x, B = 16;

        //if (Q->getKD() == 2)
        //    B = 4;
        //else
        //    B = 16;

        for (i = 0; i < B; i++) {
            x = 0;
            for (j = 0; j < 8; j++)
                if (i & (1 << M[j]))
                    x = (x << 1) | 1;
                else
                    x = (x << 1);

            tab_extend_3[i] = x << 24;
            //cout << std::hex << x << endl;
        }
    }


    // TODO: see how the extend works in DFS!
    // we need the number of 1s of the level (until the node position)
    inline uint32_t materialize_node_3(uint64_t node, uint64_t *rank_vector) {
        // get the number of 1s of that level
        return tab_extend_3[Q->get_node(node, rank_vector)];
    }


    inline uint32_t materialize_node_4(uint64_t node, uint64_t *rank_vector) {
        return tab_extend_4[Q->get_node(node, rank_vector)];
    }


    inline uint32_t materialize_node_5(uint64_t node, uint64_t *rank_vector) {
        return tab_extend_5[Q->get_node(node, rank_vector)];
    }


    inline uint32_t materialize_node_3_lastlevel(uint64_t node) {
        return tab_extend_3[Q->get_node_last_level(node)];
    }


    inline uint32_t materialize_node_4_lastlevel(uint64_t node) {
        return tab_extend_4[Q->get_node_last_level(node)];
    }


    inline uint32_t materialize_node_5_lastlevel(uint64_t node) {
        return tab_extend_5[Q->get_node_last_level(node)];
    }

    void printBv() {
        //Q->printBv();
    }

    // TODO: change it by leaf num!!! it is not the same
    uint64_t get_num_leaves(uint64_t node) {
//        return Q->subtree(node);
        return Q->leaf_num(node);
    }

    // TODO: only for testing
    se_quadtree_dfuds *getQ() {
        return Q;
    }

};



#endif //QDAGS_QDAGS_DFUDS_HPP
