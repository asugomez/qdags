#ifndef INCLUDED_QDAGS
#define INCLUDED_QDAGS
// por la linea 194 esta el EXTEND del mapeo
#include <sdsl/bit_vectors.hpp>
#include "se_quadtree.hpp"

#include<bits/stdc++.h>
#include<ratio>
#include<chrono>
#include<ctime>

using namespace std::chrono;

extern high_resolution_clock::time_point start_rank, stop_rank;
extern double total_time_rank;
extern duration<double> time_span_rank;

typedef uint8_t type_mapping_M;

bool compare_pairs(const pair<uint64_t, uint64_t> &i, const pair<uint64_t, uint64_t> &j) {
    return i.second < j.second;
}


class qdag {
public:

    typedef vector<uint64_t> att_set;

private:
    se_quadtree *Q;

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

    qdag() = default;

    uint64_t size() {
        uint64_t s = Q->size() + Msize * sizeof(uint16_t) + attribute_set.size() * sizeof(uint64_t)
                     + M_prime.size() * sizeof(vector<type_mapping_M> *) + sizeof(uint64_t);

        //for (uint64_t i = 0; i < M_prime.size(); i++)
        //    s += M_prime[i]->size()*sizeof(type_mapping_M);

        return s;
    }


    void setAtts(uint64_t att1, uint64_t att2) {
        attribute_set[0] = att1;
        attribute_set[1] = att2;
    }

    qdag(const qdag &_Q) {
        this->Q = _Q.Q;
        _Q.Q->inc_ref_count();
        this->M = _Q.M;
        for (uint64_t i = 0; i < _Q.attribute_set.size(); i++)
            this->attribute_set.push_back(_Q.attribute_set[i]);

        this->grid_side = _Q.grid_side;
        this->Msize = _Q.Msize;
        this->is_extended_qdag = _Q.is_extended_qdag;

    }

    qdag(std::vector<std::vector<uint64_t>> &points,
         att_set &_attribute_set,
         const uint64_t _grid_side,
         uint8_t k, uint8_t d
    ) {

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

        std::sort(map_sort_att.begin(), map_sort_att.end(), compare_pairs);

        for (i = 0; i < points.size(); i++) {
            //if (i%1000000==0) cout << i << endl;
            for (j = 0; j < d; j++)
                tuple_aux[j] = points[i][map_sort_att[j].first];

            for (j = 0; j < d; j++)
                points[i][j] = tuple_aux[j];

        }

        std::sort(attribute_set.begin(), attribute_set.end());

        //cout << "Construyendo el quadtree" << endl;
        Q = new se_quadtree(points, _grid_side, k, d);

        grid_side = _grid_side;
        is_extended_qdag = false;

        //M_prime.reserve(Msize);

        //for (uint64_t i = 0; i < Msize; i++)
        //    M_prime.push_back(new std::vector<type_mapping_M>());

        //for (uint64_t i = 0; i < Msize; i++)
        //    M_prime[i]->push_back(i);

    }


    qdag(vector<uint64_t> bv[],
         att_set &_attribute_set,
         const uint64_t _grid_side,
         uint8_t k, uint8_t d
    ) {
        // OJO PIOJO: aqui el constructor no está terminado: no pasa el path de enteros a bits!
        Q = new se_quadtree(bv, _grid_side, k, d);

        Msize = std::pow(k, d);

        M = new type_mapping_M[Msize];

        for (uint64_t i = 0; i < Msize; i++)
            M[i] = i;  // identity mapping

        attribute_set = _attribute_set;
        //std::sort(attribute_set.begin(), attribute_set.end());
        grid_side = _grid_side;
        is_extended_qdag = false;

        //M_prime.reserve(Msize);

        //for (uint64_t i = 0; i < Msize; i++)
        //    M_prime.push_back(new std::vector<type_mapping_M>());

        //for (uint64_t i = 0; i < Msize; i++) {
        //    M_prime[M[i]]->push_back(i);
        //}

    }


    /*qdag(qdag &q, att_set &_attribute_set)
    {
        this->Q = q.Q;
    Msize = q.Msize;
    M = new type_mapping_M[Msize];
    for (uint64_t i = 0; i < Msize; i++)
        M[i] = q.M[i];

        attribute_set = _attribute_set;
    std::sort(attribute_set.begin(), attribute_set.end());
}*/

    ~qdag() {
        //if (Q && !is_extended_qdag) {
        //    delete Q;
        //    Q = NULL;
        //}
        if (is_extended_qdag) delete M;
    }

    // extend del paper!!
    qdag *extend(att_set &attribute_set_A) //extiende el qdag a ese set de atributos
    {
        uint16_t dim = attribute_set_A.size(); // d
        uint16_t dim_prime = attribute_set.size(); // d'
        uint64_t p = std::pow(Q->getK(), dim);       // tamaño del mapeo, número de hijos por nodo en el quadtree original

        type_mapping_M *_M = new type_mapping_M[p]; // mapeo

        uint64_t mask;
        uint64_t i, i_prime;

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

        qdag *q = new qdag();

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


    uint16_t getM(uint16_t i) {
        return M[i];
    }

    rank_bv_64 *getBv() {
        return Q->getBv();
    }

    inline uint64_t rank(uint16_t level, uint64_t node){
        return Q->rank(level,node);
    }

    void get_children(uint16_t level, uint64_t node, uint64_t *children_array, uint64_t &n_children) {
        return Q->get_children(level, node, children_array, n_children);
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

    // TODO: problema? siempre materializo,,, cuando uso el mapeo?
    // TODO: y porqué de 32 bits?
    // materializa el nodo del qdag, uno llega a un nodo en el quadtree que existe, pero conceptualmente
    // corresponde a un qdags q no existe en realidad.
    // esta funcion materializa el nodo en el qdag (donde no existe)
    // por ejem, quadtree con 4 nodos
    // pero trabajo en dimension 5 --> se extiende a 32 hijos
    // a partir de esos 4 bits, te genera el de 32 inmediatamente
    inline uint32_t materialize_node_3(uint64_t level, uint64_t node, uint64_t *rank_vector) {
        // DEBUG: ver roots[i] = node
        uint64_t r = Q->rank(level, node);
        return tab_extend_3[Q->get_node(level, node, rank_vector, r)];
    }


    inline uint32_t materialize_node_4(uint64_t level, uint64_t node, uint64_t *rank_vector) {
        uint64_t r = Q->rank(level, node);
        return tab_extend_4[Q->get_node(level, node, rank_vector, r)];
    }


    inline uint32_t materialize_node_5(uint64_t level, uint64_t node, uint64_t *rank_vector) {
        uint64_t r = Q->rank(level, node);
        return tab_extend_5[Q->get_node(level, node, rank_vector, r)];
    }


    inline uint32_t materialize_node_3_lastlevel(uint64_t level, uint64_t node) {
        return tab_extend_3[Q->get_node_last_level(level, node)];
    }


    inline uint32_t materialize_node_4_lastlevel(uint64_t level, uint64_t node) {
        return tab_extend_4[Q->get_node_last_level(level, node)];
    }


    inline uint32_t materialize_node_5_lastlevel(uint64_t level, uint64_t node) {
        return tab_extend_5[Q->get_node_last_level(level, node)];
    }


    //void print(std::ofstream &ofs)
    //{
    //    Q->print(ofs);
    //}
    void printBv() {
        Q->printBv();
    }

    /**
     * @param level of the node. -1 if is the root
     * @param node the i-th node (0 or 1) of the level
     * @return number of leaves of the ith-node of the level.
     * @example:
     * 0101
     * 0100 1001
     * 1111 0001 1000
     * get_num_leaves(-1,0) --> 6
     * get_num_leaves(0,0) --> 0
     * get_num_leaves(0,1) --> 4
     * get_num_leaves(2,7) --> 1
     */
    uint64_t get_num_leaves(int16_t level, uint64_t node) {
        return Q->get_num_leaves(level, node);
    }

    /**
     * Get the range of leaves in the last level of the tree that are descendants of the node.
     * Useful for the range Maximum query
     * @param level
     * @param node
     * @param init will be modified if the node is not the root. -1 if the node is empty.
     * @param fin will be modified if the node is not the root. -1 if the node is empty.
     */
    void get_range_leaves(int16_t level, uint64_t node, int64_t& init, int64_t& end){
        return Q->get_range_leaves(level, node, init, end);
    }

    /**
     *
     * @param level of the parent
     * @param parent
     * @param child the i-th get_child_se_quadtree of the parent
     * @return the number of leaves of the i-th get_child_se_quadtree of the parent node.
     * @example:
     * 0101
     * 0100 1001
     * 1111 0001 1000
     * get_num_leaves_ith_node(-1,0,0) --> 6
     * get_num_leaves_ith_node(0,0,0) --> 0
     * get_num_leaves_ith_node(0,1,0) --> 1
     * get_num_leaves_ith_node(0,2,0) --> 0 // no existen 3 nodos en el nivel 1
     */
    uint64_t get_num_leaves_ith_node(int16_t level, uint64_t parent, uint64_t child) {
        return Q->get_num_leaves_ith_node(level, parent, child);
    }

    void test_rank(){
        return Q->test_rank();
    }

    void test_get_4_bits(){
        return Q->test_get_4_bits();
    }

    /**
     *
     * @return 1 if the grid is a single point, 0 if the grid is empty, 0.5 otherwise
     */
    double value() {
        switch (this->Q->total_ones_level(this->getHeight() - 1)){
            case 0: // leaf
                return 0;
            case 1:
                return 1;
            default:
                return 0.5;
        }
    }

    /**
     * Algorithm 2 Child(Q,i)
     * Return a qdag corresponding to the i-th get_child_se_quadtree of Q
     * @param Q a qdag Q = (Q', M)
     * @param i
     * @return
     */
    void get_child_qdag(uint64_t level,uint16_t i){
        /*qdag *q = new qdag();
        q->Q = this->Q->get_child_se_quadtree(level,i);
        q->M = this->M;
        q->attribute_set = this->attribute_set;
        q->grid_side = this->grid_side;
        q->is_extended_qdag = this->is_extended_qdag;
        q->Msize = this->Msize;
        return q;*/
        //TODO: su esta extendido, hay que hacer get(M[i]) (el mapeo)
        uint64_t* test = this->Q->get_child_se_quadtree(level,i);
        cout << "get_child_qdag" << endl;
        for(uint16_t j=0;j< getHeight()-1;j++){
            cout << test[j] << " ";
        }

    }

};

#endif
