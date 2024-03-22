//
// Created by Asunción Gómez on 14-03-24.
//

#ifndef QDAGS_SE_QUADTREE_DFUDS_HPP
#define QDAGS_SE_QUADTREE_DFUDS_HPP

#include <tuple>
#include <fstream>
#include "include/bp_support_sada_v2.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/k2_tree_helper.hpp>
//#include "rank.hpp"

using namespace std;
using namespace sdsl;

class se_quadtree_dfuds {
public:
    typedef k2_tree_ns::idx_type idx_type;
    typedef k2_tree_ns::size_type size_type_k2;
    typedef bit_vector::size_type size_type_bp;
    typedef int_vector_size_type size_type_bv;
    typedef bit_vector::difference_type difference_type;

private:
    rank_bv_64 bv_s; // TODO: maybe we can delete this!
    //bp_support_sada<> bp_s; // array S in pre-order. It contains the k^d bits that tell which of the k^d possible children of the node exist.
    const bit_vector* bit_vector_b;
    asu::bp_support_sada_v2<> bp_b; // array B in pre-order with the description of each node 1^c 0, with c the number of children.

    uint16_t height;   // number of levels of the tree [0,..., height-1]

    uint8_t k; // k_d--tree
    uint8_t d; // number of attributs
    size_type_k2 k_d; // k^d
    //vector<uint64_t> total_ones;
    uint64_t ref_count;

protected:

    /*! Get the chunk index ([0, k^d[) of a submatrix point.
   *
   * Gets a point in the global matrix and returns its corresponding chunk
   * in the submatrix specified.
   *
   * \param point vector representing the d dimensional point.
   * \param offset vector with the upper-left point of the global matrix
   * \param l size of the chunk at the submatrix.
   * \param k the k parameter from the k^2 tree.
   * \returns the index of the chunk containing the point at the submatrix.
   */
    uint16_t get_chunk_idx(vector<idx_type> &point, const idx_type *offset,
                           size_type_k2 l, uint8_t _k) {
        uint16_t i, _d = point.size();
        uint16_t total_sum = 0, _k_aux = 1;

        for (i = _d - 1; i > 0; --i) {
            total_sum += ((point[i] - offset[i]) / l) * _k_aux;
            _k_aux *= _k;
        }

        total_sum += ((point[0] - offset[0]) / l) * _k_aux;
        return total_sum;
    }

    // TODO: see if it's necessary to add 110!!

    //! Build a space efficient quadtree from a set of d-dimensional points
    /*! This method takes a vector of d-dimensional points
     *  and the grid-side size. It takes linear time over the amount of
     *  points to build the quadtree representation.
     *  \param edges A non-empty vector with the d-dimensional points
     *  \param size side of the grid containing the points, all the point
                    coordinates in vector edges must be within [0, size[
     */
    void build_from_edges_dfs(std::vector<std::vector<idx_type>> &edges,
                              const size_type_k2 size, uint8_t __k, uint8_t __d) {

        typedef std::tuple<idx_type, idx_type, size_type_k2, idx_type *> t_part_tuple;

        k = __k;
        d = __d;
        height = std::ceil(std::log(size) / std::log(k));
        height = height > 1 ? height : 1; // If size == 0

        k_d = std::pow(k, d);

        bit_vector k_t_ = bit_vector(k_d, 0);

        bit_vector k_t_b = bit_vector(3, 0);  // 1^c 0
        bit_vector k_t_s = bit_vector(k_d, 0); // init 110

        k_t_b[0] = 1;
        k_t_b[1] = 1;
        k_t_b[2] = 0;
        size_type_bv size_bv_b = 3; // starting position to write in the bit_vector
        size_type_bv size_bv_s = 0; // starting position to write in the bit_vector



        std::stack<t_part_tuple> s;
        idx_type t = 0;
        idx_type t_aux;
        idx_type i, j, it, c, r, z;
        size_type_k2 l = std::pow(k, height - 1);

        std::vector<idx_type> pos_by_chunk(k_d + 1, 0);

        idx_type *top_left_point = new idx_type[d]();
        idx_type *top_left_point_aux;

        s.emplace(0, edges.size(), l, top_left_point);

        size_type_k2 cur_l = l, cur_level = 0, n_ones = 0;

        uint64_t n_children;
        uint64_t index;
        bool is_leaves = false;

        while (!s.empty()) {
            std::vector<idx_type> amount_by_chunk(k_d, 0);
            std::tie(i, j, l, top_left_point) = s.top();
            s.pop();

            if (l != cur_l) {
                cur_l = l;
                k_t_.resize(t);
                k_t_s.resize(size_bv_s+t);
                k_t_s.set_int(size_bv_s, * k_t_.data(), t);

                //cout << k_t_ << endl;

                size_bv_s += t;

                // write 1^c 0
                if(!is_leaves){
                    n_children = __builtin_popcount(*k_t_.data());//t/k_d;
                    k_t_b.resize(size_bv_b + n_children  + 1);
                    for(index = size_bv_b; index < size_bv_b + n_children; index++){
                        k_t_b[index] = 1;
                    }
                    k_t_b[index] = 0;
                    size_bv_b+= n_children + 1;
                }
                if(is_leaves){
                    n_children = t/k_d;
                    k_t_b.resize(size_bv_b + n_children );
                    for(index = size_bv_b; index < size_bv_b + n_children; index++){
                        k_t_b[index] = 0;
                    }
                    size_bv_b+= n_children;
                }

                is_leaves = false;
                cur_level++;
                t = 0;
                k_t_.resize(0);
                k_t_ = bit_vector(k_d, 0);//bit_vector(k_d * n_ones, 0);
                n_ones = 0;
            }

            // Get size for each chunk
            for (it = i; it < j; it++)
                amount_by_chunk[get_chunk_idx(edges[it], top_left_point, l, k)] += 1;
            // l = k^h = 1
            if (l == 1) {
                for (it = 0; it < k_d; it++, t++)
                    if (amount_by_chunk[it] != 0)
                        k_t_[t] = 1;
                    else
                        k_t_[t] = 0;
                // At l == 1 no new elements are enqueued
                is_leaves = true;
                continue;
            }
            // Set starting position in the vector for each chunk
            pos_by_chunk[0] = i;

            for (it = 1; it < k_d; it++)
                pos_by_chunk[it] = pos_by_chunk[it - 1] + amount_by_chunk[it - 1];
            // To handle the last case when it = k_d - 1
            pos_by_chunk[k_d] = j;

            // Push to the queue every non zero elements chunk
            t_aux = t+k_d-1;
            //for (it = 0; it < k_d; it++, t++)
            for (it = k_d; it --> 0; t++, t_aux--)
                // If not empty chunk, set bit to 1
                if (amount_by_chunk[it] != 0) {
                    uint16_t p = std::pow(k, d - 1);
                    idx_type it_aux = it;
                    c = it % k;
                    k_t_[t_aux] = 1;
                    n_ones++;

                    top_left_point_aux = new idx_type[d]();

                    for (z = 0; z < d - 1; ++z) {
                        r = it_aux / p;
                        it_aux -= r * p;
                        p /= k;
                        top_left_point_aux[z] = top_left_point[z] + r * l;
                    }

                    top_left_point_aux[d - 1] = top_left_point[d - 1] + c * l;

                    s.emplace(pos_by_chunk[it],
                                        pos_by_chunk[it + 1],
                                        l / k,
                                        top_left_point_aux);
                } else
                    k_t_[t_aux] = 0;

            idx_type chunk;

            // Sort edges' vector
            for (unsigned ch = 0; ch < k_d; ch++) {
                idx_type be = ch == 0 ? i : pos_by_chunk[ch - 1];
                for (it = pos_by_chunk[ch]; it < be + amount_by_chunk[ch];) {
                    chunk = get_chunk_idx(edges[it], top_left_point, l, k);

                    if (pos_by_chunk[chunk] != it)
                        std::iter_swap(edges.begin() + it,
                                       edges.begin() + pos_by_chunk[chunk]);
                    else
                        it++;
                    pos_by_chunk[chunk]++;
                }
            }
            delete[] top_left_point;

        }

        k_t_.resize(t);
        k_t_s.resize(size_bv_s+t);
        k_t_s.set_int(size_bv_s, * k_t_.data(), t);

        //cout << k_t_ << endl;
        // write 1^c 0
        n_children = t/k_d;//
        k_t_b.resize(size_bv_b + n_children );
        for(index = size_bv_b; index < size_bv_b + n_children; index++){
            k_t_b[index] = 0;
        }

        // bitvectors
        bv_s = rank_bv_64(k_t_s); // TODO: fix this!! max 64 bits!!
        bit_vector_b = new bit_vector(k_t_b);

        // TODO: fix problem size uint64_t
        // construct bp
        bp_b = asu::bp_support_sada_v2<>(bit_vector_b);

        cout << "finish" << endl;
    }

public:

    se_quadtree_dfuds() = default;

    uint64_t size() {
        //TODO
        //bp_b->size();// + bp_s->size();
        return 0;
    }

    se_quadtree_dfuds(std::vector<std::vector<idx_type>> &edges,
                      const size_type_k2 grid_side, uint8_t __k, uint8_t __d){
        assert(grid_side > 0);
        assert(edges.size() > 0);

        build_from_edges_dfs(edges, grid_side, __k, __d);
        ref_count = 1;
    }

    ~se_quadtree_dfuds() {
        ref_count--;
        if(ref_count == 0){
            delete bit_vector_b;
            //delete bv_s;
        }
    }

    void inc_ref_count() {
        ref_count++;
    }

    uint8_t getK() const{
        return k;
    }


    uint8_t getD() const{
        return d;
    }

    uint64_t getKD() const{
        return k_d;
    }

    uint64_t getHeight() const{
        return height;
    }

    inline uint8_t get_node_bits(uint64_t start_pos) {
        return bv_s.get_kd_bits(start_pos, k_d);
    }


    /// ----------------- Operations from Compact Data Structures (pg. 341) ----------------- ///

    size_type_bp fwd_search(size_type_bp start_pos, difference_type diff_excess)const{
        return this->bp_b.fwd_search(start_pos, diff_excess);
    }

//    size_type_bp bwd_search(size_type_bp start_pos, difference_type diff_excess)const{
//        return this->bp_b.bwd_search(start_pos, diff_excess);
//    }

//    size_type_bp find_open(size_type_bp node_v)const{
//        return bp_b.find_open(node_v);
//    }

    size_type_bv find_close(size_type_bp node_v)const{
        return bp_b.find_close(node_v);
    }

    size_type_bp root()const{
        return 3;
    }

    uint64_t rank_one(size_type_bv node_v){
        return bv_s.rank(node_v);
    }

    // suc_0(B,v)
    size_type_bp succ_zero(size_type_bp node_v) const{
        return bp_b.succ_zero(node_v);
    }

//    // select_0(B,v)
//    size_type_bp select_zero(size_type_bp node_v){
//        return bp_b.select_zero(node_v); // TODO: select_zero(0) --> 3
//    }

    /*
    // pred_0(B,v)
    size_type_bp pred_zero(size_type_bp node_v) const{
        // a partir de una posición, encontrar la posición de la anterior ocurrencia de 0
        return bp_b.pred_zero(node_v);
    }



    // rank_0(B,v) (de B[0,v-1])
    size_type_bp rank_zero(size_type_bp node_v) const{
        return node_v - bp_b.rank(node_v);
    }

    // rank_00(B,v)
    size_type_bp rank_zero_zero(size_type_bp node_v) const{
        return bp_b.rank_zero_zero(node_v);
    }

     //no need of these operations?


    // select_00(B,v)
    size_type_bp select_zero_zero(size_type_bp node_v){

    }


    size_type_bp first_child(size_type_bp node_v)const{
        assert(bp_b.is_open(node_v));
        return succ_zero(node_v) + 1;
    }

    size_type_bp last_child(size_type_bp node_v)const{
        assert(bp_b.is_open(node_v));
        return bp_b.find_close(node_v) + 1;
    }

    size_type_bp next_sibling(size_type_bp node_v)const{
        // B[open(B,v-1) - 1] == 1
        assert(bp_b.is_open(bp_b->find_open(node_v-1) - 1));
        return fwd_search(node_v-1,-1) + 1;
    }

    size_type_bp preceding_sibling(size_type_bp node_v)const{
        // B[v-2, v-1] != [1,0]
        assert(bp_b.is_open(node_v - 2)); // 1
        assert(! bp_b.is_open(node_v - 1)); // 0
        return bp_b.find_close(bp_b.find_open(node_v-1) + 1) + 1;
    }

    size_type_bp parent(size_type_bp node_v)const{
        assert(node_v != 3);
        return pred_zero(node_v - 1) + 1;
    }

    bool is_leaf(size_type_bp node_v)const{
        return bit_vector_b->get_int(node_v,1) == 0;
    }*/

    // size of the subtree (counting the node_v)
    size_type_bv subtree(size_type_bp node_v)const{
        // (fwdsearch(B, v-1, -1) - v)/2 + 1
        assert(node_v >= 3);
        return (fwd_search(node_v -1, -1) - node_v)/2 + 1;
    }

    size_type_bv children(size_type_bp node_v) const{
        return succ_zero(node_v) - node_v;
    }
/*
    size_type_bv child(size_type_bp node_v, size_type_bp t) const{
        return bp_b.find_close(succ_zero(node_v));// - t) + 1;
    }

    size_type_bv childrank(size_type_bp node_v) const{
        size_type_bv p = bp_b.find_open(node_v - 1);
        return succ_zero(p) - p;
    }


    // numbers of 1s per level til node_v

    size_type_bv leaf_rank(size_type_bp node_v) const{
        return rank_zero_zero(node_v-1) + 1;
    }

    size_type_bv leaf_num(size_type_bp node_v) const{
        return leaf_rank(fwd_search(node_v-1, -1) + 1) - leaf_rank(node_v);
    }
     */


    // ---- other operations for join ---- //
    // TODO: see how it works without level (i dont think we need it)!
    /**
     *
     * @param node the starting position in the bitvector S in preorder
     * @param rank_array
     * @param r the starting position in the bitvector B in preorder
     * @return
     */
    inline uint8_t get_node(uint64_t node, uint64_t *rank_array, uint64_t r) {
        uint64_t start_pos = (bp_b.rank_zero(node) - 1) * k_d;
        uint8_t nd;
        // r + n_children + 1 (1^c 0)
        size_type_bv n_children = children(node);
        r += n_children + 1;
        if(k_d == 4){
            nd = bv_s.get_4_bits(start_pos);
            switch (nd) {
                case 0:
                    break;
                case 1:
                    rank_array[0] = r;
                    break;
                case 2:
                    rank_array[1] = r;
                    break;
                case 3:
                    rank_array[0] = r;
                    r+= subtree(r) + 1;
                    rank_array[1] = r;
                    break;
                case 4:
                    rank_array[2] = r;
                    break;
                case 5:
                    rank_array[0] = r;
                    r+= subtree(r) + 1;
                    rank_array[2] = r;
                    break;
                case 6:
                    rank_array[1] = r;
                    r+= subtree(r) + 1;
                    rank_array[2] = r;
                    break;
                case 7:
                    rank_array[0] = r;
                    r+= subtree(r) + 1;
                    rank_array[1] = r;
                    r+= subtree(r) + 1;
                    rank_array[2] = r;
                    break;
                case 8:
                    rank_array[3] = r;
                    break;
                case 9:
                    rank_array[0] = r;
                    r+= subtree(r) + 1;
                    rank_array[3] = r;
                    break;
                case 10:
                    rank_array[1] = r;
                    r+= subtree(r) + 1;
                    rank_array[3] = r;
                    break;
                case 11:
                    rank_array[0] = r;
                    r+= subtree(r) + 1;
                    rank_array[1] = r;
                    r+= subtree(r) + 1;
                    rank_array[3] = r;
                    break;
                case 12:
                    rank_array[2] = r;
                    r+= subtree(r) + 1;
                    rank_array[3] = r;
                    break;
                case 13:
                    rank_array[0] = r;
                    r+= subtree(r) + 1;
                    rank_array[2] = r;
                    r+= subtree(r) + 1;
                    rank_array[3] = r;
                    break;
                case 14:
                    rank_array[1] = r;
                    r+= subtree(r) + 1;
                    rank_array[2] = r;
                    r+= subtree(r) + 1;
                    rank_array[3] = r;
                    r+= subtree(r) + 1;
                    break;
                case 15:
                    rank_array[0] = r;
                    r+= subtree(r) + 1;
                    rank_array[1] = r;
                    r+= subtree(r) + 1;
                    rank_array[2] = r;
                    r+= subtree(r) + 1;
                    rank_array[3] = r;
                    break;
            }
        } else {
            nd = bv_s.get_2_bits(start_pos);
            switch (nd) {
                case 0:
                    break;
                case 1:
                    rank_array[0] = r;
                    break;
                case 2:
                    rank_array[1] = r;
                    break;
                case 3:
                    rank_array[0] = r;
                    r+= subtree(r) + 1;
                    rank_array[1] = r;
                    break;
            }
        }

        return nd;
    }

    inline uint8_t get_node_last_level(uint64_t node) {
        if(k_d == 4){
            return bv_s.get_4_bits(node);
        } else {
            return bv_s.get_2_bits(node);
        }
    }




};
#endif //QDAGS_SE_QUADTREE_DFUDS_HPP
