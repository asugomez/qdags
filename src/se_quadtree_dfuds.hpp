//
// Created by Asunción Gómez on 14-03-24.
//

#ifndef QDAGS_SE_QUADTREE_DFUDS_HPP
#define QDAGS_SE_QUADTREE_DFUDS_HPP

#include <tuple>
#include <fstream>
#include "../include/bp_support_sada_v2.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/k2_tree_helper.hpp>
#include "rank.hpp"

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
    rank_bv_64 bv_s;
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
        bit_vector k_t_s = bit_vector(k_d, 0);

		// init 110
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
                k_t_s.bit_resize(size_bv_s+t);
                k_t_s.bit_resize(size_bv_s+k_d);
                k_t_s.set_int(size_bv_s, * k_t_.data(), t);


                size_bv_s += t;

                // write 1^c0
                if(!is_leaves){
                    n_children = __builtin_popcount(*k_t_.data());//t/k_d;
                    k_t_b.bit_resize(size_bv_b + n_children  + 1);
                    for(index = size_bv_b; index < size_bv_b + n_children; index++){
                        k_t_b[index] = 1;
                    }
                    k_t_b[index] = 0;
                    size_bv_b+= n_children + 1;
                } else {
                    // separar en chunks de k_t=4
                    uint64_t n_leaves = t/k_d;
                    for(uint64_t leaf = 0; leaf < n_leaves; leaf++){
                        n_children = __builtin_popcount(k_t_.get_int(leaf*k_d,k_d));//t/k_d;
                        k_t_b.bit_resize(size_bv_b + 2*n_children  + 1);
                        for(index = size_bv_b; index < size_bv_b + n_children; index++){
                            k_t_b[index] = 1;
                        }
                        k_t_b[index] = 0;
                        size_bv_b += n_children + 1;
                        for(index = index+1; index < size_bv_b + n_children; index++){
                            k_t_b[index] = 0;
                        }
                        size_bv_b += n_children;
                    }
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

        // write bv S
        k_t_.resize(t );
        k_t_s.bit_resize(size_bv_s+t);
        k_t_s.set_int(size_bv_s, * k_t_.data(), t);

        // write bv B
        // write 1^c 0
        uint64_t n_leaves = t/k_d;
        for(uint64_t leaf = 0; leaf < n_leaves; leaf++){
            n_children = __builtin_popcount(k_t_.get_int(leaf*k_d,k_d));//t/k_d;
            k_t_b.bit_resize(size_bv_b + 2*n_children  + 1);
            for(index = size_bv_b; index < size_bv_b + n_children; index++){
                k_t_b[index] = 1;
            }
            k_t_b[index] = 0;
            size_bv_b += n_children + 1;
            for(index = index+1; index < size_bv_b + n_children; index++){
                k_t_b[index] = 0;
            }
            size_bv_b += n_children;
        }

        // bitvectors
        bv_s = rank_bv_64(k_t_s);
        bit_vector_b = new bit_vector(k_t_b);

        // construct bp
        bp_b = asu::bp_support_sada_v2<>(bit_vector_b);
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
        //bit_vector_b
        //bp_b
        //bv_s
        if(ref_count == 0){
            //delete bit_vector_b;
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

    inline uint32_t get_node_bits(uint64_t start_pos) {
        return bv_s.get_kd_bits(start_pos, k_d);
    }


    /// ----------------- Operations from Compact Data Structures (pg. 341) ----------------- ///


//    size_type_bv find_open(size_type_bv node_v)const{
//        return bp_b.find_open(node_v);
//    }

    size_type_bv find_close(size_type_bv node_v)const{
        return bp_b.find_close(node_v);
    }

    size_type_bv root()const{
        return 3;
    }

    uint64_t rank_one(size_type_bv node_v){
        return bv_s.rank(node_v);
    }

    // suc_0(B,v)
    size_type_bv succ_zero(size_type_bv node_v) const{
        return bp_b.succ_zero(node_v);
    }

    /**
     *
     * @param node_v
     * @return The next sibling of node v, if it exists.
     */
    size_type_bv next_sibling(size_type_bv node_v)const{
        // B[open(B,v-1) - 1] == 1
        assert(bp_b.is_open(bp_b->find_open(node_v-1) - 1));
        return bp_b.next_sibling(node_v);
    }

    /**
     *
     * @param node_v
     * @return The previous sibling of node v, if it exists.
     */
    size_type_bv preceding_sibling(size_type_bv node_v)const{
        // B[v-2, v-1] != [1,0]
        assert(!bp_b.is_open(node_v - 2)); // 1
        assert(bp_b.is_open(node_v - 1)); // 0
        return bp_b.find_close(bp_b.find_open(node_v-1) + 1) + 1;
    }

    // pred_0(B,v)
    size_type_bv pred_zero(size_type_bv node_v) const{
        // a partir de una posición, encontrar la posición de la anterior ocurrencia de 0
        return bp_b.pred_zero(node_v);
    }


    /**
     *
     * @param node_v
     * @return Number of nodes in the subtree of node_v, counting node_v.
     */
    size_type_bv subtree(size_type_bv node_v)const{
        // (fwdsearch(B, v-1, -1) - v)/2 + 1
        assert(node_v >= 3);
        return bp_b.subtree(node_v);
    }

    /**
     *
     * @param node_v
     * @return Number of children of node v.
     */
    size_type_bv children(size_type_bv node_v) const{
        assert(node_v >= 3);
        return succ_zero(node_v) - node_v;
    }

    /**
     *
     * @param node_v
     * @param t the t-th child (first child should be t=1)
     * @return The t-th child of node v, if it exists
     * Asuming t <= children(node_v)
     */
    size_type_bv child(size_type_bv node_v, size_type_bv t) const{
        return bp_b.find_close(succ_zero(node_v) - t) + 1;
    }

    /**
     *
     * @param node_v
     * @return The t such that node v is the tth child of its parent.
     */
    size_type_bv childrank(size_type_bv node_v) const{
        size_type_bp p = bp_b.find_open(node_v - 1);
        return succ_zero(p) - p;
    }


    /**
     * @param node_v
     * @return The number of leaves to the left of node v, plus 1.
     */
    size_type_bv leaf_rank(size_type_bv node_v) const{
        return bp_b.rank_zero_zero(node_v);
    }

    /**
     *
     * @param node_v
     * @return The number of leaves in the subtree of node v.
     */
    size_type_bv leaf_num(size_type_bv node_v) const{
        return leaf_rank(bp_b.next_sibling(node_v)) - leaf_rank(node_v);
    }

    // rank_0(B,v) (de B[0,v-1])
    size_type_bv rank_zero(size_type_bv node_v) const{
        return node_v - bp_b.rank(node_v);
    }

    size_type_bv first_child(size_type_bv node_v)const{
        assert(bp_b.is_open(node_v));
        return succ_zero(node_v) + 1;
    }

    size_type_bv last_child(size_type_bv node_v)const{
        assert(bp_b.is_open(node_v));
        return bp_b.find_close(node_v) + 1;
    }

    size_type_bv parent(size_type_bv node_v)const{
        assert(node_v != 3);
        return pred_zero(bp_b.find_open(node_v - 1)) + 1;
    }

    bool is_leaf(size_type_bv node_v)const{
        return bit_vector_b->get_int(node_v-1,2) == 0;
    }

    // ---- other operations for join ---- //

    /**
     *
     * @param node the starting position in the bitvector B in preorder
     * @param rank_array
     * @return
     */
    inline uint8_t get_node(uint64_t node, uint64_t *rank_array) {
        uint64_t start_pos = (bp_b.rank_zero(node) - 1 - leaf_rank(node))*k_d;
        uint8_t nd;
        // node + n_children + 1 (1^c 0)
        node = first_child(node);
        if(k_d == 4){
            nd = bv_s.get_4_bits(start_pos);
            switch (nd) {
                case 0:
                    break;
                case 1:
                    rank_array[0] = node;
                    break;
                case 2:
                    rank_array[1] = node;
                    break;
                case 3:
                    rank_array[0] = node;
                    node = next_sibling(node);
                    rank_array[1] = node;
                    break;
                case 4:
                    rank_array[2] = node;
                    break;
                case 5:
                    rank_array[0] = node;
                    node = next_sibling(node);
                    rank_array[2] = node;
                    break;
                case 6:
                    rank_array[1] = node;
                    node = next_sibling(node);
                    rank_array[2] = node;
                    break;
                case 7:
                    rank_array[0] = node;
                    node = next_sibling(node);
                    rank_array[1] = node;
                    node = next_sibling(node);
                    rank_array[2] = node;
                    break;
                case 8:
                    rank_array[3] = node;
                    break;
                case 9:
                    rank_array[0] = node;
                    node = next_sibling(node);
                    rank_array[3] = node;
                    break;
                case 10:
                    rank_array[1] = node;
                    node = next_sibling(node);
                    rank_array[3] = node;
                    break;
                case 11:
                    rank_array[0] = node;
                    node = next_sibling(node);
                    rank_array[1] = node;
                    node = next_sibling(node);
                    rank_array[3] = node;
                    break;
                case 12:
                    rank_array[2] = node;
                    node = next_sibling(node);
                    rank_array[3] = node;
                    break;
                case 13:
                    rank_array[0] = node;
                    node = next_sibling(node);
                    rank_array[2] = node;
                    node = next_sibling(node);
                    rank_array[3] = node;
                    break;
                case 14:
                    rank_array[1] = node;
                    node = next_sibling(node);
                    rank_array[2] = node;
                    node = next_sibling(node);
                    rank_array[3] = node;
                    break;
                case 15:
                    rank_array[0] = node;
                    node = next_sibling(node);
                    rank_array[1] = node;
                    node = next_sibling(node);
                    rank_array[2] = node;
                    node = next_sibling(node);
                    rank_array[3] = node;
                    break;
            }
        }
        else {
            nd = bv_s.get_2_bits(start_pos);
            switch (nd) {
                case 0:
                    break;
                case 1:
                    rank_array[0] = node;
                    break;
                case 2:
                    rank_array[1] = node;
                    break;
                case 3:
                    rank_array[0] = node;
                    node = next_sibling(node);
                    rank_array[1] = node;
                    break;
            }
        }
        return nd;
    }

    /**
     *
     * @param node
     * @return The bits of the node of the tree that are descendants of the node.
     */
    inline uint8_t get_node_last_level(uint64_t node) {
        uint64_t start_pos = (bp_b.rank_zero(node) - 1 - leaf_rank(node))*k_d;
        if(k_d == 4){
            return bv_s.get_4_bits(start_pos);
        } else {
            return bv_s.get_2_bits(start_pos);
        }
    }

    /**
     * Get the range of leaves in the last level of the tree that are descendants of the node.
     * Useful for the range Maximum query
     * @param node
     * @param init will be modified if the node is not the root. -1 if the node is empty.
     * @param fin will be modified if the node is not the root. -1 if the node is empty.
     */
    bool get_range_leaves(uint64_t node, uint64_t& init, uint64_t& end){
        init = leaf_rank(node);
        end = leaf_rank(next_sibling(node)) - 1;
        return init <= end;
    }

    /**
     *
     * @param node
     * @param t the position of the child
     * @return the rank of the t-th child of node.
     * If the node is 1001, and t = 3, the index will be 1
     */
    uint64_t get_rank_node_child(uint64_t node, uint64_t t){
        return bp_b.rank(node+t) - bp_b.rank(node);
    }


};
#endif //QDAGS_SE_QUADTREE_DFUDS_HPP
