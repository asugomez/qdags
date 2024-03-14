//
// Created by Asunción Gómez on 14-03-24.
//

#ifndef QDAGS_SE_QUADTREE_DFUDS_HPP
#define QDAGS_SE_QUADTREE_DFUDS_HPP

#include <tuple>
#include <fstream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/k2_tree_helper.hpp>
#include <sdsl/bp_support.hpp>
//#include "rank.hpp"

using namespace std;
using namespace sdsl;

class se_quadtree_dfuds {
public:
    typedef k2_tree_ns::idx_type idx_type;
    typedef k2_tree_ns::size_type size_type;

private:
    // TODO: comoo funcionará close, etc (operaciones BV) si están separados?
    bp_support_g<> *bp_s; // array S in pre-order. It contains the k^d bits that tell which of the k^d possible children of the node exist.
    bp_support_g<> *bp_b; // array B in pre-order with the description of each node 1^c 0, with c the number of children.

    uint16_t height;   // number of levels of the tree [0,..., height-1]

    uint8_t k; // k_d--tree
    uint8_t d; // number of attributs
    size_type k_d; // k^d
    //vector<uint64_t> total_ones;
    uint64_t ref_count;

protected:

    /*! Get the chunk index ([0, k^d[) of a submatrix point.
   **** TODO: PREGUNTA: como están indexados? cual es el orden
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
    uint16_t get_chunk_idx(vector<idx_type> &point, idx_type *offset,
                           size_type l, uint8_t _k) {
        uint16_t i, _d = point.size();
        uint16_t total_sum = 0, _k_aux = 1;

        for (i = _d - 1; i > 0; --i) {
            total_sum += ((point[i] - offset[i]) / l) * _k_aux;
            _k_aux *= _k;
        }

        total_sum += ((point[0] - offset[0]) / l) * _k_aux;
        return total_sum;
    }


    // TODO: do it xd
    //! Build a space efficient quadtree from a set of d-dimensional points
    /*! This method takes a vector of d-dimensional points
     *  and the grid-side size. It takes linear time over the amount of
     *  points to build the quadtree representation.
     *  \param edges A non-empty vector with the d-dimensional points
     *  \param size side of the grid containing the points, all the point
                    coordinates in vector edges must be within [0, size[
     */
    void build_from_edges(std::vector<std::vector<idx_type>> &edges,
                          const size_type size, uint8_t __k, uint8_t __d) {

        typedef std::tuple<idx_type, idx_type, size_type, idx_type *> t_part_tuple;

        k = __k;
        d = __d;
        height = std::ceil(std::log(size) / std::log(k));
        height = height > 1 ? height : 1; // If size == 0

        //bv = new rank_bv_64[height];
        //total_ones.reserve(height);

        k_d = std::pow(k, d);

        bit_vector k_t_ = bit_vector(k_d, 0);  // OJO, cuidado con esto

        std::queue<t_part_tuple> q;
        idx_type t = 0, last_level = 0;
        idx_type i, j, r_0, c_0, it, c, r, z;
        size_type l = std::pow(k, height - 1);

        std::vector<idx_type> pos_by_chunk(k_d + 1, 0);

        idx_type *top_left_point = new idx_type[d]();
        idx_type *top_left_point_aux;

        q.push(t_part_tuple(0, edges.size(), l, top_left_point));

        size_type cur_l = l, cur_level = 0, n_ones = 0;

        while (!q.empty()) {

            std::vector<idx_type> amount_by_chunk(k_d, 0);
            std::tie(i, j, l, top_left_point) = q.front();
            q.pop();

            if (l != cur_l) {
                cur_l = l;
                k_t_.resize(t);
                //bv[cur_level] = rank_bv_64(k_t_);
                //total_ones[cur_level] = bv[cur_level].n_ones();
                cur_level++;
                t = 0;
                k_t_.resize(0);
                k_t_ = bit_vector(k_d * n_ones, 0);
                n_ones = 0;
            }

            // Get size for each chunk
            for (it = i; it < j; it++)
                amount_by_chunk[get_chunk_idx(edges[it], top_left_point, l, k)] += 1;

            if (l == 1) {
                for (it = 0; it < k_d; it++, t++)
                    if (amount_by_chunk[it] != 0)
                        k_t_[t] = 1;
                    else
                        k_t_[t] = 0;
                // At l == 1 no new elements are enqueued
                continue;
            }

            // Set starting position in the vector for each chunk
            pos_by_chunk[0] = i;

            for (it = 1; it < k_d; it++)
                pos_by_chunk[it] = pos_by_chunk[it - 1] + amount_by_chunk[it - 1];
            // To handle the last case when it = k_d - 1
            pos_by_chunk[k_d] = j;

            // Push to the queue every non zero elements chunk
            for (it = 0; it < k_d; it++, t++)
                // If not empty chunk, set bit to 1
                if (amount_by_chunk[it] != 0) {
                    uint16_t p = std::pow(k, d - 1);
                    idx_type it_aux = it;
                    c = it % k;
                    k_t_[t] = 1;
                    n_ones++;

                    top_left_point_aux = new idx_type[d]();

                    for (z = 0; z < d - 1; ++z) {
                        r = it_aux / p;
                        it_aux -= r * p;
                        p /= k;
                        top_left_point_aux[z] = top_left_point[z] + r * l;
                    }

                    top_left_point_aux[d - 1] = top_left_point[d - 1] + c * l;

                    q.push(t_part_tuple(pos_by_chunk[it],
                                        pos_by_chunk[it + 1],
                                        l / k,
                                        top_left_point_aux));
                } else
                    k_t_[t] = 0;

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
        bv[height - 1] = rank_bv_64(k_t_);

        //total_ones[height - 1] = bv[height - 1].n_ones();

    }

public:

    se_quadtree_dfuds() = default;

    uint64_t size() {
        // TODO: not implemented
        uint64_t i, s = 0;
        /*for (i = 0; i < X; i++){
            s += bp_s[i].size();
        for (i = 0; i < X; i++){
            s += bp_b[i].size();*/

        return s;
    }

    //! Constructor
    /*! This constructor takes a vector of edges describing the graph
     *  and the graph size. It takes linear time over the amount of
     *  edges to build the k_d representation.
     *  \param edges A vector with all the edges of the graph, it can
     *               not be empty.
     *  \param size grid side, all the point coordinates in edges
                    vector must be within 0 and size ([0, size[).
     */

    se_quadtree_dfuds(std::vector<std::vector<idx_type>> &edges,
                        const size_type grid_side, uint8_t __k, uint8_t __d) {
        assert(grid_side > 0);
        assert(edges.size() > 0);

        build_from_edges(edges, grid_side, __k, __d);
        ref_count = 1;
    }

    ~se_quadtree_dfuds() {
        ref_count--;
        if (ref_count == 0){
            delete[] bp_s;
            delete[] bp_b;
        }
    }

    void inc_ref_count() {
        ref_count++;
    }

    uint8_t getK() {
        return k;
    }

    uint8_t getD() {
        return d;
    }

    uint64_t getKD() {
        return k_d;
    }


    uint16_t getHeight() {
        return height;
    }


    // TODO: get_node!

};

#endif //QDAGS_SE_QUADTREE_DFUDS_HPP
