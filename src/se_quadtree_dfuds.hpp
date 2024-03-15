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
    typedef bit_vector::size_type size_type_bp;

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
    void build_from_edges_dfs(std::vector<std::vector<idx_type>> &edges,
                              const size_type size, uint8_t __k, uint8_t __d) {

        typedef std::tuple<idx_type, idx_type, size_type, idx_type *> t_part_tuple;

        k = __k;
        d = __d;
        height = std::ceil(std::log(size) / std::log(k));
        height = height > 1 ? height : 1; // If size == 0

        //bv = new rank_bv_64[height];
        //total_ones.reserve(height);

        k_d = std::pow(k, d);

        //bit_vector k_t_ = bit_vector(k_d, 0);  //
        //bit_vector k_t_b = bit_vector(k_d, 0);  // 1^c 0
        bit_vector k_t_s = bit_vector(edges.size()* k_d, 0);

        k_t_s[0] = 1;
        k_t_s[1] = 1;
        k_t_s[2] = 0;

        //last_pos !! PROBLEME de tamaño


    }


};
#endif //QDAGS_SE_QUADTREE_DFUDS_HPP
