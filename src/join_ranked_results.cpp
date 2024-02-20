#include <algorithm>
#include "qdags.hpp"
#include "./MinMaxHeap.hpp"
//#include <sdsl/rmq_support_sparse_table.hpp>
#include <sdsl/rmq_support.hpp>
#include "utils.hpp"

using namespace sdsl;
using namespace std;

// TODO declare const int 0,1,2 for the type fun
// TODO: see improvement in root and root_temp
/**
 *
 * @param Q set of qdags.
 * @param nQ number of Qdags.
 * @param max_level of the qdag.
 * @param nAtt number of attributes of the final quadtree.
 * @param bounded_result if we have to stop computing the results after a certain number of tuples.
 * @param UPPER_BOUND the number of tuples to compute. Only used if bounded_result is true.
 * @param pq the priority queue to store the tuples.
 * @param type_priority_fun 0 sum , 1 maximum. Only taking into account if partial_results is true.
 * @param priorities the priorities of the points.
 * @return true if we accomplished succesfully the join. Otherwise, return false.
 */
bool AND_ranked(qdag *Q[], uint16_t nQ,
                 uint64_t max_level, uint64_t nAtt, bool bounded_result, uint64_t UPPER_BOUND,
                 priority_queue<qdagWeight> &pq, uint8_t type_priority_fun, int_vector<>& priorities ) {
    uint64_t p = Q[0]->nChildren(); // number of children of the qdag extended
    uint64_t k_d[nQ];
    uint32_t children;
    uint16_t children_to_recurse[p];
    uint64_t children_to_recurse_size;
    uint16_t cur_level;
    uint16_t l = (uint16_t) log2(p); // bits number to define the node's children
    uint64_t results = 0;
    while(!pq.empty()){
        children = 0xffffffff;
        qdagWeight tupleQdags = pq.top();
        pq.pop();
        cur_level = tupleQdags.level;
        uint64_t roots[nQ];
        // HERE we have a problem
        for(uint64_t i = 0; i < nQ; i++){
            roots[i] = tupleQdags.roots[i];
        }
        // if it's a leaf, output the point coordenates
        uint64_t rank_vector[nQ][64];
        for (uint64_t i = 0; i < nQ && children; ++i){
            k_d[i] = Q[i]->getKD(); // k^d del i-esimo qdag
            if (nAtt == 3) {
                if (cur_level == max_level) {
                    children &= Q[i]->materialize_node_3_lastlevel(cur_level, roots[i]);
                } else {
                    children &= Q[i]->materialize_node_3(cur_level, roots[i], rank_vector[i]);
                }
            }
            else if (nAtt == 4) {
                if(cur_level == max_level){
                    children &= Q[i]->materialize_node_4_lastlevel(cur_level, roots[i]);
                }
                else {
                    children &= Q[i]->materialize_node_4(cur_level, roots[i], rank_vector[i]);
                }
            }
            else if (nAtt == 5) {
                if(cur_level == max_level){
                    children &= Q[i]->materialize_node_5_lastlevel(cur_level, roots[i]);
                }
                else {
                    children &= Q[i]->materialize_node_5(cur_level, roots[i], rank_vector[i]);
                }
            }
            else {
                cout << "Code only works for queries of up to 5 attributes..." << endl;
                return false;
            }
        }

        // in children we have the actual non empty nodes
        // in children_to_recurse we store the index of these non emtpy nodes
        children_to_recurse_size = bits::cnt((uint64_t) children);
        uint64_t i = 0;
        uint64_t msb; // most significant bit
        while (i < children_to_recurse_size) {
            msb = __builtin_clz(children);
            children_to_recurse[i] = msb;
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1));
        }

        uint16_t child;
        uint16_t diff_level = max_level-cur_level;
        uint16_t next_level = cur_level + 1;
        uint64_t path;

        for (i = 0; i < children_to_recurse_size; ++i) {
            uint64_t* root_temp= new uint64_t[nQ]; // TODO: when to liberar la memoria
            child = children_to_recurse[i];
            path = 0;

            // compute the weight of the tuple (ONLY if it's not a leaf)
            double total_weight = 0;
            if(cur_level != max_level) {
                // calculate the weight of the tuple
                for (uint64_t j = 0; j < nQ; j++) {
                    // we store the parent node that corresponds in the original quadtree of each qdag
                    root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);

                    // TODO: do the ranked results
                    uint64_t priority_ith_node =0; // TODO_ define priority ith node (rMq)
                    if (type_priority_fun == 0) // sum
                        total_weight += priority_ith_node;
                    else if (type_priority_fun == 1) { // max
                        if (total_weight < priority_ith_node) {
                            total_weight = priority_ith_node;
                        }
                    }
                }
            }
            // --> add the child to the path
            //uint16_t cur_child_qdag = Q[0]->getM(child);
            path = child << (diff_level * l); // height * bits to represent the children
            path += tupleQdags.path; // add the bits to the bitvector
            //tupleQdags.path += path; // add the bits to the bitvector
            // compute the coordinates if it's a leaf
            if(cur_level == max_level){
                uint32_t coordinates[nAtt];
                for(uint32_t k = 0; k < nAtt; k++){
                    coordinates[k] = 0;
                }
                getCoordinates(path, l, max_level, coordinates);
                cout << "results nro ° " << results << endl;
                cout << "point output: " << path << endl;
                for(uint64_t k = 0; k < nAtt; k++){
                    cout << " coord " << k << " = " << coordinates[k] << endl;
                }
                if(bounded_result && ++results >= UPPER_BOUND)
                    return true;
            }
            else{ // insert the tuple
                // TODO_ pq no se pide memoria?
                qdagWeight this_node = {next_level, root_temp, total_weight, path} ;
                pq.push(this_node); // add the tuple to the queue
            }
        }
    }
    return true;
}

/**
 * Computing the join in a certain order. After computing size_queue results, we will stop computing the results.
 * @param Q
 * @param nQ number of Qdags
 * @param max_level of the qdag
 * @param nAtt number of attributes of the final quadtree.
 * @param pq the priority queue to store the tuples.
 * @param size_queue the size of the priority queue.
 * @param type_priority_fun 0 sum , 1 maximum. Only taking into account if partial_results is true.
 * @param priorities the priorities of the points.
 * @return true if we accomplished succesfully the join. Otherwise, return false.
 */
bool AND_ranked_fixed_queue(qdag *Q[], uint16_t nQ,
                 uint64_t max_level, uint64_t nAtt,
                 minmax::MinMaxHeap<qdagWeight>& pq, uint64_t size_queue,
                 uint8_t type_priority_fun, int_vector<>& priorities) {
    uint64_t p = Q[0]->nChildren(); // number of children of the qdag extended
    uint64_t k_d[nQ];
    uint32_t children;
    uint16_t children_to_recurse[p];
    uint64_t children_to_recurse_size;
    uint16_t cur_level;
    uint16_t l = (uint16_t) log2(p); // bits number to define the node's children
    uint64_t min_value_queue = INT_MAX;
    uint64_t results = 0;
    while(!pq.empty()){
        children = 0xffffffff;
        qdagWeight tupleQdags = pq.popMax();

        // update the queue properties
        if(min_value_queue == tupleQdags.weight){ // TODO: see when we pop the last element of the pq
            min_value_queue = pq.size() ? pq.findMin().weight : INT_MAX;
        }

        cur_level = tupleQdags.level;
        uint64_t roots[nQ];
        // HERE we have a problem
        for(uint64_t i = 0; i < nQ; i++){
            roots[i] = tupleQdags.roots[i];
        }
        // if it's a leaf, output the point coordenates
        uint64_t rank_vector[nQ][64];
        for (uint64_t i = 0; i < nQ && children; ++i){
            k_d[i] = Q[i]->getKD(); // k^d del i-esimo qdag
            if (nAtt == 3) {
                if (cur_level == max_level) {
                    children &= Q[i]->materialize_node_3_lastlevel(cur_level, roots[i]);
                } else {
                    children &= Q[i]->materialize_node_3(cur_level, roots[i], rank_vector[i]);
                }
            }
            else if (nAtt == 4) {
                if(cur_level == max_level){
                    children &= Q[i]->materialize_node_4_lastlevel(cur_level, roots[i]);
                }
                else {
                    children &= Q[i]->materialize_node_4(cur_level, roots[i], rank_vector[i]);
                }
            }
            else if (nAtt == 5) {
                if(cur_level == max_level){
                    children &= Q[i]->materialize_node_5_lastlevel(cur_level, roots[i]);
                }
                else {
                    children &= Q[i]->materialize_node_5(cur_level, roots[i], rank_vector[i]);
                }
            }
            else {
                cout << "Code only works for queries of up to 5 attributes..." << endl;
                return false;
            }
        }

        // in children we have the actual non empty nodes
        // in children_to_recurse we store the index of these non emtpy nodes
        children_to_recurse_size = bits::cnt((uint64_t) children);
        uint64_t i = 0;
        uint64_t msb; // most significant bit
        while (i < children_to_recurse_size) {
            msb = __builtin_clz(children);
            children_to_recurse[i] = msb;
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1));
        }

        uint16_t child;
        uint16_t diff_level = max_level-cur_level;
        uint16_t next_level = cur_level + 1;
        uint64_t path;

        for (i = 0; i < children_to_recurse_size; ++i) {
            uint64_t* root_temp= new uint64_t[nQ];
            child = children_to_recurse[i];
            path = 0;

            // compute the weight of the tuple (ONLY if it's not a leaf)
            double total_weight = 0;
            if(cur_level != max_level) {
                // calculate the weight of the tuple
                for (uint64_t j = 0; j < nQ; j++) {
                    // we store the parent node that corresponds in the original quadtree of each qdag
                    root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);
                    // TODO: do the ranked results
                    uint64_t priority_ith_node =0; // TODO_ define priority ith node (rMq)
                    if (type_priority_fun == 0) // sum
                        total_weight += priority_ith_node;
                    else if (type_priority_fun == 1) { // max
                        if (total_weight < priority_ith_node) {
                            total_weight = priority_ith_node;
                        }
                    }

                }
            }
            // --> add the child to the path
            path = child << (diff_level * l); // height * bits to represent the children
            path += tupleQdags.path; // add the bits to the bitvector
            // compute the coordinates if it's a leaf
            if(cur_level == max_level){
                // DEBUG: see the tupleQdags and path
                uint32_t coordinates[nAtt];
                for(uint32_t k = 0; k < nAtt; k++){
                    coordinates[k] = 0;
                }
                getCoordinates(path, l, max_level, coordinates);
                cout << "results nro ° " << results << endl;
                cout << "point output: " << path << endl;
                for(uint64_t k = 0; k < nAtt; k++){
                    cout << " coord " << k << " = " << coordinates[k] << endl;
                }
                if(++results >= size_queue) // we've already output the top size_queue results
                    return true;
            }
            else{ // insert the tuple
                if(pq.size() < size_queue) {
                    qdagWeight this_node = {next_level, root_temp, total_weight, path};
                    pq.push(this_node); // add the tuple to the queue
                    if(total_weight < min_value_queue){ // update the min value
                        min_value_queue = total_weight;
                    }
                } else { // the queue is full
                    if(total_weight > min_value_queue){ // this tuple has a higher priority --> need to be inserted
                        // pop the least value of the queue
                        pq.popMin();
                        // update the min_value_queue
                        qdagWeight minTuple = pq.findMin();
                        min_value_queue = minTuple.weight;
                        // add the new tuple
                        qdagWeight this_node = {next_level, root_temp, total_weight, path};
                        pq.push(this_node); // add the tuple to the queue
                    }
                }
            }
        }
    }
    return true;
}


/**
 *
 * @param Q vector of qdags
 * @param bounded_result if we have to stop computing the results after a certain number of tuples.
 * @param UPPER_BOUND the number of tuples to compute. Only used if bounded_result is true.
 * @param type_priority_fun 0 sum , 1 maximum. (p1+p2+...+pn) or max(p1,p2,...,pn)
 * @param size_queue the size of the priority queue. -1 if we don't want to limit the size of the queue.
 * @param priorities the priorities of the points.
 * @return
 */
bool multiJoinRankedResults(vector<qdag> &Q, bool bounded_result, uint64_t UPPER_BOUND,
                             uint8_t type_priority_fun, int64_t size_queue, int_vector<>& priorities) {
    qdag::att_set A;
    map<uint64_t, uint8_t> attr_map;
    // iterar por el vector de los qdags
    // computes the union of the attribute sets
    // cada uno de los qdags guarda los atributos q le corresponden (con getAttr los retorna, y con nAttr te devuelve el nro de los atributos)
    for (uint64_t i = 0; i < Q.size(); i++) {
        uint64_t nAttr = Q[i].nAttr();
        for (uint64_t j = 0; j < nAttr; j++)
            attr_map[Q[i].getAttr(j)] = 1;
    }

    for (map<uint64_t, uint8_t>::iterator it = attr_map.begin(); it != attr_map.end(); it++)
        A.push_back(it->first); // conjunto de atributos A

    qdag *Q_star[Q.size()]; // arreglo de qdags Q start
    uint64_t Q_roots[Q.size()]; // arreglo de enteros: mantiene el nodo actual en el q estas en cada uno de los qdags (el arreglo te dice en donde estas)

    for (uint64_t i = 0; i < Q.size(); i++){ // iteras por los qdags de la query y se eextiende cada uno de los qdags
        // preprocesamiento
        Q_star[i] = Q[i].extend(A); // extender los qdags en base al conjunto A
        if (A.size() == 3)
            Q_star[i]->createTableExtend3();
        else if (A.size() == 4)
            Q_star[i]->createTableExtend4();
        else if (A.size() == 5)
            Q_star[i]->createTableExtend5();
        else {
            cout << "Code only works for queries of up to 5 attributes..." << endl;
            exit(1);
        }
        // esta en 0 en todos lados al inicio
        Q_roots[i] = 0; // root of every qdag
    }


    // ---------------  everything is the same up to here --------------- //

    uint64_t max_level = Q_star[0]->getHeight() - 1;

    rmq_succinct_sct<false> rMq = rmq_succinct_sct<false>(&priorities);




    if(size_queue == -1) {
        priority_queue<qdagWeight> pq;
        pq.push({0, Q_roots, 1, 0}); // insert the root of the qdag
        AND_ranked(Q_star, Q.size(),
                    max_level, A.size(),
                    bounded_result, UPPER_BOUND,
                    pq, type_priority_fun, priorities);
    }
    else {
        minmax::MinMaxHeap<qdagWeight> pq;
        pq.push({0, Q_roots, 1, 0}); // insert the root of the qdag
        AND_ranked_fixed_queue(Q_star, Q.size(),
                                max_level, A.size(),
                                pq, size_queue, type_priority_fun, priorities);
    }
    return true;
}

