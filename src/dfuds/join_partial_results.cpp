//
// Created by Asunción Gómez on 20-03-24.
//
#include "../qdags_dfuds.hpp"
#include "../utils.hpp"

const uint8_t TYPE_FUN_NUM_LEAVES_DFUDS = 0;
const uint8_t TYPE_FUN_DENSITY_LEAVES_DFUDS = 1;

/**
 *
 * @param Q set of qdags.
 * @param nQ number of Qdags
 * @param max_level of the qdag
 * @param nAtt number of attributes of the final quadtree.
 * @param bounded_result if we have to stop computing the results after a certain number of tuples.
 * @param UPPER_BOUND the number of tuples to compute. Only used if bounded_result is true.
 * @param pq the priority queue to store the tuples.
 * @param type_order_fun 0 num leaves, 1 density. Only taking into account if partial_results is false.
 * @param grid_size the size of the grid. Needed when type_order_fun is 1.
 * @return true if we accomplished succesfully the join. Otherwise, return false.
 */
bool AND_partial_dfuds(qdag_dfuds *Q[], uint16_t nQ, uint64_t max_level, uint64_t nAtt, bool bounded_result, uint64_t UPPER_BOUND,
                 priority_queue<qdagWeight> &pq, uint8_t type_order_fun, uint64_t grid_size, vector<uint16_t*>& results_points) {
    uint64_t p = Q[0]->nChildren(); // number of children of the qdag extended
    uint64_t k_d[nQ];
    uint32_t children;
    uint16_t children_to_recurse[p];
    uint64_t children_to_recurse_size;
    int16_t cur_level;
    uint16_t l = (uint16_t) log2(p); // bits number to define the node's children
    uint64_t results = 0;
    while(!pq.empty()){
        children = 0xffffffff;
        qdagWeight tupleQdags = pq.top();
        pq.pop();
        cur_level = tupleQdags.level;
        // if it's a leaf, output the point coordenates
        uint64_t rank_vector[nQ][64];
        for (uint64_t i = 0; i < nQ && children; ++i){
            k_d[i] = Q[i]->getKD(); // k^d del i-esimo qdag
            if (nAtt == 3) {
                if (cur_level == max_level) {
                    children &= Q[i]->materialize_node_3_lastlevel(tupleQdags.roots[i]);
                } else {
                    children &= Q[i]->materialize_node_3(tupleQdags.roots[i], rank_vector[i]);
                }
            }
            else if (nAtt == 4) {
                if(cur_level == max_level){
                    children &= Q[i]->materialize_node_4_lastlevel(tupleQdags.roots[i]);
                }
                else {
                    children &= Q[i]->materialize_node_4(tupleQdags.roots[i], rank_vector[i]);
                }
            }
            else if (nAtt == 5) {
                if(cur_level == max_level){
                    children &= Q[i]->materialize_node_5_lastlevel(tupleQdags.roots[i]);
                }
                else {
                    children &= Q[i]->materialize_node_5(tupleQdags.roots[i], rank_vector[i]);
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

        for (i = 0; i < children_to_recurse_size; ++i) {
            uint64_t* root_temp = new uint64_t[nQ];
            uint16_t* coordinatesTemp = new uint16_t[l];
            // ex: 1011 --> children_to_recurse = [0,2,3]
            child = children_to_recurse[i]; // the index of children where we have 1

            for(uint16_t k = 0; k < l; k++)
                coordinatesTemp[k] = tupleQdags.coordinates[k];
            transformCoordinates(coordinatesTemp, l, diff_level, child);

            // compute the coordinates if it's a leaf
            if(cur_level == max_level){
                results_points.push_back(coordinatesTemp);
                delete[] root_temp;
                if(bounded_result && ++results >= UPPER_BOUND)
                    return true;
            } else{
                // compute the weight of the tuple (ONLY if it's not a leaf)
                double total_weight = DBL_MAX;
                // calculate the weight of the tuple
                for (uint64_t j = 0; j < nQ; j++) {
                    // we store the parent node that corresponds in the original quadtree of each qdag
                    root_temp[j] = (rank_vector[j][Q[j]->getM(child)]);
                    uint64_t n_leaves_child_node = Q[j]->get_num_leaves(root_temp[j]);
                    if (n_leaves_child_node < total_weight) {
                        total_weight = n_leaves_child_node;
                    }
                }
                if(type_order_fun == TYPE_FUN_DENSITY_LEAVES_DFUDS) // density estimator, otherwise it's the number of leaves (min of the tuple)
                    total_weight /= grid_size;
                // insert the tuple
                qdagWeight this_node = {next_level, root_temp, total_weight, coordinatesTemp} ;
                pq.push(this_node); // add the tuple to the queue
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
 * @param type_order_fun 0 num leaves, 1 density. (n1+n2+...+nn) or (n1/100+n2/100+...+nn/100)
 * @param grid_size the size of the grid. Needed when type_order_fun is 1.
 * @param size_queue the size of the priority queue. -1 if we don't want to limit the size of the queue.
 * @return
 */
bool multiJoinPartialResultsDfuds(vector<qdag_dfuds> &Q, bool bounded_result, uint64_t UPPER_BOUND, uint8_t type_order_fun,
                             uint64_t grid_size, vector<uint16_t*>& results_points) {
    qdag_dfuds::att_set A;
    map<uint64_t, uint8_t> attr_map;
    //iterar por el vector de los qdags
    // computes the union of the attribute sets
    // cad auno de los qdags guarda los atributos q le corresponden (con getAttr los retorna, y con nAttr te devuelve el nro de los atributos)
    for (uint64_t i = 0; i < Q.size(); i++) {
        uint64_t nAttr = Q[i].nAttr();
        for (uint64_t j = 0; j < nAttr; j++)
            attr_map[Q[i].getAttr(j)] = 1;
    }

    for (map<uint64_t, uint8_t>::iterator it = attr_map.begin(); it != attr_map.end(); it++)
        A.push_back(it->first); // conjunto de atributos A

    qdag_dfuds *Q_star[Q.size()]; // arreglo de qdags Q start
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
        Q_roots[i] = 3; // root of every qdag in preorder!!
    }


    // ---------------  everything is the same up to here --------------- //

    uint64_t max_level = Q_star[0]->getHeight() - 1;

    uint16_t coordinates[A.size()];
    for(uint16_t i = 0; i < A.size(); i++)
        coordinates[i] = 0;

    priority_queue<qdagWeight> pq; // maxHeap
    pq.push({0, Q_roots, 1, coordinates}); // insert the root of the qda

    AND_partial_dfuds(Q_star, Q.size(),
                max_level, A.size(),
                bounded_result, UPPER_BOUND,
                pq, type_order_fun, grid_size, results_points);


    cout << "number of results: " << results_points.size() << endl;
    uint64_t i=0;
    while(i < results_points.size()){
        uint16_t* coordinates = results_points[i];
        for(uint64_t j = 0; j< A.size(); j++){
            cout << coordinates[j] << " ";
        }
        cout << endl;
        i++;
    }

    for (uint64_t i = 0; i < Q.size(); i++)
        delete Q_star[i];

    return true;
}