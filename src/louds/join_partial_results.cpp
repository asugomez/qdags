#include <algorithm>
#include "../qdags.hpp"
#include "../utils.hpp"

const uint8_t TYPE_FUN_NUM_LEAVES = 0;
const uint8_t TYPE_FUN_DENSITY_LEAVES = 1;


/**
 *
 * @param Q set of qdags.
 * @param nQ number of Qdags
 * @param nAtt number of attributes of the final quadtree.
 * @param bounded_result if we have to stop computing the results after a certain number of tuples.
 * @param UPPER_BOUND the number of tuples to compute. Only used if bounded_result is true.
 * @param max_level of the qdag
 * @param grid_size the size of the grid. Needed when type_order_fun is 1.
 * @param type_order_fun 0 num leaves, 1 density. Only taking into account if partial_results is false.
 * @param pq the priority queue to store the tuples.
 * @param results_points vector with the coordinates of the output.
 * @return true if we accomplished succesfully the join. Otherwise, return false.
 */
bool AND_partial(
        qdag *Q[],
        uint16_t nQ,
        uint64_t nAtt,
        bool bounded_result,
        uint64_t UPPER_BOUND,
        uint64_t max_level,
        uint64_t grid_size,
        uint8_t type_order_fun,
        priority_queue<qdagWeight> &pq,
        vector<uint16_t *> &results_points) {

    uint64_t p = Q[0]->nChildren(); // number of children of the qdag extended
    uint64_t k_d[nQ];
    uint32_t children;
    uint16_t children_to_recurse[p];
    uint64_t children_to_recurse_size;
    uint16_t cur_level;
    uint16_t l = (uint16_t) log2(p); // bits number to define the quadtree_formula's children
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
                    children &= Q[i]->materialize_node_3_lastlevel(cur_level, tupleQdags.roots[i]);
                } else {
                    children &= Q[i]->materialize_node_3(cur_level, tupleQdags.roots[i], rank_vector[i]);
                }
            }
            else if (nAtt == 4) {
                if(cur_level == max_level){
                    children &= Q[i]->materialize_node_4_lastlevel(cur_level, tupleQdags.roots[i]);
                }
                else {
                    children &= Q[i]->materialize_node_4(cur_level, tupleQdags.roots[i], rank_vector[i]);
                }
            }
            else if (nAtt == 5) {
                if(cur_level == max_level){
                    children &= Q[i]->materialize_node_5_lastlevel(cur_level, tupleQdags.roots[i]);
                }
                else {
                    children &= Q[i]->materialize_node_5(cur_level, tupleQdags.roots[i], rank_vector[i]);
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
            uint64_t* root_temp= new uint64_t[nQ];
            uint16_t* coordinatesTemp = new uint16_t[l];
            child = children_to_recurse[i];

            for(uint16_t k = 0; k < l; k++)
                coordinatesTemp[k] = tupleQdags.coordinates[k];
            transformCoordinates(coordinatesTemp, l, diff_level, child);

            // compute the coordinates if it's a leaf
            if(cur_level == max_level){
                delete[] root_temp;
                results_points.push_back(coordinatesTemp);
                if(bounded_result && ++results >= UPPER_BOUND)
                    return true;
            } else{
                // compute the weight of the tuple (ONLY if it's not a leaf)
                double total_weight = DBL_MAX;
                // calculate the weight of the tuple
                for (uint64_t j = 0; j < nQ; j++) {
                    // we store the parent quadtree_formula that corresponds in the original quadtree of each qdag
                    root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);
                    uint64_t n_leaves_child_node = Q[j]->get_num_leaves(cur_level, Q[j]->getM(child));
                    //cout << "n_leaves_child_node: " << n_leaves_child_node << endl;
                    if (n_leaves_child_node < total_weight) {
                        total_weight = n_leaves_child_node;
                    }
                }
                if(type_order_fun == TYPE_FUN_DENSITY_LEAVES) // density estimator, otherwise it's the number of leaves (min of the tuple)
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
 * @param grid_size the size of the grid. Needed when type_order_fun is 1.
 * @param type_order_fun 0 num leaves, 1 density. (n1+n2+...+nn) or (n1/100+n2/100+...+nn/100)
 * @param results_points vector with the coordinates of the output.
 * @return
 */
bool multiJoinPartialResults(
        vector<qdag> &Q,
        bool bounded_result,
        uint64_t UPPER_BOUND,
        uint64_t grid_size,
        uint8_t type_order_fun,
        vector<uint16_t *> &results_points) {

    qdag::att_set A;
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
    uint16_t coordinates[A.size()];
    for(uint16_t i = 0; i < A.size(); i++)
        coordinates[i] = 0;

    priority_queue<qdagWeight> pq; // maxHeap
    pq.push({0, Q_roots, 1, coordinates}); // insert the root of the qda
    AND_partial(Q_star, Q.size(), A.size(),
                bounded_result, UPPER_BOUND,
                max_level, grid_size, type_order_fun,
                pq, results_points);

    cout << "number of results: " << results_points.size() << endl;
//    uint64_t i=0;
//    while(i < results_points.size()){
//        uint16_t* coordinates = results_points[i];
//        for(uint64_t j = 0; j< A.size(); j++){
//            cout << coordinates[j] << " ";
//        }
//        cout << endl;
//        i++;
//    }

    for (uint64_t i = 0; i < Q.size(); i++)
        delete Q_star[i];

    return true;
}


/**
 * AND algorithm that will return the top k results.
 * We only visit the first k results.
 * @param Q
 * @param nQ number of Qdags.
 * @param nAtt number of attributes of the final quadtree.
 * @param roots
 * @param cur_level
 * @param max_level maximum level of the qdag.
 * @param grid_size
 * @param type_order_fun 0 sum , 1 maximum. Only taking into account if partial_results is true.
 * @param coordinates
 * @param size_queue the size of the priority queue.
 * @param top_results a vector with the points coordinates of the output.
 * @return
 */
bool
AND_partial_backtracking(
        qdag *Q[],
        uint16_t nQ,
        uint64_t nAtt,
        uint64_t *roots,
        uint16_t cur_level,
        uint16_t max_level,
        uint64_t grid_size,
        uint8_t type_order_fun,
        uint16_t *coordinates,
        uint64_t size_queue,
        vector<uint16_t*> &top_results) {

    if(top_results.size() >= size_queue){
//        cout << "number of results up: " << top_results.size() << endl;
        return true;
    }
    // number of children of the qdag (extended)
    uint64_t p = Q[0]->nChildren();
    bool just_zeroes = true;
    uint64_t k_d[nQ];
    uint16_t children_to_recurse[p];
    uint64_t i;
    uint64_t children_to_recurse_size = 0;
    uint32_t children = 0xffffffff; // each bit represent a quadtree_formula (empty or not)
    uint16_t l = (uint16_t) log2(p); // bits number to define the quadtree_formula's children

    // last level --> add result to the priority queue
    if (cur_level == max_level){

        for (i = 0; i < nQ && children; ++i){
            //k_d[i] = Q[i]->getKD();
            if (nAtt == 3)
                children &= Q[i]->materialize_node_3_lastlevel(cur_level, roots[i]);
            else if (nAtt == 4)
                children &= Q[i]->materialize_node_4_lastlevel(cur_level, roots[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5_lastlevel(cur_level, roots[i]);
        }

        // number of results
        children_to_recurse_size = bits::cnt((uint64_t) children);
        i = 0;
        uint64_t msb;

        while (i < children_to_recurse_size) {
            msb = __builtin_clz(children);
            children_to_recurse[i] = msb; // write the position of the 1s
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1));
        }

        uint16_t child;
        uint16_t diff_level = max_level-cur_level;
        // we do not call recursively the function AND as we do in the other levels
        // add output to the priority queue of results
        for (i = 0; i < children_to_recurse_size; ++i){
            child = children_to_recurse[i];
            uint16_t* coordinatesTemp = new uint16_t[l];

            for(uint16_t k = 0; k < l; k++)
                coordinatesTemp[k] = coordinates[k];
            transformCoordinates(coordinatesTemp, l, diff_level, child);

            top_results.push_back(coordinatesTemp);
            just_zeroes = false;
            // queue full
            if(top_results.size() >= size_queue){
                break;
            }
        }
    }
    // call recursively in DFS order
    else {
        uint64_t rank_vector[nQ][64];

        for (i = 0; i < nQ && children; ++i){
            k_d[i] = Q[i]->getKD(); // k^d del i-esimo quadtree original
            if (nAtt == 3) {
                children &= Q[i]->materialize_node_3(cur_level, roots[i],rank_vector[i]); // entero de 32 bits, y se hace &,
            }else if (nAtt == 4)
                children &= Q[i]->materialize_node_4(cur_level, roots[i], rank_vector[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5(cur_level, roots[i], rank_vector[i]);
        }
        // number of 1s in the children
        children_to_recurse_size = bits::cnt((uint64_t) children);
        i = 0;
        uint64_t msb;

        while (i < children_to_recurse_size) {
            msb = __builtin_clz(children); // the most significant bit of children
            children_to_recurse[i] = msb; // we store the position of the msb in children_to_recurse
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1)); // delete the msb of children
        }

        uint16_t child;

        uint16_t diff_level = max_level-cur_level;
        uint16_t next_level = cur_level + 1;

        priority_queue<orderJoinQdag> order_to_traverse;
        // veo cada hijo en común que tiene el nodo actual
        uint64_t root_temp[children_to_recurse_size][nQ];
        for (i = 0; i < children_to_recurse_size; ++i) {
            child = children_to_recurse[i]; // the position of the 1s in children
            uint16_t* coordinatesTemp = new uint16_t[l];
            for(uint16_t k = 0; k < l; k++)
                coordinatesTemp[k] = coordinates[k];
            transformCoordinates(coordinatesTemp, l, diff_level, child);

            // compute the weight of the tuple
            double total_weight = DBL_MAX;
            for (uint64_t j = 0; j < nQ; j++) {
                root_temp[i][j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);
                uint64_t n_leaves_child_node = Q[j]->get_num_leaves(cur_level,Q[j]->getM(child));
                if (n_leaves_child_node < total_weight) {
                    total_weight = n_leaves_child_node;
                }
            }

            if(type_order_fun == TYPE_FUN_DENSITY_LEAVES) // density estimator, otherwise it's the number of leaves (min of the tuple)
                total_weight /= grid_size;

            orderJoinQdag this_node = {i, coordinatesTemp, total_weight} ;
            order_to_traverse.push(this_node); // add the tuple to the queue

        }
        // recursive call in order according to the n_leaves
        //while(!order_to_traverse.empty()){
        for (i = 0; i < children_to_recurse_size; ++i) {
            orderJoinQdag order = order_to_traverse.top();
            order_to_traverse.pop();
            AND_partial_backtracking(Q, nQ, nAtt, root_temp[order.index], next_level, max_level, grid_size,
                                     type_order_fun,
                                     order.coordinates, size_queue, top_results);
        }

    }

    return !just_zeroes;
}


/**
 *
 * @param Q
 * @param grid_size
 * @param type_order_fun
 * @param size_queue the size of the priority queue. In this case is a vector.
 * @param top_results the vector to store the top results (output).
 * @return
 */
bool
multiJoinPartialResultsBacktracking(
        vector<qdag> &Q,
        uint64_t grid_size,
        uint8_t type_order_fun,
        int64_t size_queue,
        vector<uint16_t *> &top_results) {
    qdag::att_set A;
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

    uint16_t coordinates[A.size()];
    for(uint16_t i = 0; i < A.size(); i++)
        coordinates[i] = 0;

    AND_partial_backtracking(Q_star, Q.size(), A.size(), Q_roots, 0, max_level, grid_size, type_order_fun,
                             coordinates, size_queue,
                             top_results);

//    cout << "number of results: " << top_results.size() << endl;
//    uint64_t i=0;
//    while(i < top_results.size()){
//        uint16_t* coordinates = top_results[i];
//        for(uint64_t j = 0; j< A.size(); j++){
//            cout << coordinates[j] << " ";
//        }
//        cout << endl;
//        i++;
//    }
    for (uint64_t i = 0; i < Q.size(); i++)
        delete Q_star[i];

    return true;
}


