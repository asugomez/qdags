#include <algorithm>
#include "qdags.hpp"

// antes de llamar a AND_ordered
// hacer esto: (poner la primera tupla con las raices
//        qdagWeight tupleQdag = {cur_level, 0, 1};
//        pq.push(tupleQdag);
//    }

const int MAX_ROOTS = 100; // Máximo número de elementos en roots

struct qdagWeight {
    int16_t level; // the level of the node
    uint64_t* roots; // the parent of the subtree of each qdag that conform the tuple.
    double weight; // priority or number of leaves of the tuple
    uint64_t path;  // the bits that encode the path down the leaf in the first qdag.
    bool operator<(const qdagWeight& qd) const {
        return weight < qd.weight;
    }
};

/**
 * Get the coordinates of a point in the grid given the path from the root to the leaf in the quadtree.
 * This function is O(m), with m the number of 1s in the path.
 * @param path see the paper for the definition (page 5).
 * @param l number of bits to define a child on each level
 * @param height height of the quadtree
 * @param pointCoord the coordinates of the point in the grid. We are going to write on this array.
 */
void getCoordinates(uint64_t path, uint16_t l, uint64_t height, uint32_t *pointCoord) {

    uint32_t bits_path = (height + 1) * l;
    path = (uint32_t) path << (31 - bits_path + 1); // move the beginning of the path to the begining of the 32 bits
    uint32_t n_ones = bits::cnt(path);
    uint32_t msb, level;
    // we only pass across the 1s to update the coordinates.
    for (uint32_t k = 0; k < n_ones; k++) {
        msb = __builtin_clz(path);
        level = msb/l;
        path &= (((uint32_t) 0xffffffff) >> (msb + 1)); // delete the msb of children
        uint32_t num_point = msb%l;
        pointCoord[num_point] += 1 << (height - level);
    }

}


// TODO declare const int 0,1,2 for the type fun
/**
 *
 * @param Q
 * @param nQ number of Qdags
 * @param max_level
 * @param pq
 * @param partial_results if we are computing partial or ranked results.
 * @param type_priority_fun 0 sum , 1 maximum. Only taking into account if partial_results is true.
 * @param type_order_fun 0 num leaves, 1 density. Only taking into account if partial_results is false.
 * @return
 */
bool AND_ordered(qdag *Q[], uint16_t nQ,
                 uint64_t max_level, uint64_t nAtt,
                 bool bounded_result, uint64_t UPPER_BOUND,
                 priority_queue<qdagWeight>& pq,
                 bool partial_results, uint8_t type_priority_fun, uint8_t type_order_fun, uint64_t grid_size) {
    uint64_t p = Q[0]->nChildren(); // number of children of the qdag extended
    uint64_t k_d[nQ];
    uint32_t children;
    uint16_t children_to_recurse[p];
    uint64_t children_to_recurse_size;
    int16_t cur_level;
    uint16_t l = (uint64_t) log2(p); // bits number to define the node's children
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
        // DEBUG: see if this works roots tiene como un arreglo
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
        int16_t next_level = cur_level + 1;
        uint64_t path;

        for (i = 0; i < children_to_recurse_size; ++i) {
            uint64_t* root_temp= new uint64_t[nQ];
            child = children_to_recurse[i];
            path = 0;

            // compute the weight of the tuple (ONLY if it's not a leaf)
            double total_weight = partial_results ? DBL_MAX : 0;
            if(cur_level != max_level) {
                // calculate the weight of the tuple
                // DEBUG: ver como funciona esto si es para cad atupla correctamente
                for (uint64_t j = 0; j < nQ; j++) {
                    // we store the parent node that corresponds in the original quadtree of each qdag
                    uint64_t tt = Q[j]->getM(child);
                    uint64_t jiji = rank_vector[j][Q[j]->getM(child)];
                    // HERE we have a problem
                    root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);
                    if (partial_results) {
                        uint64_t n_leaves_child_node = Q[j]->get_num_leaves(cur_level, Q[j]->getM(child));
                        if (n_leaves_child_node < total_weight) {
                            total_weight = n_leaves_child_node;
                        }
                    }
                    else { // ranked results
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
            }
            if(partial_results && type_order_fun == 1) // density estimator, otherwise it's the number of leaves (min of the tuple)
                    total_weight /= grid_size;
            // --> add the child to the path
            //uint16_t cur_child_qdag = Q[0]->getM(child);
            path = child << (diff_level * l); // height * bits to represent the children
            path += tupleQdags.path; // add the bits to the bitvector
            //tupleQdags.path += path; // add the bits to the bitvector
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
                if(bounded_result && ++results >= UPPER_BOUND)
                    return true;
            }
            else{ // insert the tuple
                /*cout << "level " << cur_level << endl;
                for(uint64_t j = 0; j < nQ; j++){
                    cout << "curr root " << root_temp[j] << endl;
                }*/
                cout << endl;
                qdagWeight this_node = {next_level, root_temp, total_weight, path} ;
                pq.push(this_node); // add the tuple to the queue
                // DEBUG: see if the tuple is added to the queue
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
 * @param partial_results whether we are computing partial or ranked results.
 * @param type_priority_fun 0 sum , 1 maximum. (p1+p2+...+pn) or max(p1,p2,...,pn)
 * @param type_order_fun 0 num leaves, 1 density. (n1+n2+...+nn) or (n1/100+n2/100+...+nn/100)
 * @return
 */
bool multiJoinPartialResultsSize(vector<qdag> &Q, bool bounded_result, uint64_t UPPER_BOUND,
                                  bool partial_results, uint8_t type_priority_fun, uint8_t type_order_fun, uint64_t grid_size) {
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
    uint64_t Q_roots[MAX_ROOTS]; // arreglo de enteros: mantiene el nodo actual en el q estas en cada uno de los qdags (el arreglo te dice en donde estas)

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
    priority_queue<qdagWeight> pq;
    pq.push({0, Q_roots, 1, 0}); // insert the root of the qdag

    AND_ordered(Q_star, Q.size(),
                max_level, A.size(),
                bounded_result, UPPER_BOUND,
                pq, partial_results, type_priority_fun,  type_order_fun, grid_size);

    return true;
}

