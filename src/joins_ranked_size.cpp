#include <algorithm>
#include "qdags.hpp"

// antes de llamar a AND_ordered
// hacer esto: (poner la primera tupla con las raices
//        qdagWeight tupleQdag = {cur_level, 0, 1};
//        pq.push(tupleQdag);
//    }

struct qdagWeight {
    int16_t level; // the level of the node
    uint64_t* roots; // the parent of the subtree of each qdag that conform the tuple.
    uint64_t weight; // priority or number of leaves of the tuple
    uint64_t path;  // the bits that encode the path down the leaf in the first qdag.
    bool operator<(const qdagWeight& qd) const {
        return weight < qd.weight;
    }
};


void getCoordinates(uint64_t path, uint64_t l, uint64_t height, uint64_t *pointCoord, uint64_t nAttr) {

    uint64_t diff_level;
    for(uint64_t level=0; level<height; level++) {
        diff_level = height - level;
        uint32_t node_cur_level = (uint32_t) path >> (diff_level * l);
        uint64_t n_ones = bits::cnt((uint64_t) node_cur_level);

        uint64_t msb, num_point;
        for (uint64_t k = 0; k < n_ones; k++) {
            msb = __builtin_ctz(node_cur_level);
            node_cur_level &= ~(1 << msb); // delete the msb of children
            msb %= l;
            num_point = (l - 1) - msb; // the index is the inverse
            pointCoord[num_point] += 1 << diff_level;
        }
    }

}



//TODO: hacerla recursiva o no?
// TODO declare const int 0,1,2 for the type fun
/**
 *
 * @param Q
 * @param nQ number of Qdags
 * @param max_level
 * @param pq
 * @param partial_results if we are computing partial or ranked results.
 * @param type_priority_fun 0 sum , 1 maximum
 * @param type_order_fun 0 num leaves, 1 density
 * @return
 */
bool AND_ordered(qdag *Q[], uint16_t nQ,
                 uint64_t max_level, uint64_t nAtt,
                 bool bounded_result, uint64_t UPPER_BOUND, uint64_t grid_size,
                 priority_queue<qdagWeight>& pq, bool partial_results, uint16_t type_priority_fun, uint16_t type_order_fun) {
    // TODO: see if p is correct or it has to be the k_d
    uint64_t p = Q[0]->nChildren(); // number of children of the qdag extended
    uint64_t k_d[nQ];
    uint16_t children_to_recurse[p];
    uint64_t children_to_recurse_size = 0;
    uint32_t children = 0xffffffff;
    uint64_t results = 0;
    while(!pq.empty()){
        qdagWeight tupleQdags = pq.top(); //  level, el nodo (indice por el cual partir el join), su prioridad
        pq.pop();
        int16_t cur_level = tupleQdags.level;
        uint64_t *roots = tupleQdags.roots; // TODO: see if this works
        uint64_t l = (uint64_t) log2(p); // bits number to define the node's children
        // if it's a leaf, output the point coordenates
        uint64_t root_temp[nQ];
        uint64_t rank_vector[nQ][64];
        for (uint64_t i = 0; i < nQ && children; ++i){
            k_d[i] = Q[i]->getKD(); // k^d del i-esimo qdag
            if (nAtt == 3) {
                if (cur_level == max_level) {
                    children &= Q[i]->materialize_node_3_lastlevel(cur_level, roots[i]);
                } else {
                    children &= Q[i]->materialize_node_3(cur_level, roots[i], rank_vector[i]);
                }
                children &= Q[i]->materialize_node_3(cur_level, roots[i], rank_vector[i]);
            } else if (nAtt == 4) {
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

        for (i = 0; i < children_to_recurse_size; ++i) {
            child = children_to_recurse[i];

            uint64_t total_weight = partial_results ? UINT64_MAX : 0;
            if(cur_level != max_level) { // compute the weight of the tuple
                for (uint64_t j = 0; j < nQ; j++) {
                    // we store the parent node that corresponds in the original quadtree of each qdag

                    root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);
                    if (partial_results) {
                        uint64_t n_leaves_child_node = Q[j]->get_num_leaves(cur_level, Q[j]->getM(child));
                        if (n_leaves_child_node < total_weight) {
                            if (type_order_fun == 0) // min num leaves estimator among the qdags tuple
                                total_weight = n_leaves_child_node;
                            else if (type_order_fun == 1) // density estimator
                                total_weight = n_leaves_child_node / grid_size;
                        }
                    } else { // ranked results
                        // TODO: do the ranked results
                        uint64_t priority_ith_node; // TODO_ define priority ith node (rMq)
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
            // --> add the child to the path
            uint16_t cur_child_qdag = Q[0]->getM(child);
            uint16_t diff_level = max_level-cur_level;
            uint64_t path = cur_child_qdag << (diff_level * (uint64_t) log2(k_d[0])); // height * bits to represent the children
            tupleQdags.path += path; // add the bits to the bitvector
            // compute the coordinates
            if(cur_level == max_level){
                // TODO: output coordinates
                // TODO: las coordenadas de un qdag seran las de los otros iguales?
                // DEBUG: see the tupleQdags and path
                // TOODO: esta operacion es lenta, podríamos tener una tabla guardada
                uint64_t coordinates[nAtt];
                for(uint64_t i = 0; i < nAtt; i++){
                    coordinates[i] = 0;
                }
                getCoordinates(tupleQdags.path, l, max_level, coordinates, nAtt);
                cout << "point output: " << tupleQdags.path << endl;
                // TODO: maybe we have to materialize the last level
                if(bounded_result && ++results >= UPPER_BOUND)
                    return true;
            }
            // else: insert the tuple
            cur_level++;
            qdagWeight this_node = {cur_level, root_temp, total_weight, tupleQdags.path} ;
            pq.push(this_node); // add the tuple to the queue
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
                                  bool partial_results, uint16_t type_priority_fun, uint16_t type_order_fun) {
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

    // vector de enteros path bitvector, que tienee tantas entradas como altura de los qdags: tiene un vector por cada nivel
    vector<uint64_t> bv[Q_star[0]->getHeight()]; // OJO, asume que todos los qdags son de la misma altura
    // arreglo de la ultima posicion en donde escribi en cada nivel. Necesario pq uno va llenando los bitvector por cada nivel,
    // entonces para ver donde tengo q escribir en cada nivel, el siguiente elemnto
    uint64_t last_pos[Q_star[0]->getHeight()];

    for (uint64_t i = 0; i < Q_star[0]->getHeight(); i++)
        last_pos[i] = 0; // al inicio todo esta en 0


    // ---------------  everything is the same up to here --------------- //
    // arreglo de qdags extendidos, arreglo de roots, también se necesita el nivel, la altura, path es la respuesta, arreglo de ult. posicion, cantidad de atributos, si debe tener cota, nro de la cota
    uint64_t roots[Q.size()];
    for (uint64_t i = 0; i < Q.size(); i++){
        roots[i] = 0;
    }
    uint64_t grid_size = 100;
    uint64_t max_level = Q_star[0]->getHeight() - 1;
    priority_queue<qdagWeight> pq;
    pq.push({0, roots, 1, 0}); // insert the root of the qdag

    AND_ordered(Q_star, Q.size(),
                max_level, A.size(),
                bounded_result, UPPER_BOUND, grid_size,
                pq, partial_results, type_priority_fun,  type_order_fun);

    return true;
}

