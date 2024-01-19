#include <algorithm>
#include <queue>
#include <cstdint>
#include <cmath>
#include "qdags.hpp"
#include "parallel_for.hpp"
#include "qdagWeight.cpp"
#include "utils.cpp"


/**
 * Compute the priority
 * @param pri a list of qdags weight of the same node
 * @param type the function type which we are going to compute the priorities.
 * @return the output of f(pri[0], pri[1], ....) and
 */
uint64_t computeNodeWeight(qdagWeight* qdagsWeights, uint64_t type){
    uint64_t total_weight = 0;
    // TODO: change it por alguna funcion de c que lo haga al tiro.
    for(uint64_t i = 0; i < sizeof(qdagsWeights); i++){
        if(type == 0){ // sum
            total_weight += qdagsWeights[i].weight;
        }
        if(type == 1){ // max
            if(qdagsWeights[i].weight > total_weight){
                total_weight = qdagsWeights[i].weight;
            }
        }

    }
    return total_weight;
}

uint64_t computeWeight(uint64_t* qdagsWeights, uint64_t type){
    uint64_t total_weight = 0;
    // TODO: change it por alguna funcion de c que lo haga al tiro.
    for(uint64_t i = 0; i < sizeof(qdagsWeights); i++){
        if(type == 0){ // sum
            total_weight += qdagsWeights[i];
        }
        if(type == 1){ // max
            if(qdagsWeights[i] > total_weight){
                total_weight = qdagsWeights[i];
            }
        }

    }
    return total_weight;
}


// antes de llamar a AND_ordered
// hacer esto: (poner la primera tupla con las raices
//        qdagWeight tupleQdag = {cur_level, 0, 1};
//        pq.push(tupleQdag);
//    }

//TODO: hacerla recursiva o no?
// TODO declare const int 0,1,2 for the type fun
/**
 *
 * @param Q
 * @param nQ number of Qdags
 * @param max_level
 * @param pq
 * @param type_priority_fun 0 sum , 1 maximum
 * @param type_order_fun 0 num leaves, 1 density
 * @param partial_results if we are computing partial or ranked results.
 * @return
 */
bool AND_ordered(qdag *Q[], uint64_t *roots, uint16_t nQ,
                 uint16_t max_level,
                 vector<uint64_t> bv[], uint64_t last_pos[], uint64_t nAtt,
                 bool bounded_result, uint64_t UPPER_BOUND,
                 priority_queue<qdagWeight>& pq, bool partial_results, uint16_t type_priority_fun, uint16_t type_order_fun) {
    // TODO: see if p is correct or it has to be the k_d
    uint64_t p = Q[0]->nChildren(); // aridad del quadtree
    while(!pq.empty()){
        qdagWeight tupleQdags = pq.top(); //  level, el nodo (indice por el cual partir el join), su prioridad
        pq.pop();
        uint16_t cur_level = tupleQdags.level;
        uint16_t cur_node = tupleQdags.level;
        // if it's a leaf, output the point coordenates
        if(cur_level == max_level){


        } else{
            // calcular la prioridad de cada tupla e insertar la struct a la priority queue
            for(uint64_t i = 0; i < p ; i++){
                bool insert = true;
                uint64_t total_weight = partial_results ? UINT64_MAX : 0;
                // TODO: maybe before doing this, see if there is an empty node (evittar hacer este calculo)
                // TODO: hacer esto mas rapido
                for(uint16_t j = 0; j < nQ; j++){
                    // for partial results, compute the num leaves
                    if(partial_results){
                        uint64_t n_leaves_ith_node = Q[j]->get_num_leaves_ith_node(cur_level, cur_node, i);
                        // empty node
                        if(n_leaves_ith_node == 0){
                            insert = false;
                            break;
                        }
                        if (n_leaves_ith_node < total_weight) {
                            // TODO declare const int 0,1,2 for the type fun
                            if(type_order_fun == 0) // min num leaves estimator
                                total_weight = n_leaves_ith_node;
                            else if(type_order_fun == 1) // density estimator
                                // TODO: replace 100 for the grid size
                                total_weight = n_leaves_ith_node / 100;
                        }
                    } else{ // ranked results
                        // TODO: do the ranked results
                        uint64_t priority_ith_node ; // TODO_ define priority ith node (rMq)
                        if(type_priority_fun == 0) // sum
                            total_weight += priority_ith_node;
                        else if(type_priority_fun == 1) { // max
                            if (total_weight < priority_ith_node) {
                                total_weight = priority_ith_node;
                            }
                        } else {
                            // TODO: another monotone function?
                        }
                    }
                }
                // insert tre struct in the pri queue
                if(insert){
                    uint16_t next_level = cur_level+1;
                    // insert the i-th node of this level
                    vector<bool> nodePath;// = new vector<bool>(ceil(log2(p)));
                    decToBinary(i, nodePath);
                    // TODO: check the result
                    for(uint16_t j = 0; j < nodePath.size(); j++){
                        tupleQdags.bv.push_back(nodePath[j]);
                    }
                    qdagWeight this_node = {next_level, i, total_weight, &bv} ;
                    pq.push(this_node);
                }

            }
        }


    }

    return true;
}
/**
 * operations:
 * first_child
 * last_child
 * count nro 1 entre un rango de bitvector
 * 
 * Louds_tree.hpp --> size_type (number of children)
 *                  --> child(v,i)
*/

// multi join:
// recibe un vector de qdags
// si se quiere una cota (las primeras 1000 tuplas por ejemplo) --> true
// UPPER_BOUND --> 1000
qdag *multiJoinPartialResultsSize(vector<qdag> &Q, bool bounded_result, uint64_t UPPER_BOUND) {
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

    for (uint64_t i = 0; i < Q.size(); i++) // iteras por los qdags de la query y se eextiende cada uno de los qdags
    {
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

    // vector de enteros bv bitvector, que tienee tantas entradas como altura de los qdags: tiene un vector por cada nivel
    vector<uint64_t> bv[Q_star[0]->getHeight()]; // OJO, asume que todos los qdags son de la misma altura
    // arreglo de la ultima posicion en donde escribi en cada nivel. Necesario pq uno va llenando los bitvector por cada nivel,
    // entonces para ver donde tengo q escribir en cada nivel, el siguiente elemnto
    uint64_t last_pos[Q_star[0]->getHeight()];

    for (uint64_t i = 0; i < Q_star[0]->getHeight(); i++)
        last_pos[i] = 0; // al inicio todo esta en 0


    // ---------------  everything is the same up to here --------------- //
    // arreglo de qdags extendidos, arreglo de roots, también se necesita el nivel, la altura, bv es la respuesta, arreglo de ult. posicion, cantidad de atributos, si debe tener cota, nro de la cota
    AND_ordered(Q_star, Q_roots, Q.size(), 0, Q_star[0]->getHeight() - 1, bv, last_pos, A.size(), bounded_result, UPPER_BOUND);
    // constuyo el qdag con bv
    qdag *qResult = new qdag(bv, A, Q_star[0]->getGridSide(), Q_star[0]->getK(), (uint8_t) A.size());

    return qResult;
}

