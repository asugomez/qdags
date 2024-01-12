#include <algorithm>
#include <queue>
#include "qdags.hpp"
#include "parallel_for.hpp"
#include "priority.cpp"


/**
 * Compute the priority
 * @param pri a list of priorities
 * @param type the function type which we are going to compute the priorities.
 * @return the output of f(pri[0], pri[1], ....)
 */
uint64_t computePriority(uint64_t* pri, uint64_t type){
    uint64_t total_pri = 0;
    // TODO: change it por alguna funcion de c que lo haga al tiro.
    for(uint64_t i = 0; i < sizeof(pri); i++){
        if(type == 0){ // sum
            total_pri += pri[i];
        }
        if(type == 1){ // max
            if(pri[i]> total_pri){
                total_pri = pri[i];
            }
        }

    }
    return total_pri;
}


//TODO: hacerla recursiva o no?
/**
 *
 * @param Q
 * @param nQ number of Qdags
 * @param max_level
 * @param pq
 * @param type_priority_fun 0 sum , 1 maximum
 * @param partial_results if we are computing partial or ranked results.
 * @return
 */
bool AND_ordered(qdag *Q[], uint64_t *roots, uint16_t nQ,
                 uint16_t cur_level, uint16_t max_level,
                 vector<uint64_t> bv[], uint64_t last_pos[], uint64_t nAtt,
                 bool bounded_result, uint64_t UPPER_BOUND,
                 priority_queue<qdagPri>& pq, uint16_t type_priority_fun, bool partial_results) {
    uint64_t p = Q[0]->nChildren(); // aridad del quadtree
    if(pq.empty() && cur_level == max_level){
        return true;
    }
    if(pq.empty() && cur_level == 0){ // cuando se está al inicio, meto la tupla de raices
        qdagPri tupleQdag = {cur_level, 0, 1};
        pq.push(tupleQdag);
    }
    while(!pq.empty()){
        qdagPri tupleQdags = pq.top(); //  level, el nodo (indice por el cual partir el join), su prioridad
        pq.pop();
        cur_level = tupleQdags.level;
        // if level == max level --> computer join (copiar codigo de AND)
        bool just_zeroes = true;
        uint64_t k_d[16 /*nQ*/]; //CUIDADO, solo hasta 16 relaciones por query (solo hasta 16 tablas)

        uint16_t children_to_recurse[512 /*p*/]; // CUIDADO, solo hasta 9 atributos distintos por query

        uint64_t i;
        uint64_t children_to_recurse_size = 0;

        uint32_t children = 0xffffffff;
        if(cur_level == max_level) {
            for (i = 0; i < nQ && children; ++i) // materializar nodos
            {
                //k_d[i] = Q[i]->getKD();
                if (nAtt == 3)
                    children &= Q[i]->materialize_node_3_lastlevel(cur_level, roots[i]);
                else if (nAtt == 4)
                    children &= Q[i]->materialize_node_4_lastlevel(cur_level, roots[i]);
                else if (nAtt == 5)
                    children &= Q[i]->materialize_node_5_lastlevel(cur_level, roots[i]);
            }

            children_to_recurse_size = bits::cnt((uint64_t) children); // contar nro de reusltados
            i = 0;
            uint64_t msb;

            while (/*children &&*/ i < children_to_recurse_size) {
                msb = __builtin_clz(children);
                children_to_recurse[i] = msb;
                ++i;
                children &= (((uint32_t) 0xffffffff) >> (msb + 1));
            }

            int64_t last_child = -1;
            uint16_t child;

            for (i = 0; i <
                        children_to_recurse_size; ++i) // aqui NO se hace invocacion recursiva! a diferencia de cuando no se esta en  el ultimo nivel
            {
                child = children_to_recurse[i];

                if (child - last_child > 1)
                    last_pos[cur_level] += (child - last_child - 1);

                last_child = child;
                if (bounded_result && bv[max_level].size() >= UPPER_BOUND)
                    return false;
                else {
                    bv[cur_level].push_back(last_pos[cur_level]++);
                    just_zeroes = false;
                }
            }

            if (p - last_child > 1)
                last_pos[cur_level] += (p - last_child - 1);
        }
        else {
            // TODO: sumar las prioridades
            cur_level += 1;
            for(uint64_t i = 0 ; i < p; i++){ // recorriendo las nuevas p posibles tuplas
                bool insertTuple = true;
                uint64_t priorities[p];
                for(uint16_t j = 0; j < nQ; j++){ // recorrer la tupla de qdas
                    // TOOD: que pasa cuando accedo a un nodo que no exisste: ejemplo 3 hijo de R[level] y solo tiene 2
                    uint64_t start_pos_next_level = Q[j]->getBv()[cur_level].rank(i*Q[j]->getKD());//bv[level].rank(node*k_d);
                    uint8_t existNode = Q[j]->getBv()[cur_level].get_4_bits(start_pos_next_level);
                    if(existNode){
                        if(partial_results){
                            priorities[j] = Q[j]->get_num_leaves(cur_level,i);
                        } else{ // ranked_results
                            //TODO: need a rMq
                        }
                    } else{ // si un nodo de la tupla es vacío, no debo insertar esos hijos
                        insertTuple = false;
                        break;
                    }
                }
                if(insertTuple){
                    uint64_t pri = computePriority(priorities, type_priority_fun);
                    pq.push({cur_level, i, pri});
                }
            }
        }

        // else: insertar las nuevas p tuplas a la pri queue segun su prioridad
    }
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
bool AND(qdag *Q[], uint64_t *roots, uint16_t nQ,
         uint16_t cur_level, uint16_t max_level,
         vector<uint64_t> bv[], uint64_t last_pos[], uint64_t nAtt,
         bool bounded_result, uint64_t UPPER_BOUND) {
    uint64_t p = Q[0]->nChildren(); // aridad del quadtree
    bool result = false;
    //uint64_t root_temp[nQ];
    bool just_zeroes = true;
    uint64_t k_d[16 /*nQ*/]; //CUIDADO, solo hasta 16 relaciones por query (solo hasta 16 tablas)

    uint16_t children_to_recurse[512 /*p*/]; // CUIDADO, solo hasta 9 atributos distintos por query

    uint64_t i;
    uint64_t children_to_recurse_size = 0;

    uint32_t children = 0xffffffff;

    if (cur_level == max_level) // ver si llego al ultimo nivel o no, y finalizo:
    {
        for (i = 0; i < nQ && children; ++i) // materializar nodos
        {
            //k_d[i] = Q[i]->getKD();
            if (nAtt == 3)
                children &= Q[i]->materialize_node_3_lastlevel(cur_level, roots[i]);
            else if (nAtt == 4)
                children &= Q[i]->materialize_node_4_lastlevel(cur_level, roots[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5_lastlevel(cur_level, roots[i]);
        }

        children_to_recurse_size = bits::cnt((uint64_t) children); // contar nro de reusltados
        i = 0;
        uint64_t msb;

        while (/*children &&*/ i < children_to_recurse_size) {
            msb = __builtin_clz(children);
            children_to_recurse[i] = msb;
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1));
        }

        int64_t last_child = -1;
        uint16_t child;

        for (i = 0; i <
                    children_to_recurse_size; ++i) // aqui NO se hace invocacion recursiva! a diferencia de cuando no se esta en  el ultimo nivel
        {
            child = children_to_recurse[i];

            if (child - last_child > 1)
                last_pos[cur_level] += (child - last_child - 1);

            last_child = child;
            if (bounded_result && bv[max_level].size() >= UPPER_BOUND)
                return false;
            else {
                bv[cur_level].push_back(last_pos[cur_level]++);
                just_zeroes = false;
            }
        }

        if (p - last_child > 1)
            last_pos[cur_level] += (p - last_child - 1);
    } else // nivel q no es el ultimo
    {
        uint64_t root_temp[16 /*nQ*/]; // CUIDADO, solo hasta 16 relaciones por query
        uint64_t rank_vector[16][64];

        // TODO: pero esto va qdag por qdag, y no atributo x atributo del qdag
        //mientras hayan qdags, por cada qdag de lquery
        for (i = 0; i < nQ && children; ++i) // iterar, nQ tamaño de la query (cantidad de tablas),
        {
            k_d[i] = Q[i]->getKD(); // k^d del i-esimo qdag
            if (nAtt == 3)
                // le materializo el nodo actual del qdag i.
                // te retorna un entero de 32 bits (materializacion) y se le hace AND con el children --> sobreviven solo los 1s que corresponde al nodo del primer qdag
                children &= Q[i]->materialize_node_3(cur_level, roots[i],
                                                     rank_vector[i]); // entero de 32 bits, y se hace &,
                // sobreviven solo los 1s
            else if (nAtt == 4)
                children &= Q[i]->materialize_node_4(cur_level, roots[i], rank_vector[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5(cur_level, roots[i], rank_vector[i]);
            // despues de hacer eso con todos los qdags, en children quedaran aquellos hijos por los que queremos descender
            // bajar por aquellos hijos q esten presentes en TODOS los qdags!!!!
        }
        // en children estan aquellos 1 por aquellos hijos por donde necesitamos descender --> estamos buscando la interseccion rapida

        children_to_recurse_size = bits::cnt(
                (uint64_t) children); // pop:count : cantidad de 1s q tiene un arreglo de bits
        // entonces nos dira el nro de hijos q tendremos q recorrer! children_to_recurse_size
        i = 0;
        uint64_t msb;

        while (/*children &&*/ i <
                               children_to_recurse_size) {   /// dame el 1 mas significativo dado un entero --> tiempo constante obtenemos cada uno de los 1s
            msb = __builtin_clz(children); // el 1 mas significativo
            children_to_recurse[i] = msb; // aqui lo guardo para decir que tengo q bajar por ese hijo
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1)); // borro ese 1 mas significativo
        }

        int64_t last_child = -1;
        uint16_t child;

        for (i = 0; i < children_to_recurse_size; ++i) {

            child = children_to_recurse[i];

            for (uint64_t j = 0; j < nQ; j++) // la posicion del hijo del nodo actual (donde esta del quadtree)
                root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1); // la raiz de cada uno de los qdags
            // por donde me tengo q mover en cada uno de los qdags

            if (child - last_child > 1) // ademas de recorrer el quadtree, estoy generando el resultado
                last_pos[cur_level] += (child - last_child -
                                        1); // ultima posicion del nivel en el q estoy (para el resultado)

            last_child = child;

            if (bounded_result && bv[max_level].size() >=
                                  UPPER_BOUND) // si ya encontre los resultados que me estaban pidiendo --> termino el computo (mato las ramas)
                return false;
            else if (cur_level == max_level ||
                     AND(Q, root_temp, nQ, cur_level + 1, max_level, bv, last_pos, nAtt, bounded_result,
                         UPPER_BOUND)) // llamada recursiva
            {
                bv[cur_level].push_back(last_pos[cur_level]++); // colocamos un 1 un hijo en el resultado
                just_zeroes = false;
            } else {
                if (cur_level < max_level)
                    last_pos[cur_level + 1] -= p; // me devuelvo lo q habia avanzado
                last_pos[cur_level]++;
            }
        }

        if (p - last_child > 1)
            last_pos[cur_level] += (p - last_child - 1);
    }

    return !just_zeroes;
}

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
    AND(Q_star, Q_roots, Q.size(), 0, Q_star[0]->getHeight() - 1, bv, last_pos, A.size(), bounded_result, UPPER_BOUND);
    // constuyo el qdag con bv
    qdag *qResult = new qdag(bv, A, Q_star[0]->getGridSide(), Q_star[0]->getK(), (uint8_t) A.size());

    return qResult;
}

