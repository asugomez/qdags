#include <algorithm>
#include "qdags.hpp"
//#include "parallel_for.hpp"

void ANDCount(qdag *Q[], uint64_t *roots, uint16_t nQ,
              uint16_t cur_level, uint16_t max_level,
              uint64_t &ntuples, uint64_t nAtt) {

    uint64_t p = Q[0]->nChildren();

    uint64_t i;
    uint64_t children_to_recurse_size = 0;

    if (cur_level == max_level) {
        uint32_t children = 0xffffffff;
        for (i = 0; i < nQ && children; ++i) {
            if (nAtt == 3)
                children &= Q[i]->materialize_node_3_lastlevel(cur_level, roots[i]);
            else if (nAtt == 4)
                children &= Q[i]->materialize_node_4_lastlevel(cur_level, roots[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5_lastlevel(cur_level, roots[i]);
        }
        ntuples += bits::cnt((uint64_t) children);
        return;
    } else {
        uint64_t k_d[16 /*nQ*/];       // CUIDADO, solo hasta 16 relaciones por query
        uint64_t root_temp[16 /*nQ*/]; // CUIDADO, solo hasta 16 relaciones por query
        uint64_t rank_vector[16][64];
        uint16_t children_to_recurse[512]; // CUIDADO, solo hasta 9 atributos distintos

        uint32_t children = 0xffffffff;
        for (i = 0; i < nQ && children; ++i) {
            k_d[i] = Q[i]->getKD();
            if (nAtt == 3)
                children &= Q[i]->materialize_node_3(cur_level, roots[i], rank_vector[i]);
            else if (nAtt == 4)
                children &= Q[i]->materialize_node_4(cur_level, roots[i], rank_vector[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5(cur_level, roots[i], rank_vector[i]);
        }

        children_to_recurse_size = bits::cnt((uint64_t) children);
        i = 0;
        uint64_t msb;

        while (/*children &&*/ i < children_to_recurse_size) {
            msb = __builtin_clz(children);
            children_to_recurse[i] = msb;
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1));
        }

        uint16_t child;

        for (i = 0; i < children_to_recurse_size; ++i) {
            child = children_to_recurse[i];

            for (uint64_t j = 0; j < nQ; j++)
                root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);

            ANDCount(Q, root_temp, nQ, cur_level + 1, max_level, ntuples, nAtt);
        }
    }
}

void parANDCount(uint16_t totalThreads, uint16_t threadId, uint16_t levelOfCut,
                 qdag *Q[], uint64_t *roots, uint16_t nQ,
                 uint16_t cur_level, uint16_t max_level,
                 uint64_t &ntuples, uint64_t nAtt, uint64_t ancestor[]) {

    uint64_t p = Q[0]->nChildren();

    uint64_t i;
    uint64_t children_to_recurse_size = 0;

    if (cur_level == max_level) {
        uint32_t children = 0xffffffff;
        for (i = 0; i < nQ && children; ++i) {
            if (nAtt == 3)
                children &= Q[i]->materialize_node_3_lastlevel(cur_level, roots[i]);
            else if (nAtt == 4)
                children &= Q[i]->materialize_node_4_lastlevel(cur_level, roots[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5_lastlevel(cur_level, roots[i]);
        }
        ntuples += bits::cnt((uint64_t) children);
        return;
    } else {
        uint64_t k_d[16 /*nQ*/];       // CUIDADO, solo hasta 16 relaciones por query
        uint64_t root_temp[16 /*nQ*/]; // CUIDADO, solo hasta 16 relaciones por query
        uint64_t rank_vector[16][64];
        uint16_t children_to_recurse[512]; // CUIDADO, solo hasta 9 atributos distintos

        uint32_t children = 0xffffffff;
        for (i = 0; i < nQ && children; ++i) {
            k_d[i] = Q[i]->getKD();
            if (nAtt == 3)
                children &= Q[i]->materialize_node_3(cur_level, roots[i], rank_vector[i]);
            else if (nAtt == 4)
                children &= Q[i]->materialize_node_4(cur_level, roots[i], rank_vector[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5(cur_level, roots[i], rank_vector[i]);
        }

        children_to_recurse_size = bits::cnt((uint64_t) children);
        i = 0;
        uint64_t msb;
        uint64_t firstNodeId;

        if (cur_level == levelOfCut) {
            firstNodeId = 0;
            uint64_t nodesPerLevel = 1 << nAtt;
            uint64_t multiplier = nodesPerLevel;

            for (int level = cur_level - 1; level >= 0; level--) {
                firstNodeId += ancestor[level] * multiplier;
                multiplier *= nodesPerLevel;
            }
        }

        while (/*children &&*/ i < children_to_recurse_size) {
            msb = __builtin_clz(children);
            children_to_recurse[i] = msb;
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1));
        }

        uint16_t child;

        for (i = 0; i < children_to_recurse_size; ++i) {
            child = children_to_recurse[i];
            // cout << get_child_se_quadtree << ' ';

            if (cur_level == levelOfCut) {
                // CHECK WHETHER TO SKIP THIS NODE
                // THIS ASSUMES THE FIRST LEVELS ARE FULL
                if ((firstNodeId + child) % totalThreads != threadId)
                    // This node corresponds to another thread, skip it!
                    continue;
            }


            for (uint64_t j = 0; j < nQ; j++)
                root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);

            if (cur_level <= levelOfCut) {
                ancestor[cur_level] = child;
                parANDCount(totalThreads, threadId, levelOfCut, Q, root_temp, nQ, cur_level + 1, max_level, ntuples,
                            nAtt, ancestor);
            } else
                ANDCount(Q, root_temp, nQ, cur_level + 1, max_level, ntuples, nAtt);
        }
    }
}


/**
 * Travers recursively the qdags and if the intersection is not empty, write the output (leaf) on path
 * @param Q
 * @param roots
 * @param nQ
 * @param cur_level
 * @param max_level
 * @param bv int vector to store the position of 1 on each level of the output bitvector.
 * @param last_pos
 * @param nAtt
 * @param bounded_result
 * @param UPPER_BOUND
 * @return
 */
bool AND(qdag *Q[], uint64_t *roots, uint16_t nQ,
         uint16_t cur_level, uint16_t max_level,
         vector<uint64_t> bv[], uint64_t last_pos[], uint64_t nAtt,
         bool bounded_result, uint64_t UPPER_BOUND) {
    // number of children of the qdag (extended)
    uint64_t p = Q[0]->nChildren();
    bool result = false;
    //uint64_t root_temp[nQ];
    bool just_zeroes = true;
    uint64_t k_d[16 /*nQ*/]; //CUIDADO, solo hasta 16 relaciones por query (solo hasta 16 tablas)

    uint16_t children_to_recurse[512 /*p*/]; // CUIDADO, solo hasta 9 atributos distintos por query

    uint64_t i;
    uint64_t children_to_recurse_size = 0;

    uint32_t children = 0xffffffff; // each bit represent a node (empty or not)

    // in the last level we materialize the node
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

        while (/*children &&*/ i < children_to_recurse_size) {
            msb = __builtin_clz(children);
            children_to_recurse[i] = msb; // write the position of the 1s
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1));
        }

        int64_t last_child = -1;
        uint16_t child;

        // we do not call recursively the function AND as we do in the other levels
        for (i = 0; i < children_to_recurse_size; ++i){
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
    } else {
        uint64_t root_temp[16 /*nQ*/]; // CUIDADO, solo hasta 16 relaciones por query
        uint64_t rank_vector[16][64];

        for (i = 0; i < nQ && children; ++i){
            k_d[i] = Q[i]->getKD(); // k^d del i-esimo quadtree original
            if (nAtt == 3) {
                //uint32_t testing_node = Q[i]->materialize_node_3(cur_level, roots[i], rank_vector[i]);
                //std::bitset<32> t_node(testing_node); // DEBUG: ver los bits del materialize node
                // te retorna un entero de 32 bits (materializacion) y se le hace AND con el children --> sobreviven solo los 1s que corresponde al nodo del primer qdag
                children &= Q[i]->materialize_node_3(cur_level, roots[i],rank_vector[i]); // entero de 32 bits, y se hace &,
                // sobreviven solo los 1s
            }else if (nAtt == 4)
                children &= Q[i]->materialize_node_4(cur_level, roots[i], rank_vector[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5(cur_level, roots[i], rank_vector[i]);

            // despues de hacer eso con todos los qdags, en children quedaran aquellos hijos por los que queremos descender
            // bajar por aquellos hijos q esten presentes en TODOS los qdags!!!!
        }
        // en children estan aquellos 1 por aquellos hijos por donde necesitamos descender --> estamos buscando la interseccion rapida
        // number of 1s in the children
        children_to_recurse_size = bits::cnt((uint64_t) children);
        i = 0;
        uint64_t msb;

        while (/*children &&*/ i < children_to_recurse_size) {
            msb = __builtin_clz(children); // the most significant bit of children
            children_to_recurse[i] = msb; // we store the position of the msb in children_to_recurse
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1)); // delete the msb of children
        }

        int64_t last_child = -1;
        uint16_t child;

        for (i = 0; i < children_to_recurse_size; ++i) {
            child = children_to_recurse[i]; // the position of the 1s in children
            for (uint64_t j = 0; j < nQ; j++) {
                uint16_t cur_node_test = Q[j]->getM(child);
                std::bitset<64> ith_node(cur_node_test);
                // the node in the j-th quadtree where we are moving
                // it's like a mapping between the get_child_se_quadtree in the output and which node in the qdag we are moving
                root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);
            }

            // we are writing the result in the output bitvector
            if (child - last_child > 1)
                last_pos[cur_level] += (child - last_child - 1); // last position we wrote in the current level.

            last_child = child;

            if (bounded_result && bv[max_level].size() >= UPPER_BOUND)
                return false;
            else if (cur_level == max_level ||
                     AND(Q, root_temp, nQ, cur_level + 1, max_level, bv, last_pos, nAtt, bounded_result,
                         UPPER_BOUND)) // recursive call
            {
                bv[cur_level].push_back(last_pos[cur_level]++); // we write the current result in the output bitvector
                just_zeroes = false;
            } else {
                if (bounded_result && bv[max_level].size() >= UPPER_BOUND) // finish the build of the quadtree
                    bv[cur_level].push_back(last_pos[cur_level]++); // we write the current result in the output bitvector
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

bool parAND(uint16_t totalThreads, uint16_t threadId, uint16_t levelOfCut, std::mutex &sharedMutex,
            qdag *Q[], uint64_t *roots, uint16_t nQ,
            uint16_t cur_level, uint16_t max_level,
            vector<uint64_t> bv[], uint64_t last_pos[], uint64_t ancestor[], uint64_t nAtt,
            bool bounded_result, uint64_t UPPER_BOUND) {
    uint64_t p = Q[0]->nChildren();
    bool result = false;
    //uint64_t root_temp[nQ];
    bool just_zeroes = true;
    uint64_t k_d[16 /*nQ*/]; //CUIDADO, solo hasta 16 relaciones por query

    uint16_t children_to_recurse[512 /*p*/]; // CUIDADO, solo hasta 9 atributos distintos por query

    uint64_t i;
    uint64_t children_to_recurse_size = 0;

    uint32_t children = 0xffffffff;

    if (cur_level == max_level) {
        for (i = 0; i < nQ && children; ++i) {
            //k_d[i] = Q[i]->getKD();
            if (nAtt == 3)
                children &= Q[i]->materialize_node_3_lastlevel(cur_level, roots[i]);
            else if (nAtt == 4)
                children &= Q[i]->materialize_node_4_lastlevel(cur_level, roots[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5_lastlevel(cur_level, roots[i]);
        }

        children_to_recurse_size = bits::cnt((uint64_t) children);
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

        for (i = 0; i < children_to_recurse_size; ++i) {
            child = children_to_recurse[i];

            if (child - last_child > 1)
                last_pos[cur_level] += (child - last_child - 1);

            last_child = child;
            if (bounded_result && bv[max_level].size() >= UPPER_BOUND)
                return false;
            else {
                sharedMutex.lock();
                bv[cur_level].push_back(last_pos[cur_level]++);
                sharedMutex.unlock();

                just_zeroes = false;
            }
        }

        if (p - last_child > 1)
            last_pos[cur_level] += (p - last_child - 1);
    } else {
        uint64_t root_temp[16 /*nQ*/]; // CUIDADO, solo hasta 16 relaciones por query
        uint64_t rank_vector[16][64];

        for (i = 0; i < nQ && children; ++i) {
            k_d[i] = Q[i]->getKD();
            if (nAtt == 3)
                children &= Q[i]->materialize_node_3(cur_level, roots[i], rank_vector[i]);
            else if (nAtt == 4)
                children &= Q[i]->materialize_node_4(cur_level, roots[i], rank_vector[i]);
            else if (nAtt == 5)
                children &= Q[i]->materialize_node_5(cur_level, roots[i], rank_vector[i]);
        }

        children_to_recurse_size = bits::cnt((uint64_t) children);
        i = 0;
        uint64_t msb;

        uint64_t firstNodeId;

        if (cur_level == levelOfCut) {
            firstNodeId = 0;
            uint64_t nodesPerLevel = 1 << nAtt;
            uint64_t multiplier = nodesPerLevel;

            for (int level = cur_level - 1; level >= 0; level--) {
                firstNodeId += ancestor[level] * multiplier;
                multiplier *= nodesPerLevel;
            }
        }

        while (/*children &&*/ i < children_to_recurse_size) {
            msb = __builtin_clz(children);
            children_to_recurse[i] = msb;
            ++i;
            children &= (((uint32_t) 0xffffffff) >> (msb + 1));
        }

        int64_t last_child = -1;
        uint16_t child;

        for (i = 0; i < children_to_recurse_size; ++i) {

            child = children_to_recurse[i];

            if (cur_level == levelOfCut) {
                // CHECK WHETHER TO SKIP THIS NODE
                // THIS ASSUMES THE FIRST LEVELS ARE FULL
                if ((firstNodeId + child) % totalThreads != threadId)
                    // This node corresponds to another thread, skip it!
                    continue;
            }

            for (uint64_t j = 0; j < nQ; j++)
                root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);

            if (child - last_child > 1)
                last_pos[cur_level] += (child - last_child - 1);

            last_child = child;

            if (bounded_result && bv[max_level].size() >= UPPER_BOUND)
                return false;

            ancestor[cur_level] = child;
            if (cur_level == max_level ||
                //  (cur_level > levelOfCut
                //       ? AND(Q, root_temp, nQ, cur_level + 1, max_level, path, last_pos, nAtt, bounded_result, UPPER_BOUND)
                //       : parAND(totalThreads, threadId, levelOfCut, Q, root_temp, nQ, cur_level + 1, max_level, path, last_pos, nAtt, bounded_result, UPPER_BOUND))
                parAND(totalThreads, threadId, levelOfCut, sharedMutex, Q, root_temp, nQ, cur_level + 1, max_level, bv,
                       last_pos, ancestor, nAtt, bounded_result, UPPER_BOUND)) {
                sharedMutex.lock();
                bv[cur_level].push_back(last_pos[cur_level]++);
                sharedMutex.unlock();

                just_zeroes = false;
            } else {
                if (cur_level < max_level)
                    last_pos[cur_level + 1] -= p;
                last_pos[cur_level]++;
            }
        }

        if (p - last_child > 1)
            last_pos[cur_level] += (p - last_child - 1);
    }

    return !just_zeroes;
}

uint64_t multiJoinCount(vector<qdag> &Q) {
    qdag::att_set A;
    map<uint64_t, uint8_t> attr_map;

    // computes the union of the attribute sets
    for (uint64_t i = 0; i < Q.size(); i++) {
        uint64_t nAttr = Q[i].nAttr();
        for (uint64_t j = 0; j < nAttr; j++)
            attr_map[Q[i].getAttr(j)] = 1;
    }

    for (map<uint64_t, uint8_t>::iterator it = attr_map.begin(); it != attr_map.end(); it++)
        A.push_back(it->first);

    qdag *Q_star[Q.size()];
    uint64_t Q_roots[Q.size()];

    for (uint64_t i = 0; i < Q.size(); i++) {
        Q_star[i] = Q[i].extend(A);
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
        Q_roots[i] = 0; // root of every qdag
    }

    uint64_t ntuples = 0;
    ANDCount(Q_star, Q_roots, Q.size(), 0, Q_star[0]->getHeight() - 1, ntuples, A.size());

    for (uint64_t i = 0; i < Q.size(); i++)
        delete Q_star[i];

    return ntuples;
}

/*
uint64_t parMultiJoinCount(vector<qdag> &Q) {
    qdag::att_set A;
    map<uint64_t, uint8_t> attr_map;

    // computes the union of the attribute sets
    for (uint64_t i = 0; i < Q.size(); i++) {
        uint64_t nAttr = Q[i].nAttr();
        for (uint64_t j = 0; j < nAttr; j++)
            attr_map[Q[i].getAttr(j)] = 1;
    }

    for (map<uint64_t, uint8_t>::iterator it = attr_map.begin(); it != attr_map.end(); it++)
        A.push_back(it->first);

    qdag *Q_star[Q.size()];
    uint64_t Q_roots[Q.size()];

    for (uint64_t i = 0; i < Q.size(); i++) {
        Q_star[i] = Q[i].extend(A);
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
        Q_roots[i] = 0; // root of every qdag
    }

    uint64_t ntuples = 0;

    // -------
    unsigned nb_threads_hint = THREADS_BY_CORE * std::thread::hardware_concurrency();
    unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);
    uint16_t levelOfCut = 5;

    std::mutex tuplesMutex;
    auto height = Q_star[0]->getHeight();

    parallel_for(nb_threads, [&](int start, int end) {
        for (int i = start; i < end; ++i) {
            uint64_t ntuples_i = 0;
            uint64_t ancestor[height];

            parANDCount(nb_threads, i, levelOfCut, Q_star, Q_roots, Q.size(), 0, height - 1, ntuples_i, A.size(),
                        ancestor);

            tuplesMutex.lock();
            ntuples += ntuples_i;
            tuplesMutex.unlock();
        }
    });

    for (uint64_t i = 0; i < Q.size(); i++)
        delete Q_star[i];

    return ntuples;
}
*/

/**
 * Multi join between the qdags in Q.
 * @param Q a vector of qdags. Every qdag has the same height.
 * @param bounded_result if the result is bounded.
 * @param UPPER_BOUND the limit of the result if bounded.
 * @return a quadtree (qdag not extended) with the result of the multi join.
 */
qdag *multiJoin(vector<qdag> &Q, bool bounded_result, uint64_t UPPER_BOUND) {
    qdag::att_set A;
    map<uint64_t, uint8_t> attr_map;

    for (uint64_t i = 0; i < Q.size(); i++) {
        uint64_t nAttr = Q[i].nAttr();
        for (uint64_t j = 0; j < nAttr; j++)
            attr_map[Q[i].getAttr(j)] = 1;
    }

    // set of attributs
    for (map<uint64_t, uint8_t>::iterator it = attr_map.begin(); it != attr_map.end(); it++)
        A.push_back(it->first);

    qdag *Q_star[Q.size()]; // Q star array
    uint64_t Q_roots[Q.size()]; // int array that tells you the current node of each qdag we are.

    // extension of the qdags
    for (uint64_t i = 0; i < Q.size(); i++){
        // preprocessing of the qdags
        // only for binary relations
        Q_star[i] = Q[i].extend(A);
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
        Q_roots[i] = 0; // root of every qdag
    }

    // An int vector where we are going to write the output. It has as many entries as the height of the qdags.
    // On each level, we store the index where we have nonempty nodes.
    vector<uint64_t> bv[Q_star[0]->getHeight()];

    // array of the last position where we wrote on each level. It is useful to write the output on the bitvector abd knowing where to write.
    uint64_t last_pos[Q_star[0]->getHeight()];

    for (uint64_t i = 0; i < Q_star[0]->getHeight(); i++)
        last_pos[i] = 0; // initialize the last position

    // array of extended qdags, array of roots, current level, height, the output bitvector path, arreglo de ult. posicion, number of attributs,
    AND(Q_star, Q_roots, Q.size(), 0, Q_star[0]->getHeight() - 1, bv, last_pos, A.size(), bounded_result, UPPER_BOUND);

    // we create a new qdag with the output store in path (position of 1s on each level)
    qdag *qResult = new qdag(bv, A, Q_star[0]->getGridSide(), Q_star[0]->getK(), (uint8_t) A.size());

//    cout << "positions" << endl;
//    for (uint64_t i = 0; i < Q_star[0]->getHeight(); i++) {
//        for (uint64_t j = 0; j < bv[i].size(); j++) {
//            cout << bv[i][j] << " ";
//        }
//        cout << endl;
//    }

//    for (uint64_t j = 0; j < bv[Q_star[0]->getHeight()-1].size(); j++) {
//        cout << bv[Q_star[0]->getHeight()-1][j] << " ";
//    }
//    cout << endl;

//    cout << " number of results: " << bv[Q_star[0]->getHeight()-1].size() << endl;
//    cout << bv[Q_star[0]->getHeight()-1].size() << endl;

    //cout << endl;

    return qResult;
}

// version paralela
/*
qdag *parMultiJoin(vector<qdag> &Q, bool bounded_result, uint64_t UPPER_BOUND) {
    qdag::att_set A;
    map<uint64_t, uint8_t> attr_map;

    // computes the union of the attribute sets
    for (uint64_t i = 0; i < Q.size(); i++) {
        uint64_t nAttr = Q[i].nAttr();
        for (uint64_t j = 0; j < nAttr; j++)
            attr_map[Q[i].getAttr(j)] = 1;
    }

    for (map<uint64_t, uint8_t>::iterator it = attr_map.begin(); it != attr_map.end(); it++)
        A.push_back(it->first);

    unsigned nb_threads_hint = THREADS_BY_CORE * std::thread::hardware_concurrency();
    unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);
    uint16_t levelOfCut = 5;

    std::mutex tuplesMutex;

    qdag *Q_star[Q.size()];
    uint64_t Q_roots[Q.size()];

    for (uint64_t i = 0; i < Q.size(); i++) {
        Q_star[i] = Q[i].extend(A);
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
        Q_roots[i] = 0; // root of every qdag
    }

    auto height = Q_star[0]->getHeight();
    vector<uint64_t> path[height]; // OJO, asume que todos los qdags son de la misma altura

    parallel_for(nb_threads, [&](int start, int end) {
        for (int threadId = start; threadId < end; ++threadId) {

            uint64_t last_pos[height];
            uint64_t ancestor[height];

            for (uint64_t i = 0; i < height; i++)
                last_pos[i] = 0;

            parAND(nb_threads, threadId, levelOfCut, tuplesMutex, Q_star, Q_roots, Q.size(), 0, height - 1, path,
                   last_pos, ancestor, A.size(), bounded_result, UPPER_BOUND);
        }
    });

    qdag *qResult = new qdag(path, A, Q_star[0]->getGridSide(), Q_star[0]->getK(), (uint8_t) A.size());

    return qResult;
}
*/
