#include <algorithm>
#include "../qdags.hpp"
//#include <sdsl/rmq_support_sparse_table.hpp>
#include <sdsl/rmq_support.hpp>
#include "../utils.hpp"

using namespace sdsl;
using namespace std;

const uint8_t TYPE_FUN_PRI_SUM = 0;
const uint8_t TYPE_FUN_PRI_MAX = 1;

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
 * @param rMq the rmq structure to get the minimum priority in a range.
 * @param results_points the output.
 * @return true if we accomplished succesfully the join. Otherwise, return false.
 */
bool AND_ranked(
        qdag *Q[],
        uint16_t nQ,
        uint64_t max_level,
        uint64_t nAtt,
        bool bounded_result,
        uint64_t UPPER_BOUND,
        priority_queue<qdagWeight> &pq,
        uint8_t type_priority_fun,
        vector<int_vector<>> &priorities,
        vector<rmq_succinct_sct<false>> &rMq,
        vector<uint256_t>& results_points ) {

    uint64_t p = Q[0]->nChildren(); // number of children of the qdag extended
    uint64_t k_d[nQ];
    uint32_t children;
    uint16_t children_to_recurse[p];
    uint64_t children_to_recurse_size;
	uint64_t i;
    uint16_t cur_level;
    uint64_t results = 0;
    while(!pq.empty()){
        children = 0xffffffff;
        qdagWeight tupleQdags = pq.top();
        pq.pop();
        cur_level = tupleQdags.level;
        // if it's a leaf, output the point coordenates
        if(cur_level > max_level){
            results_points.push_back(tupleQdags.path);
            if(bounded_result && ++results >= UPPER_BOUND)
                return true;
        }
        else{
            uint64_t rank_vector[16 /*nQ*/][64];
            for (i = 0; i < nQ && children; ++i){
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
			if(children_to_recurse_size == 0){
				delete[] tupleQdags.roots;
			}
			else{
				i = 0;
				uint64_t msb; // most significant bit
				while (i < children_to_recurse_size) {
					msb = __builtin_clz(children);
					children_to_recurse[i] = msb;
					++i;
					children &= (((uint32_t) 0xffffffff) >> (msb + 1));
				}

				uint16_t child;
				uint16_t next_level = cur_level + 1;

				double total_weight = 0;
				// compute the coordinates if it's a leaf
				if(cur_level == max_level){
					for (i = 0; i < children_to_recurse_size; ++i) {
						child = children_to_recurse[i];
						uint256_t newPath = tupleQdags.path;
						getNewMortonCodePath(newPath, nAtt, cur_level, (uint256_t) child);

						// priority
						uint64_t priority_ith_node;
						for (uint64_t j = 0; j < nQ; j++) {
							uint64_t min_idx = Q[j]->get_index_pri(tupleQdags.roots[j] + Q[j]->getM(child));
							priority_ith_node = priorities[j][min_idx];

							if (type_priority_fun == TYPE_FUN_PRI_SUM) // sum
								total_weight += priority_ith_node;
							else if (type_priority_fun == TYPE_FUN_PRI_MAX) { // max
								if (total_weight < priority_ith_node) {
									total_weight = priority_ith_node;
								}
							}
						}
						// insert the last level node
						qdagWeight this_node = {next_level, nullptr, total_weight, newPath};
						pq.push(this_node); // add the tuple to the queue
					}
				}
				else {
					// Allocate 2D array for root_temp
					for (i = 0; i < children_to_recurse_size; ++i) {
						uint64_t* root_temp = new uint64_t[nQ];
						child = children_to_recurse[i];
						uint256_t newPath = tupleQdags.path;
						getNewMortonCodePath(newPath, nAtt, cur_level, (uint256_t) child);


						// compute the weight of the tuple (ONLY if it's not a leaf)
						// calculate the weight of the tuple
						uint64_t init, end, priority_ith_node;
						bool success;
						for (uint64_t j = 0; j < nQ; j++) {
							// we store the parent node that corresponds in the original quadtree of each qdag
							root_temp[j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);
							init = 0;
							end = priorities[j].size() - 1;
							priority_ith_node = 0;
							success = Q[j]->get_range_leaves(cur_level + 1, root_temp[j], init, end);
							if (success) {
								bit_vector::size_type min_idx = rMq[j](init, end);
								priority_ith_node = priorities[j][min_idx];
							} else {
								cout << "error in get range leaves and_ranked (louds)" << endl;
								exit(1);
							}
							if (type_priority_fun == TYPE_FUN_PRI_SUM) // sum
								total_weight += priority_ith_node;
							else if (type_priority_fun == TYPE_FUN_PRI_MAX) { // max
								if (total_weight < priority_ith_node) {
									total_weight = priority_ith_node;
								}
							}
						}
						// insert the tuple
						qdagWeight this_node = {next_level, root_temp, total_weight, newPath};
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
 * @param size_queue the size of the priority queue.
 * @param priorities the priorities of the points. We have a vector of priorities for each qdag.
 * @param rMq the rmq structure to get the minimum priority in a range.
 * @param results_points vector with the output
 * @return
 */
bool multiJoinRankedResults(
        vector<qdag> &Q,
        bool bounded_result,
        uint64_t UPPER_BOUND,
        uint8_t type_priority_fun,
        vector<int_vector<>> &priorities,
        vector<rmq_succinct_sct<false>> &rMq,
        vector<uint256_t> &results_points) {

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

    priority_queue<qdagWeight> pq;
    pq.push({0, Q_roots, 1, 0}); // insert the root of the qdag

    AND_ranked(Q_star, Q.size(),
               max_level, A.size(),
               bounded_result, UPPER_BOUND,
               pq, type_priority_fun,
               priorities, rMq, results_points);

//    cout << "number of results: " << results_points.size() << endl;
//    uint64_t i=0;
//    while(i < results_points.size()){
//		cout << results_points[i] << endl;
////        uint16_t* coordinates = results_points[i];
////        for(uint64_t j = 0; j< A.size(); j++){
////            cout << coordinates[j] << " ";
////        }
////        cout << endl;
//        i++;
//    }

    for (uint64_t i = 0; i < Q.size(); i++)
        delete Q_star[i];

    return true;
}



/**
 * * Computing the join in a certain order. After computing size_queue results, we will stop computing the results.
 * @param Q
 * @param nQ number of Qdags
 * @param nAtt number of attributes of the final quadtree.
 * @param roots
 * @param cur_level
 * @param max_level maximum level of the qdag
 * @param type_priority_fun 0 sum , 1 maximum. Only taking into account if partial_results is true.
 * @param coordinates
 * @param size_queue the size of the priority queue.
 * @param priorities the priorities of the points.
 * @param rMq the rmq structure to get the minimum priority in a range.
 * @param top_results the priority queue to store the top results (output).
 * @return true if we accomplished succesfully the join. Otherwise, return false.
 */
bool
AND_ranked_backtracking(
        qdag *Q[],
        uint16_t nQ,
        uint64_t nAtt,
        uint64_t *roots,
        uint16_t cur_level,
        uint16_t max_level,
        uint8_t type_priority_fun,
		uint256_t path,
        uint64_t size_queue,
        vector<int_vector<>> &priorities,
        vector<rmq_succinct_sct<false>> &rMq,
        priority_queue<qdagResults> &top_results) {

    uint64_t p = Q[0]->nChildren(); // number of children of the qdag (extended)
    bool just_zeroes = true;
    uint64_t k_d[nQ];
    uint16_t children_to_recurse[p];
    uint64_t i;
    uint64_t children_to_recurse_size = 0;
    uint32_t children = 0xffffffff; // each bit represent a node (empty or not)
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
		// we do not call recursively the function AND as we do in the other levels
		// add output to the priority queue of results
		for (i = 0; i < children_to_recurse_size; ++i){
			child = children_to_recurse[i];
			uint256_t newPath = path;
			getNewMortonCodePath(newPath, nAtt, cur_level, (uint256_t) child);

			// priority
			double this_weight = 0;
			uint64_t min_idx, priority_ith_node;
			for (uint64_t j = 0; j < nQ; j++) {
				min_idx = Q[j]->get_index_pri(roots[j] + Q[j]->getM(child));
				priority_ith_node = priorities[j][min_idx];

				if (type_priority_fun == TYPE_FUN_PRI_SUM) // sum
					this_weight += priority_ith_node;
				else if (type_priority_fun == TYPE_FUN_PRI_MAX) { // max
					if (this_weight < priority_ith_node) {
						this_weight = priority_ith_node;
					}
				}
			}
			// queue full --> compare priorities
			if(top_results.size() >= size_queue ){
				qdagResults minResult = top_results.top();
				// add result if the priority is higher than the minimum priority in the queue
				if(this_weight > minResult.weight){
					top_results.pop();
					top_results.push({newPath, this_weight});
				}
			}
			else{
				top_results.push({newPath, this_weight});
				just_zeroes = false;
			}
		}

    }
	// call recursively in DFS order
    else {
		uint64_t rank_vector[16 /*nQ*/][64];

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
		uint16_t next_level = cur_level + 1;

		uint64_t root_temp[children_to_recurse_size][16 /*nQ*/]; // CUIDADO, solo hasta 16 relaciones por query

		std::vector<std::pair<uint64_t, uint64_t>> order_to_traverse;
		uint256_t newPath[children_to_recurse_size];
		for (i = 0; i < children_to_recurse_size; ++i) {
			child = children_to_recurse[i]; // the position of the 1s in children
			newPath[i] = path;
			getNewMortonCodePath(newPath[i], nAtt, cur_level, (uint256_t) child);

			// compute the weight of the tuple
			double total_weight = 0;
			uint64_t init, end, priority_ith_node;
			bool success;
			for (uint64_t j = 0; j < nQ; j++) {
				root_temp[i][j] = k_d[j] * (rank_vector[j][Q[j]->getM(child)] - 1);
				init = 0;
				end = priorities[j].size()-1;
				priority_ith_node = 0;
				success = Q[j]->get_range_leaves(cur_level+1,root_temp[i][j],init,end);
				if(success){
					bit_vector::size_type min_idx = rMq[j](init, end);
					priority_ith_node = priorities[j][min_idx];
				} else {
					cout << "error in get range leaves and_backtracking (louds)" << endl;
					exit(1);
				}
				if (type_priority_fun == TYPE_FUN_PRI_SUM) // sum
					total_weight += priority_ith_node;
				else if (type_priority_fun == TYPE_FUN_PRI_MAX) { // max
					if (total_weight < priority_ith_node) {
						total_weight = priority_ith_node;
					}
				}
			}
			// queue full --> compare priorities
			if(top_results.size() >= size_queue ){
				qdagResults minResult = top_results.top();
				if(total_weight >= minResult.weight){ // if it has a higher priority, then we enter into this branch of the tree
					order_to_traverse.push_back({i,total_weight}); // add the tuple to the queue
				}
			}
			else{
				order_to_traverse.push_back({i,total_weight}); // add the tuple to the queue
			}
		}
		sortBySecond(order_to_traverse);
		for(i=0; i<children_to_recurse_size; i++){
				AND_ranked_backtracking(Q,
										nQ,
										nAtt,
										root_temp[order_to_traverse[i].first],
										next_level, max_level, type_priority_fun,
										newPath[order_to_traverse[i].first],
										size_queue,
										priorities, rMq,
										top_results);
			}

    }
    return !just_zeroes;
}

/**
 *
 * @param Q
 * @param type_priority_fun 0 sum , 1 maximum. Only taking into account if partial_results is true.
 * @param size_queue the size of the priority queue.
 * @param rMq the rmq structure to get the minimum priority in a range.
 * @param priorities the priorities of the points.
 * @param top_results the priority queue to store the top results (output).
 * @return
 */
bool multiJoinRankedResultsBacktracking(
        vector<qdag> &Q,
        uint8_t type_priority_fun,
        int64_t size_queue,
        vector<int_vector<>> &priorities,
        vector<rmq_succinct_sct<false>> &rMq,
        priority_queue<qdagResults>& top_results) {
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


    AND_ranked_backtracking(Q_star, Q.size(),
                            A.size(),
							Q_roots,
                            0, max_level,
                            type_priority_fun,
							0,
							size_queue,
                            priorities,
							rMq,
                            top_results);

//    uint64_t size_queue_top = top_results.size();
//    cout << "number of results: " << top_results.size() << endl;
//    for(uint64_t i=0; i<size_queue_top; i++){
//        qdagResults res = top_results.top();
//        top_results.pop();
//		cout << res.path << endl;
////        for(uint64_t k=0; k<A.size(); k++) {
////            cout << res.coordinates[k] << " ";
////        }
////        cout << endl;
//    }


    for (uint64_t i = 0; i < Q.size(); i++)
        delete Q_star[i];

    return true;
}



