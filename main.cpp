// cada una de las query tiene su codigo
// j --> una forma q sale en el paper
// p --> caminos
// s --> cuadrados
// t --> triangulos
#include<fstream>
#include<bits/stdc++.h>
#include<ratio>
#include<chrono>
#include<ctime>

using namespace std::chrono;

#include <algorithm>
#include "src/qdags.hpp"
#include "src/qdags_dfuds.hpp"
#include "src/joins.cpp"
#include "src/louds/join_partial_results.cpp"
#include "src/louds/join_ranked_results.cpp"
#include "src/dfuds/join_partial_results.cpp"
#include "src/dfuds/join_ranked_results.cpp"
#include "src/lqdag.hpp"
#include "src/lqdags/operations.cpp"


high_resolution_clock::time_point start_select, stop_select;
double total_time_select = 0.0;
duration<double> time_span_select;

#define AT_X1 0
#define AT_X2 1
#define AT_X3 2
#define AT_X4 3
#define AT_X5 4

#define AT_X 0
#define AT_Y 1
#define AT_Z 2
#define AT_V 3
#define AT_U 4


std::vector<std::vector<uint64_t>> *
read_relation(const std::string filename, uint16_t n_Atts)
{
    std::ifstream input_stream(filename);
    uint64_t x;
    uint16_t i, j = 0;

    std::vector<std::vector<uint64_t>> *relation;
    std::vector<uint64_t> tuple;

    relation = new std::vector<std::vector<uint64_t>>();

    input_stream >> x;
    while (!input_stream.eof()) {
        tuple.clear();
        for (i = 0; i < n_Atts; i++) {
            tuple.push_back(x);
            input_stream >> x;
        }

        relation->push_back(tuple);
    }

    return relation;
}


uint64_t maximum_in_table(std::vector<std::vector<uint64_t>> &table, uint16_t n_columns, uint64_t max_temp) {
    uint64_t i, j;

    for (i = 0; i < table.size(); i++)
        for (j = 0; j < n_columns; j++)
            if (table[i][j] > max_temp)
                max_temp = table[i][j];


    return max_temp;
}


int main(int argc, char **argv) {

    qdag::att_set att_R;
    qdag::att_set att_S;
    qdag::att_set att_T;
    qdag::att_set att_U;


//    att_R.push_back(AT_X1); att_R.push_back(AT_X2);
//    att_S.push_back(AT_X2); att_S.push_back(AT_X3);
//    att_T.push_back(AT_X3); att_T.push_back(AT_X4);
//    att_U.push_back(AT_X4); att_U.push_back(AT_X5);

	att_R.push_back(AT_X); att_R.push_back(AT_Y);
	att_S.push_back(AT_X); att_S.push_back(AT_Z);
	att_T.push_back(AT_X); att_T.push_back(AT_U);
	att_U.push_back(AT_X); att_U.push_back(AT_V);

	std::string strRel_R(argv[1]), strRel_S(argv[2]), strRel_T(argv[3]), strRel_U(argv[4]);


    // lee desde el disco la relacion R que tiene tal cantidad de atributoss --> con eso genero la relación r rel_R
    std::vector<std::vector<uint64_t>> *rel_R = read_relation(strRel_R,att_R.size()); // att_R.sizecantidad de atributos que tiene la relacion
    std::vector<std::vector<uint64_t>> *rel_R_2 = read_relation(strRel_R,att_R.size()); // att_R.sizecantidad de atributos que tiene la relacion
    std::vector<std::vector<uint64_t>>* rel_S = read_relation(strRel_S, att_S.size());
    std::vector<std::vector<uint64_t>>* rel_S_2 = read_relation(strRel_S, att_S.size());
    std::vector<std::vector<uint64_t>>* rel_T = read_relation(strRel_T, att_T.size());
    std::vector<std::vector<uint64_t>>* rel_T_2 = read_relation(strRel_T, att_T.size());
    std::vector<std::vector<uint64_t>> *rel_U = read_relation(strRel_U, att_U.size());
    std::vector<std::vector<uint64_t>> *rel_U_2 = read_relation(strRel_U, att_U.size());

//    uint64_t grid_side = 52000000; // es como +infty para wikidata
	uint64_t grid_side = 53000000;
//    uint64_t grid_side = 8;

    qdag qdag_rel_R(*rel_R, att_R, grid_side, 2, att_R.size()); // construyo los qdags
    qdag qdag_rel_S(*rel_S, att_S, grid_side, 2, att_S.size());
    qdag qdag_rel_T(*rel_T, att_T, grid_side, 2, att_T.size());
    qdag qdag_rel_U(*rel_U, att_U, grid_side, 2, att_U.size());
    qdag_dfuds qdag_rel_R_dfuds(*rel_R_2, att_R, grid_side, 2, att_R.size());
    qdag_dfuds qdag_rel_S_dfuds(*rel_S_2, att_S, grid_side, 2, att_S.size());
    qdag_dfuds qdag_rel_T_dfuds(*rel_T_2, att_T, grid_side, 2, att_T.size());
    qdag_dfuds qdag_rel_U_dfuds(*rel_U_2, att_U, grid_side, 2, att_U.size());

//     print the tree
//    cout << endl << "rel R" << endl;
//    qdag_rel_R.printBv();
//    cout << endl << "rel S" << endl;
//    qdag_rel_S.printBv();
//    cout << endl << "rel T" << endl;
//    qdag_rel_T.printBv();
//    cout << endl << "rel U" << endl;
//    qdag_rel_U.printBv();


    vector<qdag> Q(4);
//    vector<qdag> Q(2);
    vector<qdag_dfuds> Q_dfuds(4);
//    vector<qdag_dfuds> Q_dfuds(2);

    Q[0] = qdag_rel_R;
    Q[1] = qdag_rel_S;
    Q[2] = qdag_rel_T;
    Q[3] = qdag_rel_U;

    Q_dfuds[0] = qdag_rel_R_dfuds;
    Q_dfuds[1] = qdag_rel_S_dfuds;
    Q_dfuds[2] = qdag_rel_T_dfuds;
    Q_dfuds[3] = qdag_rel_U_dfuds;

    high_resolution_clock::time_point start, stop;
    double total_time = 0.0;
    duration<double> time_span;

    // read priorities from file
    std::ifstream data_file_R(argv[5]); // Abrir el archivo de datos
    std::ifstream data_file_S(argv[6]); // Abrir el archivo de datos
    std::ifstream data_file_T(argv[7]); // Abrir el archivo de datos
    std::ifstream data_file_U(argv[8]); // Abrir el archivo de datos
    if (!data_file_R.is_open() || !data_file_S.is_open() || !data_file_T.is_open() || !data_file_U.is_open()) {
        std::cerr << "No se pudo abrir el archivo de datos." << std::endl;
        return 1;
    }

    // get number of priorities of each file
    int number_of_lines_R = 0, number_of_lines_S = 0, number_of_lines_T = 0, number_of_lines_U = 0;
    std::string line;
    while(std::getline(data_file_R, line))
        ++number_of_lines_R;
    while(std::getline(data_file_S, line))
        ++number_of_lines_S;
    while(std::getline(data_file_T, line))
        ++number_of_lines_T;
    while(std::getline(data_file_U, line))
        ++number_of_lines_U;
    int_vector<> priorities_R(number_of_lines_R,0);
    int_vector<> priorities_S(number_of_lines_S,0);
    int_vector<> priorities_T(number_of_lines_T,0);
    int_vector<> priorities_U(number_of_lines_U,0);

    data_file_R.clear();
    data_file_R.seekg(0, std::ios::beg);
    data_file_S.clear();
    data_file_S.seekg(0, std::ios::beg);
    data_file_T.clear();
    data_file_T.seekg(0, std::ios::beg);
    data_file_U.clear();
    data_file_U.seekg(0, std::ios::beg);

    // put the priorities in the int_vector

    int value;
    int i=0;
    while(data_file_R >> value){
        priorities_R[i]=value;
        i++;
    }
    i=0;
    while(data_file_S >> value){
        priorities_S[i]=value;
        i++;
    }
    i=0;
    while(data_file_T >> value){
        priorities_T[i]=value;
        i++;
    }
    i=0;
    while(data_file_U >> value){
        priorities_U[i]=value;
        i++;
    }

    vector<int_vector<>> vector_pri;
    vector_pri.push_back(priorities_R);
    vector_pri.push_back(priorities_S);
    vector_pri.push_back(priorities_T);
    vector_pri.push_back(priorities_U);

    // function type
    uint8_t type_fun = argv[9] ? atoi(argv[9]) : 1;
    // size queue
    int64_t k = argv[10] ? atoi(argv[10]) : 1000;

    vector<rmq_succinct_sct<false>> rMq_louds; // vector of rMq for each qdag
    for(uint64_t i = 0; i < Q.size(); i++){
        rMq_louds.push_back(rmq_succinct_sct<false>(&vector_pri[i]));
    }

    vector<uint16_t*> results_partial_dfuds;
    vector<uint16_t*> results_partial_dfuds_back;
    vector<uint16_t*> results_ranked_dfuds;
    priority_queue<qdagResults> results_ranked_dfuds_back;
    vector<uint16_t*> results_partial_louds;
    vector<uint16_t*> results_partial_louds_back;
    vector<uint16_t*> results_ranked_louds;
    priority_queue<qdagResults> results_ranked_louds_back;

//	cout << "----- MULTI JOIN PARTIAL RESULTS BACKTRACKING ------" << endl;
////    multiJoinPartialResultsBacktracking(Q, grid_side, type_fun, k, results_partial_louds_back);
//    results_partial_louds_back.clear();
//    start = high_resolution_clock::now();
//    multiJoinPartialResultsBacktracking(Q, grid_side, type_fun, k, results_partial_louds_back);
//    stop = high_resolution_clock::now();
//    time_span = duration_cast<microseconds>(stop - start);
//    total_time = time_span.count();
//    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;
//
//    cout << "----- MULTI JOIN TRADICIONAL ------" << endl;
////    multiJoin(Q, true, k);
//    start = high_resolution_clock::now();
//    multiJoin(Q, true, k);
//    stop = high_resolution_clock::now();
//    time_span = duration_cast<microseconds>(stop - start);
//    total_time = time_span.count();
//    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;


    cout << "----- MULTI JOIN PARTIAL RESULTS DFUDS------" << endl;
//    multiJoinPartialResultsDfuds(Q_dfuds, true, k, grid_side, type_fun, results_partial_dfuds);
    results_partial_dfuds.clear();
    start = high_resolution_clock::now();
    multiJoinPartialResultsDfuds(Q_dfuds, true, k, grid_side, type_fun, results_partial_dfuds);
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;


//    cout << "----- MULTI JOIN PARTIAL RESULTS BACKTRACKING DFUDS------" << endl;
//    multiJoinPartialResultsDfudsBacktracking(Q_dfuds, grid_side, type_fun, k, results_partial_dfuds_back);
//    results_partial_dfuds_back.clear();
//    start = high_resolution_clock::now();
//    multiJoinPartialResultsDfudsBacktracking(Q_dfuds, grid_side, type_fun, k, results_partial_dfuds_back);
//    stop = high_resolution_clock::now();
//    time_span = duration_cast<microseconds>(stop - start);
//    total_time = time_span.count();
//    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;
//
//
//    cout << "----- MULTI JOIN RANKED RESULTS DFUDS------" << endl;
//    multiJoinRankedResultsDfuds(Q_dfuds,true, k, type_fun, vector_pri, rMq_louds, results_ranked_dfuds);
//    results_ranked_dfuds.clear();
//    start = high_resolution_clock::now();
//    multiJoinRankedResultsDfuds(Q_dfuds,true, k, type_fun, vector_pri, rMq_louds, results_ranked_dfuds);
//    stop = high_resolution_clock::now();
//    time_span = duration_cast<microseconds>(stop - start);
//    total_time = time_span.count();
//    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;


//    cout << "----- MULTI JOIN RANKED RESULTS BACKTRACKING DFUDS------" << endl;
//    multiJoinRankedResultsDfudsBacktracking(Q_dfuds, type_fun, k,  vector_pri, rMq_louds, results_ranked_dfuds_back);
//    results_ranked_dfuds_back = priority_queue<qdagResults>();
//    start = high_resolution_clock::now();
//    multiJoinRankedResultsDfudsBacktracking(Q_dfuds, type_fun, k,  vector_pri, rMq_louds, results_ranked_dfuds_back);
//    stop = high_resolution_clock::now();
//    time_span = duration_cast<microseconds>(stop - start);
//    total_time = time_span.count();
//    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;
//
//
//    cout << "----- MULTI JOIN PARTIAL RESULTS ------" << endl;
//    multiJoinPartialResults(Q, true, k, grid_side, type_fun, results_partial_louds); // warmup join -> activar el caché
//    results_partial_louds.clear();
//    start = high_resolution_clock::now();
//    multiJoinPartialResults(Q, true, k, grid_side, type_fun, results_partial_louds); // warmup join -> activar el caché
//    stop = high_resolution_clock::now();
//    time_span = duration_cast<microseconds>(stop - start);
//    total_time = time_span.count();
//    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;
//
//    cout << "----- MULTI JOIN PARTIAL RESULTS BACKTRACKING ------" << endl;
//    multiJoinPartialResultsBacktracking(Q, grid_side, type_fun, k, results_partial_louds_back);
//    results_partial_louds_back.clear();
//    start = high_resolution_clock::now();
//    multiJoinPartialResultsBacktracking(Q, grid_side, type_fun, k, results_partial_louds_back);
//    stop = high_resolution_clock::now();
//    time_span = duration_cast<microseconds>(stop - start);
//    total_time = time_span.count();
//    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

//    cout << "----- MULTI JOIN RANKED RESULTS ------" << endl;
//    multiJoinRankedResults(Q, true, k, type_fun, vector_pri, rMq_louds, results_ranked_louds); // warmup join -> activar el caché
//    results_ranked_louds.clear();
//    start = high_resolution_clock::now();
//    multiJoinRankedResults(Q, true, k, type_fun, vector_pri, rMq_louds, results_ranked_louds); // warmup join -> activar el caché
//    stop = high_resolution_clock::now();
//    time_span = duration_cast<microseconds>(stop - start);
//    total_time = time_span.count();
//    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    cout << "----- MULTI JOIN RANKED RESULTS BACKTRACKING------" << endl;
//    multiJoinRankedResultsBacktracking(Q, type_fun, k, vector_pri, rMq_louds, results_ranked_louds_back); // warmup join -> activar el caché
    results_ranked_louds_back = priority_queue<qdagResults>();
    start = high_resolution_clock::now();
    multiJoinRankedResultsBacktracking(Q, type_fun, k, vector_pri, rMq_louds, results_ranked_louds_back); // warmup join -> activar el caché
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

//    cout << "----- MULTI JOIN LQDAGS ------" << endl;
//
//    uint64_t res = 0; //x1,x2,x3
//    predicate pred1 = {AT_X1, AT_X3, 0, OP_GREATER,TYPE_ATT1_ATT2};
//    predicate pred2 = {AT_X1, AT_X3, 3, OP_GREATER_EQUAL,TYPE_ATT1_ATT2};
//    predicate pred3 = {AT_X1, AT_X3, 3, OP_EQUAL,TYPE_ATT1_ATT2};
//    predicate pred4 = {AT_X2, AT_X3, 3, OP_GREATER,TYPE_ATT1_CONST};
//    predicate pred5 = {AT_X3, AT_X3, 3, OP_GREATER,TYPE_ATT1_CONST};
//    predicate pred6 = {AT_X3, AT_X3, 0, OP_EQUAL,TYPE_ATT1_CONST};
//    predicate pred7 = {AT_X3, AT_X3, 1, OP_EQUAL,TYPE_ATT1_CONST};
//    predicate pred12 = {AT_X1, AT_X3, 0, OP_LESS_EQUAL,TYPE_ATT1_ATT2};
//    start = high_resolution_clock::now();
//
//    quadtree_formula* join_r_s_t = compute_lqdag_join(Q, k, res);
//    cout << "number of results join: " << res << endl;
//
//    stop = high_resolution_clock::now();
//    time_span = duration_cast<microseconds>(stop - start);
//    total_time = time_span.count();
//    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;
//
//    res = 0;
//    uint64_t p = std::pow(Q[0].getK(), 3);
//    quadtree_formula result{};
//    lqdag* join_r_s_lqdag = lqdag_join(Q, k);
//    lqdag* child_0 = join_r_s_lqdag->lazy_child_completion(p,0,result);
//    lqdag* child_0_0 = child_0->lazy_child_completion(p,0,*result.children[0]);
//    lqdag* test_fin = join_r_s_lqdag->lazy_child_completion(p,0,*join_r_s_t);
//
//    res = 0;
//    compute_pred_lqdag_join(Q,k,res,&pred1);
//    cout << "number of results join x1 > x3: " << res << endl;
//
//    res = 0;
//    compute_pred_lqdag_join(Q,k,res,&pred2);
//    cout << "number of results join x1 >= x3: " << res << endl;
//
//    res = 0;
//    compute_pred_lqdag_join(Q,k,res,&pred3);
//    cout << "number of results join x1 = x3: " << res << endl;
//
//    res = 0;
//    compute_pred_lqdag_join(Q,k,res,&pred4);
//    cout << "number of results join x2 > 4: " << res << endl;
//
//    res = 0;
//    compute_pred_lqdag_join(Q,k,res,&pred5);
//    cout << "number of results join x3 > 4: " << res << endl;
//
//    res = 0;
//    compute_pred_lqdag_join(Q,k,res,&pred6);
//    cout << "number of results join x3 = 0: " << res << endl;
//
//    res = 0;
//    compute_pred_lqdag_join(Q,k,res,&pred7);
//    cout << "number of results join x3 = 1: " << res << endl;
//
//    res = 0;
//    quadtree_formula* test_pred_1 = compute_pred_quadtree(qdag_rel_R, k, res, &pred1);
//    cout << "number of results join A=B: " << res << endl;
//
//    res = 0;
//    quadtree_formula* test_pred_2 = compute_pred_quadtree(qdag_rel_R, k, res, &pred2);
//    cout << "number of results join A >= 3: " << res << endl;
//
//    res = 0;
//    quadtree_formula* test_pred_3 = compute_pred_quadtree(qdag_rel_R, k, res, &pred3);
//    cout << "number of results join A>3: " << res << endl;
//
//    res = 0;
//    quadtree_formula* test_pred_4 = compute_pred_quadtree(qdag_rel_R, k, res, &pred4);
//    cout << "number of results join B>3: " << res << endl;
//
//    res = 0;
//    quadtree_formula* compute_rel_R = compute_quadtree(qdag_rel_R, k, res);
//    cout << "number of results compute R: " << res << endl;
//
//    res = 0;
//    quadtree_formula* compute_rel_S = compute_quadtree(qdag_rel_S, k, res);
//    cout << "number of results compute S: " << res << endl;
//
//    res = 0;
//    quadtree_formula* compute_rel_T = compute_quadtree(qdag_rel_T, k, res);
//    cout << "number of results compute T: " << res << endl;
//
//    res = 0;
//    quadtree_formula* compute_rel_U = compute_quadtree(qdag_rel_U, k, res);
//    cout << "number of results compute U: " << res << endl;
//
////    res = 0;
////    quadtree_formula* pred_a_c = compute_pred_lqdag_join(Q,k,res,&pred12);
////    cout << "number of results pred a_c: " << res << endl;






    return 0;

}
