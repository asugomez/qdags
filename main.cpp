// cada una de las query tiene su codigo
// j --> una forma q sale en el paper
// p --> caminos
// s --> cuadrados
// t --> triangulos
#include <fstream>
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

#define AT_X 0
#define AT_Y 1
#define AT_Z 2
#define AT_V 3


std::vector<std::vector<uint64_t>> *
read_relation(const std::string filename, uint16_t n_Atts) // TODO: ver si hay una opcion mas optima
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
    // 3 tablas R, S y T
    qdag::att_set att_R; // R(Y,X)
    qdag::att_set att_S; // S(Z,X)
    qdag::att_set att_T; // T(X,V)

    att_R.push_back(AT_Y); att_R.push_back(AT_X);
    att_S.push_back(AT_Z); att_S.push_back(AT_X);
    att_T.push_back(AT_X); att_T.push_back(AT_V);

    // only for testing
//    att_R.push_back(AT_Y); att_R.push_back(AT_X);
//    att_S.push_back(AT_Y); att_S.push_back(AT_Z);// TODO: el orden hay q cambiarlo para el join
//    att_S.push_back(AT_X); att_S.push_back(AT_Z);// TODO: el orden hay q cambiarlo para el join



    qdag::att_set att_A;
    att_A.push_back(AT_X); att_A.push_back(AT_Y); att_A.push_back(AT_Z); att_A.push_back(AT_V);
//    att_A.push_back(AT_X); att_A.push_back(AT_Y); att_A.push_back(AT_Z);

    std::string strRel_R(argv[1]), strRel_S(argv[2]), strRel_T(argv[3]); // nombre de los archivos

    // lee desde el disco la relacion R que tiene tal cantidad de atributoss --> con eso genero la relación r rel_R
    std::vector<std::vector<uint64_t>> *rel_R = read_relation(strRel_R,att_R.size()); // att_R.sizecantidad de atributos que tiene la relacion
    std::vector<std::vector<uint64_t>> *rel_R_2 = read_relation(strRel_R,att_R.size()); // att_R.sizecantidad de atributos que tiene la relacion
    std::vector<std::vector<uint64_t>>* rel_S = read_relation(strRel_S, att_S.size());
    std::vector<std::vector<uint64_t>>* rel_S_2 = read_relation(strRel_S, att_S.size());
    std::vector<std::vector<uint64_t>>* rel_T = read_relation(strRel_T, att_T.size());
    std::vector<std::vector<uint64_t>>* rel_T_2 = read_relation(strRel_T, att_T.size());

    uint64_t grid_side = 52000000; // es como +infty para wikidata
//    uint64_t grid_side = 32;

    qdag qdag_rel_R(*rel_R, att_R, grid_side, 2, att_R.size()); // construyo los qdags
    qdag qdag_rel_S(*rel_S, att_S, grid_side, 2, att_S.size());
    qdag qdag_rel_T(*rel_T, att_T, grid_side, 2, att_T.size());
    qdag_dfuds qdag_rel_R_dfuds(*rel_R_2, att_R, grid_side, 2, att_R.size());
    qdag_dfuds qdag_rel_S_dfuds(*rel_S_2, att_S, grid_side, 2, att_S.size());
    qdag_dfuds qdag_rel_T_dfuds(*rel_T_2, att_T, grid_side, 2, att_T.size());

//     print the tree
//    cout << endl << "rel R" << endl;
//    qdag_rel_R.printBv();
//    cout << endl << "rel S" << endl;
//    qdag_rel_S.printBv();
//    cout << endl << "rel T" << endl;
//    qdag_rel_T.printBv();


    vector<qdag> Q(3);
//    vector<qdag> Q(2);
    vector<qdag_dfuds> Q_dfuds(3);
//    vector<qdag_dfuds> Q_dfuds(2);

    Q[0] = qdag_rel_R;
    Q[1] = qdag_rel_S;
    Q[2] = qdag_rel_T;


    Q_dfuds[0] = qdag_rel_R_dfuds;
    Q_dfuds[1] = qdag_rel_S_dfuds;
    Q_dfuds[2] = qdag_rel_T_dfuds;

    high_resolution_clock::time_point start, stop;
    double total_time = 0.0;
    duration<double> time_span;

    // read priorities from file
    std::ifstream data_file_R(argv[4]); // Abrir el archivo de datos
    std::ifstream data_file_S(argv[5]); // Abrir el archivo de datos
    std::ifstream data_file_T(argv[6]); // Abrir el archivo de datos
    if (!data_file_R.is_open() || !data_file_S.is_open() || !data_file_T.is_open()) {
        std::cerr << "No se pudo abrir el archivo de datos." << std::endl;
        return 1;
    }

    // get number of priorities of each file
    int number_of_lines_R = 0, number_of_lines_S = 0, number_of_lines_T = 0;
    std::string line;
    while(std::getline(data_file_R, line))
        ++number_of_lines_R;
    while(std::getline(data_file_S, line))
        ++number_of_lines_S;
    while(std::getline(data_file_T, line))
        ++number_of_lines_T;
    int_vector<> priorities_R(number_of_lines_R,0);
    int_vector<> priorities_S(number_of_lines_S,0);
    int_vector<> priorities_T(number_of_lines_T,0);

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

    vector<int_vector<>> vector_pri;
    vector_pri.push_back(priorities_R);
    vector_pri.push_back(priorities_S);
    vector_pri.push_back(priorities_T);

    // size queue
    int64_t size_queue = argv[7] ? atoi(argv[7]) : 1000;

    vector<rmq_succinct_sct<false>> rMq_louds; // vector of rMq for each qdag
//    vector<rmq_succinct_sct<false>> rMq_louds_back; // vector of rMq for each qdag
//    vector<rmq_succinct_sct<false>> rMq_dfuds; // vector of rMq for each qdag
//    vector<rmq_succinct_sct<false>> rMq_dfuds_back; // vector of rMq for each qdag
    for(uint64_t i = 0; i < Q.size(); i++){
        rMq_louds.push_back(rmq_succinct_sct<false>(&vector_pri[i]));
//        rMq_louds_back.push_back(rmq_succinct_sct<false>(&p2[i]));
//        rMq_dfuds.push_back(rmq_succinct_sct<false>(&vector_pri[i]));
//        rMq_dfuds_back.push_back(rmq_succinct_sct<false>(&p2[i]));
    }


    vector<uint16_t*> results_partial_dfuds;
    vector<uint16_t*> results_partial_dfuds_back;
    vector<uint16_t*> results_ranked_dfuds;
    priority_queue<qdagResults> results_ranked_dfuds_back;
    vector<uint16_t*> results_partial_louds;
    vector<uint16_t*> results_partial_louds_back;
    vector<uint16_t*> results_ranked_louds;
    priority_queue<qdagResults> results_ranked_louds_back;

    cout << "----- MULTI JOIN TRADICIONAL R S T------" << endl;
//    multiJoin(Q, true, 1000);
    start = high_resolution_clock::now();
    multiJoin(Q, true, 1000);
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;


    cout << "----- MULTI JOIN PARTIAL RESULTS DFUDS------" << endl;
    multiJoinPartialResultsDfuds(Q_dfuds, true, 1000, grid_side, 1, results_partial_dfuds);
    start = high_resolution_clock::now();
    multiJoinPartialResultsDfuds(Q_dfuds, true, 1000, grid_side, 1, results_partial_dfuds);
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;


    cout << "----- MULTI JOIN PARTIAL RESULTS BACKTRACKING DFUDS------" << endl;
    multiJoinPartialResultsDfudsBacktracking(Q_dfuds, grid_side, 1, size_queue, results_partial_dfuds_back);
    start = high_resolution_clock::now();
    multiJoinPartialResultsDfudsBacktracking(Q_dfuds, grid_side, 1, size_queue, results_partial_dfuds_back);
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;


    cout << "----- MULTI JOIN RANKED RESULTS DFUDS------" << endl;
    multiJoinRankedResultsDfuds(Q_dfuds,true, 1000, 1, vector_pri, rMq_louds, results_ranked_dfuds);
    start = high_resolution_clock::now();
    multiJoinRankedResultsDfuds(Q_dfuds,true, 1000, 1, vector_pri, rMq_louds, results_ranked_dfuds);
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;


    cout << "----- MULTI JOIN RANKED RESULTS BACKTRACKING DFUDS------" << endl;
    multiJoinRankedResultsDfudsBacktracking(Q_dfuds, 1, size_queue,  vector_pri, rMq_louds, results_ranked_dfuds_back);
    start = high_resolution_clock::now();
    multiJoinRankedResultsDfudsBacktracking(Q_dfuds, 1, size_queue,  vector_pri, rMq_louds, results_ranked_dfuds_back);
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;


    cout << "----- MULTI JOIN PARTIAL RESULTS ------" << endl;
    multiJoinPartialResults(Q, true, 1000, grid_side, 1, results_partial_louds); // warmup join -> activar el caché
    start = high_resolution_clock::now();
    multiJoinPartialResults(Q, true, 1000, grid_side, 1, results_partial_louds); // warmup join -> activar el caché
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    cout << "----- MULTI JOIN PARTIAL RESULTS BACKTRACKING ------" << endl;
    multiJoinPartialResultsBacktracking(Q, grid_side, 1, size_queue, results_partial_louds_back);
    start = high_resolution_clock::now();
    multiJoinPartialResultsBacktracking(Q, grid_side, 1, size_queue, results_partial_louds_back);
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    cout << "----- MULTI JOIN RANKED RESULTS ------" << endl;
    multiJoinRankedResults(Q, true, 1000, 1, vector_pri, rMq_louds, results_ranked_louds); // warmup join -> activar el caché
    start = high_resolution_clock::now();
    multiJoinRankedResults(Q, true, 1000, 1, vector_pri, rMq_louds, results_ranked_louds); // warmup join -> activar el caché
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    cout << "----- MULTI JOIN RANKED RESULTS BACKTRACKING------" << endl;
    multiJoinRankedResultsBacktracking(Q, 1, size_queue, vector_pri, rMq_louds, results_ranked_louds_back); // warmup join -> activar el caché
    start = high_resolution_clock::now();
    multiJoinRankedResultsBacktracking(Q, 1, size_queue, vector_pri, rMq_louds, results_ranked_louds_back); // warmup join -> activar el caché
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    cout << "----- MULTI JOIN LQDAGS ------" << endl;
    uint64_t res = 0;
    uint64_t k = 100;
    uint64_t p = std::pow(qdag_rel_R.getK(), att_A.size());
    lqdag* join_r_s_t = lqdag_join(Q, true, k);

    start = high_resolution_clock::now();
    quadtree_formula* test_join = join_r_s_t->completion(p,
                                                         qdag_rel_R.getHeight(),
                                                         0, k,
                                                         res);
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();
    cout << "number of results: " << res << endl;
    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    return 0;

}
