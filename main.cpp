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
#include "src/join_partial_results.cpp"
#include "src/join_ranked_results.cpp"


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

    // only for testing
    att_R.push_back(AT_Y); att_R.push_back(AT_X);
    att_S.push_back(AT_Y); att_S.push_back(AT_Z);// TODO: el orden hay q cambiarlo para el join


    /*att_R.push_back(AT_Y); att_R.push_back(AT_X);
    att_S.push_back(AT_Z); att_S.push_back(AT_X);
    att_T.push_back(AT_X); att_T.push_back(AT_V);*/

    std::string strRel_R(argv[1]), strRel_S(argv[2]), strRel_T(argv[3]); // nombre de los archivos

    // lee desde el disco la relacion R que tiene tal cantidad de atributoss --> con eso genero la relación r rel_R
    std::vector<std::vector<uint64_t>> *rel_R = read_relation(strRel_R,att_R.size()); // att_R.sizecantidad de atributos que tiene la relacion
    std::vector<std::vector<uint64_t>>* rel_S = read_relation(strRel_S, att_S.size());
    //std::vector<std::vector<uint64_t>>* rel_T = read_relation(strRel_T, att_T.size());

    uint64_t grid_side =32;// 52000000; // es como +infty para wikidata

    qdag_dfuds qdag_rel_R_dfuds(*rel_R, att_R, grid_side, 2, att_R.size());
    //cout << " grid size R : " << att_R.size() << endl;
    //cout << " grid size S : " << att_S.size() << endl;
    qdag qdag_rel_R(*rel_R, att_R, grid_side, 2, att_R.size()); // construyo los qdags
    qdag qdag_rel_S(*rel_S, att_S, grid_side, 2, att_S.size());
    //qdag qdag_rel_T(*rel_T, att_T, grid_side, 2, att_T.size());*/


    cout << endl << "rel R" << endl;
    qdag_rel_R.printBv();
    cout << endl << "rel S" << endl;
    //qdag_rel_S.printBv();

    // cout << ((((float)qdag_rel_R.size()*8) + ((float)qdag_rel_S.size()*8) + ((float)qdag_rel_T.size()*8) )/(rel_R->size()*2 + rel_S->size()*2 + rel_T->size()*2)) << "\t";
    //vector<qdag> Q(3);
    vector<qdag> Q(2);

    Q[0] = qdag_rel_R;
    Q[1] = qdag_rel_S;
    //Q[2] = qdag_rel_T;
    qdag *Join_Result;

    high_resolution_clock::time_point start, stop;
    double total_time = 0.0;
    duration<double> time_span;

    // read priorities from file
    /*std::ifstream data_file_R(argv[4]); // Abrir el archivo de datos
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

    vector<int_vector<>> p;
    p.push_back(priorities_R);
    p.push_back(priorities_S);
    p.push_back(priorities_T);*/

    /*
    int_vector<> prioritiesR={6,2,7,3,2,4,2,1,5,8,2,10,1,15,16,1,22,3,4,5,2,1,0,0,50,4,2,1,5,8,2,10,1,1,1,1,22,3,4,5,2,1,0,0,50};
    int_vector<> prioritiesS={1,1,1,1,4,1,200,3,31,22,3,4,5,2,1,15,27,50,4,5,2,1,0,4,2,1,5,8,2,10,1,1,1,1,22,3,4,5,2,1,0,0,50};
    vector<int_vector<>> p2;
    p2.push_back(prioritiesR);
    p2.push_back(prioritiesS);*/




    start = high_resolution_clock::now();
    cout << "----- MULTI JOIN TRADICIONAL ------" << endl;
    //Join_Result = multiJoin(Q, true, 1000);
    cout << "----- MULTI JOIN PARTIAL RESULTS ------" << endl;
    //multiJoinPartialResults(Q, true, 1000, 0, 0);
    cout << "----- MULTI JOIN PARTIAL RESULTS BACKTRACKING------" << endl;
    //multiJoinPartialResultsBacktracking(Q, 0, 50, 100);

    // PARTIAL JOIN
    // vector<qdag> &Q, bool bounded_result, uint64_t UPPER_BOUND,
    //                                  bool partial_results, int16_t type_priority_fun, int16_t type_order_fun, uint64_t grid_size)
    //bool join = multiJoinPartialResults(Q, true, 1000, true, 1, 1, grid_side, -1); // warmup join -> activar el caché
    //bool join = multiJoinPartialResults(Q, true, 1000, 1, grid_side, 6); // warmup join -> activar el caché
    /*int_vector<> prioritiesR={6,2,7,3,2,4,2,1,5,8,2,10,1,1,1,1,22,3,4,5,2,1,0,0,50};
    int_vector<> prioritiesS={1,1,1,1,4,1,1,1,1,1,22,3,4,5,2,1,0,0,50,4,5,2,1,0,};
    vector<int_vector<>> p;
    p.push_back(prioritiesR);
    p.push_back(prioritiesS);*/
    //qdag_rel_R.create_dfuds();
    //qdag_rel_R.printBv();

    cout << "----- MULTI JOIN RANKED RESULTS ------" << endl;
    //bool join = multiJoinRankedResults(Q, true, 1000, 1, -1, p); // warmup join -> activar el caché
    cout << "----- MULTI JOIN RANKED RESULTS BACKTRACLOMG------" << endl;
    //multiJoinRankedResultsBacktracking(Q, 1, 10, p2); // warmup join -> activar el caché

    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();



    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    return 0;

    // en runquery le paso cada una de las 3 tablas
    // no entiendo los valores de las tablas, aah si son el mismo numero entonces se comparte, pero como
}
