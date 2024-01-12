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
#include "src/joins.cpp"


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
    qdag::att_set att_R; // R(A,B)
    qdag::att_set att_S; // S(C,B)
    qdag::att_set att_T; // T(A,C)

    att_R.push_back(AT_Y); att_R.push_back(AT_X);
    att_S.push_back(AT_Z); att_S.push_back(AT_X);
    att_T.push_back(AT_Z); att_T.push_back(AT_Y); // TODO: el orden hay q cambiarlo para el join

    std::string strRel_R(argv[1]), strRel_S(argv[2]), strRel_T(argv[3]); // nombre de los archivos

    // lee desde el disco la relacion R que tiene tal cantidad de atributoss --> con eso genero la relación r rel_R
    std::vector<std::vector<uint64_t>> *rel_R = read_relation(strRel_R,
                                                              att_R.size()); // att_R.sizecantidad de atributos que tiene la relacion
    std::vector<std::vector<uint64_t>>* rel_S = read_relation(strRel_S, att_S.size());
    std::vector<std::vector<uint64_t>>* rel_T = read_relation(strRel_T, att_T.size());

    uint64_t grid_side = 8; //52000000; // es como +infty para wikidata

    //cout << "R" << endl;
    qdag qdag_rel_R(*rel_R, att_R, grid_side, 2, att_R.size()); // construyo los qdags
    //cout << "S" << endl;
    qdag qdag_rel_S(*rel_S, att_S, grid_side, 2, att_S.size());
    //cout << "T" << endl;
    qdag qdag_rel_T(*rel_T, att_T, grid_side, 2, att_T.size());

    // cout << ((((float)qdag_rel_R.size()*8) + ((float)qdag_rel_S.size()*8) + ((float)qdag_rel_T.size()*8) )/(rel_R->size()*2 + rel_S->size()*2 + rel_T->size()*2)) << "\t";
    vector<qdag> Q(3);

    Q[0] = qdag_rel_R;
    Q[1] = qdag_rel_S;
    Q[2] = qdag_rel_T;
    qdag *Join_Result;

    high_resolution_clock::time_point start, stop;
    double total_time = 0.0;
    duration<double> time_span;


    start = high_resolution_clock::now();

    //qdag_rel_R.printBv();
    cout << "-----------" << endl;
    cout << "relacion R" << endl;
    qdag_rel_R.printBv();
    qdag_rel_R.get_num_leaves(1, 1);
    /*cout << "-----------" << endl;
    cout << "relacion S" << endl;
    qdag_rel_S.printBv();
    qdag_rel_S.get_num_leaves(0, 0);
    cout << "-----------" << endl;
    cout << "relacion T" << endl;
    qdag_rel_T.printBv();
    qdag_rel_T.get_num_leaves(0, 0);*/
    //qdag_rel_R.get_children(0,0);
    cout << "-----------" << endl;
    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();

    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    return 0;

    // en runquery le paso cada una de las 3 tablas
    // no entiendo los valores de las tablas, aah si son el mismo numero entonces se comparte, pero como
}
