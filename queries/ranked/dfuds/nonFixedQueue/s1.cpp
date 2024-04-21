#include "../../../../src/dfuds/join_ranked_results.cpp"


#include <fstream>
#include<bits/stdc++.h>
#include<ratio>
#include<chrono>
#include<ctime>

using namespace std::chrono;



high_resolution_clock::time_point start_select, stop_select;
double total_time_select = 0.0;       
duration<double> time_span_select;

#define AT_X1 0
#define AT_X2 1
#define AT_X3 2
#define AT_X4 3


std::vector<std::vector<uint64_t>>* read_relation(const std::string filename, uint16_t n_Atts)
{
    std::ifstream input_stream(filename); 
    uint64_t x;
    uint16_t i, j=0;
    
    std::vector<std::vector<uint64_t>>* relation;
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


uint64_t maximum_in_table(std::vector<std::vector<uint64_t>> &table, uint16_t n_columns, uint64_t max_temp)
{
    uint64_t i, j;
    
    for (i = 0; i < table.size(); i++) 
        for (j = 0; j < n_columns; j++)
            if (table[i][j] > max_temp)
                max_temp = table[i][j];
    
    
    return max_temp;
}


int main(int argc, char** argv)
{
    qdag_dfuds::att_set att_R;
    qdag_dfuds::att_set att_S;
    qdag_dfuds::att_set att_T;
    qdag_dfuds::att_set att_U;
    
    att_R.push_back(AT_X1); att_R.push_back(AT_X2); 
    att_S.push_back(AT_X2); att_S.push_back(AT_X3); 
    att_T.push_back(AT_X3); att_T.push_back(AT_X4); 
    att_U.push_back(AT_X4); att_U.push_back(AT_X1); 
    
    std::string strRel_R(argv[1]), strRel_S(argv[2]), strRel_T(argv[3]), strRel_U(argv[4]); 
    
    std::vector<std::vector<uint64_t>>* rel_R = read_relation(strRel_R, att_R.size());
    std::vector<std::vector<uint64_t>>* rel_S = read_relation(strRel_S, att_S.size());
    std::vector<std::vector<uint64_t>>* rel_T = read_relation(strRel_T, att_T.size());
    std::vector<std::vector<uint64_t>>* rel_U = read_relation(strRel_U, att_U.size());
    
    uint64_t grid_side = 0;
    
    grid_side = maximum_in_table(*rel_R, att_R.size(), grid_side);
    grid_side = maximum_in_table(*rel_S, att_S.size(), grid_side);
    grid_side = maximum_in_table(*rel_T, att_T.size(), grid_side);
    grid_side = maximum_in_table(*rel_U, att_U.size(), grid_side);
    
    grid_side++;

    qdag_dfuds qdag_dfuds_rel_R(*rel_R, att_R, grid_side, 2, att_R.size());
    qdag_dfuds qdag_dfuds_rel_S(*rel_S, att_S, grid_side, 2, att_S.size());
    qdag_dfuds qdag_dfuds_rel_T(*rel_T, att_T, grid_side, 2, att_T.size());
    qdag_dfuds qdag_dfuds_rel_U(*rel_U, att_U, grid_side, 2, att_U.size());

    cout << ((((float)qdag_dfuds_rel_R.size()*8) + ((float)qdag_dfuds_rel_S.size()*8) + ((float)qdag_dfuds_rel_T.size()*8) + ((float)qdag_dfuds_rel_U.size()*8) )/(rel_R->size()*2 + rel_S->size()*2 + rel_T->size()*2 + rel_U->size()*2)) << "\t";
    
    vector<qdag_dfuds> Q(4);

    Q[0] = qdag_dfuds_rel_R;
    Q[1] = qdag_dfuds_rel_S;
    Q[2] = qdag_dfuds_rel_T;
    Q[3] = qdag_dfuds_rel_U;
    
    qdag_dfuds *Join_Result;

    // read priorities from file
    std::ifstream data_file_R(argv[4]); // Abrir el archivo de datos
    std::ifstream data_file_S(argv[5]); // Abrir el archivo de datos
    std::ifstream data_file_T(argv[6]); // Abrir el archivo de datos
    std::ifstream data_file_U(argv[7]); // Abrir el archivo de datos
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
    }i=0;
    while(data_file_U >> value){
        priorities_U[i]=value;
        i++;
    }

    vector<int_vector<>> p;
    p.push_back(priorities_R);
    p.push_back(priorities_S);
    p.push_back(priorities_T);
    p.push_back(priorities_U);

    uint8_t type_fun = argv[9] ? atoi(argv[9]) : 1;
    vector<rmq_succinct_sct<false>> rMq;
    for(uint64_t i = 0; i < Q.size(); i++)
        rMq.push_back(rmq_succinct_sct<false>(&p[i]));
    vector<uint16_t*> results_ranked_louds;


    high_resolution_clock::time_point start, stop;
    double total_time = 0.0;       
    duration<double> time_span;

    multiJoinRankedResultsDfuds(Q, true, 1000, type_fun, p, rMq, results_ranked_louds);
    
    start = high_resolution_clock::now();

    multiJoinRankedResultsDfuds(Q, true, 1000, type_fun, p, rMq, results_ranked_louds);

    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();    

    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    return 0;
}
