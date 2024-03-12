#include "../../../../src/join_ranked_results.cpp"


#include <fstream>
#include<bits/stdc++.h>
#include<ratio>
#include<chrono>
#include<ctime>

using namespace std::chrono;




high_resolution_clock::time_point start_select, stop_select;
double total_time_select = 0.0;       
duration<double> time_span_select;

#define AT_X 0
#define AT_Y 1
#define AT_Z 2
#define AT_U 3


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
    qdag::att_set att_R;
    qdag::att_set att_S;
    qdag::att_set att_T;
    
    att_R.push_back(AT_Y); att_R.push_back(AT_X); 
    att_S.push_back(AT_Z); att_S.push_back(AT_X); 
    att_T.push_back(AT_U); att_T.push_back(AT_X); 
    
    std::string strRel_R(argv[1]), strRel_S(argv[2]), strRel_T(argv[3]); 
    
    std::vector<std::vector<uint64_t>>* rel_R = read_relation(strRel_R, att_R.size());
    std::vector<std::vector<uint64_t>>* rel_S = read_relation(strRel_S, att_S.size());
    std::vector<std::vector<uint64_t>>* rel_T = read_relation(strRel_T, att_T.size());
    
    uint64_t grid_side = 52000000;
    
    //grid_side = maximum_in_table(*rel_R, att_R.size(), grid_side);
    //grid_side = maximum_in_table(*rel_S, att_S.size(), grid_side);
    //grid_side = maximum_in_table(*rel_T, att_T.size(), grid_side);
    
    //grid_side++;
  
    qdag qdag_rel_R(*rel_R, att_R, grid_side, 2, att_R.size());
    qdag qdag_rel_S(*rel_S, att_S, grid_side, 2, att_S.size());
    qdag qdag_rel_T(*rel_T, att_T, grid_side, 2, att_T.size());
    
    // cout << ((((float)qdag_rel_R.size()*8) + ((float)qdag_rel_S.size()*8) + ((float)qdag_rel_T.size()*8) )/(rel_R->size()*2 + rel_S->size()*2 + rel_T->size()*2)) << "\t";

    vector<qdag> Q(3);

    Q[0] = qdag_rel_R;
    Q[1] = qdag_rel_S;
    Q[2] = qdag_rel_T;
    
    qdag *Join_Result;

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

    vector<int_vector<>> p;
    p.push_back(priorities_R);
    p.push_back(priorities_S);
    p.push_back(priorities_T);

    // size queue
    int64_t size_queue = argv[7] ? atoi(argv[7]) : 100;
  
    high_resolution_clock::time_point start, stop;
    double total_time = 0.0;       
    duration<double> time_span;

    multiJoinRankedResults(Q, true, 1000, 1, size_queue, p);
 
    start = high_resolution_clock::now();

    multiJoinRankedResults(Q, true, 1000, 1, size_queue, p);

    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();    

    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    return 0;
}
