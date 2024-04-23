#include "../../../../src/louds/join_ranked_results.cpp"


#include <fstream>
#include<bits/stdc++.h>
#include<ratio>
#include<chrono>
#include<ctime>

using namespace std::chrono;



double total_time_select = 0.0;       
duration<double> time_span_select;

#define AT_X1 0
#define AT_X2 1
#define AT_X3 2


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
    
    att_R.push_back(AT_X1); att_R.push_back(AT_X2); 
    att_S.push_back(AT_X2); att_S.push_back(AT_X3); 
    
    std::string strRel_R(argv[1]), strRel_S(argv[2]); 
    
    std::vector<std::vector<uint64_t>>* rel_R = read_relation(strRel_R, att_R.size());
    std::vector<std::vector<uint64_t>>* rel_S = read_relation(strRel_S, att_S.size());
    
    uint64_t grid_side = 0;
    
    grid_side = maximum_in_table(*rel_R, att_R.size(), grid_side);
    grid_side = maximum_in_table(*rel_S, att_S.size(), grid_side);
    
    grid_side++;

    //cout << "Grid side: " << grid_side << endl;
    
    qdag qdag_rel_R(*rel_R, att_R, grid_side, 2, att_R.size());
    qdag qdag_rel_S(*rel_S, att_S, grid_side, 2, att_S.size());
    
    // cout << ((((float)qdag_rel_R.size()*8) + ((float)qdag_rel_S.size()*8) )/(rel_R->size()*2 + rel_S->size()*2)) << "\t";

    vector<qdag> Q(2);

    Q[0] = qdag_rel_R;
    Q[1] = qdag_rel_S;
    
    qdag *Join_Result;

    // read priorities from file
    std::ifstream data_file_R(argv[3]); // Abrir el archivo de datos
    std::ifstream data_file_S(argv[4]); // Abrir el archivo de datos
    if (!data_file_R.is_open() || !data_file_S.is_open()) {
        std::cerr << "No se pudo abrir el archivo de datos." << std::endl;
        return 1;
    }

    // get number of priorities of each file
    int number_of_lines_R = 0, number_of_lines_S = 0;
    std::string line;
    while(std::getline(data_file_R, line))
        ++number_of_lines_R;
    while(std::getline(data_file_S, line))
        ++number_of_lines_S;
    int_vector<> priorities_R(number_of_lines_R,0);
    int_vector<> priorities_S(number_of_lines_S,0);

    data_file_R.clear();
    data_file_R.seekg(0, std::ios::beg);
    data_file_S.clear();
    data_file_S.seekg(0, std::ios::beg);

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

    vector<int_vector<>> p;
    p.push_back(priorities_R);
    p.push_back(priorities_S);

    uint8_t type_fun = argv[5] ? atoi(argv[5]) : 1;
    vector<rmq_succinct_sct<false>> rMq;
    for(uint64_t i = 0; i < Q.size(); i++)
        rMq.push_back(rmq_succinct_sct<false>(&p[i]));
    int64_t size_queue = argv[6] ? atoi(argv[6]) : 1000;
    priority_queue<qdagResults> results_ranked_louds_back;

  
    high_resolution_clock::time_point start, stop;
    double total_time = 0.0;       
    duration<double> time_span;

    multiJoinRankedResultsBacktracking(Q, type_fun, size_queue, p, rMq, results_ranked_louds_back);

    results_ranked_louds_back = priority_queue<qdagResults>();

    start = high_resolution_clock::now();

    multiJoinRankedResultsBacktracking(Q, type_fun, size_queue, p, rMq, results_ranked_louds_back);


    stop = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(stop - start);
    total_time = time_span.count();    

    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    return 0;
}
