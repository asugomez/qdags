#include "../../../../src/dfuds/join_partial_results.cpp"

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
#define AT_V 3
#define AT_U 4

std::vector<std::vector<uint64_t>> *read_relation(const std::string filename, uint16_t n_Atts) {
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
    qdag_dfuds::att_set att_R;
    qdag_dfuds::att_set att_S;
    qdag_dfuds::att_set att_T;
    qdag_dfuds::att_set att_U;

    att_R.push_back(AT_Y);
    att_R.push_back(AT_X);
    att_S.push_back(AT_Z);
    att_S.push_back(AT_X);
    att_T.push_back(AT_X);
    att_T.push_back(AT_U);
    att_U.push_back(AT_X);
    att_U.push_back(AT_V);

    std::string strRel_R(argv[1]), strRel_S(argv[2]), strRel_T(argv[3]), strRel_U(argv[4]);

    std::vector<std::vector<uint64_t>> *rel_R = read_relation(strRel_R, att_R.size());
    std::vector<std::vector<uint64_t>> *rel_S = read_relation(strRel_S, att_S.size());
    std::vector<std::vector<uint64_t>> *rel_T = read_relation(strRel_T, att_T.size());
    std::vector<std::vector<uint64_t>> *rel_U = read_relation(strRel_U, att_U.size());

    uint64_t grid_side = 0;

    grid_side = maximum_in_table(*rel_R, att_R.size(), grid_side);
    grid_side = maximum_in_table(*rel_S, att_S.size(), grid_side);
    grid_side = maximum_in_table(*rel_T, att_T.size(), grid_side);
    grid_side = maximum_in_table(*rel_U, att_U.size(), grid_side);

    grid_side++;

    //cout << "Grid side: " << grid_side << endl;

    qdag_dfuds qdag_dfuds_rel_R(*rel_R, att_R, grid_side, 2, att_R.size());
    qdag_dfuds qdag_dfuds_rel_S(*rel_S, att_S, grid_side, 2, att_S.size());
    qdag_dfuds qdag_dfuds_rel_T(*rel_T, att_T, grid_side, 2, att_T.size());
    qdag_dfuds qdag_dfuds_rel_U(*rel_U, att_U, grid_side, 2, att_U.size());

    vector<qdag_dfuds> Q_dfuds(4);

    Q_dfuds[0] = qdag_dfuds_rel_R;
    Q_dfuds[1] = qdag_dfuds_rel_S;
    Q_dfuds[2] = qdag_dfuds_rel_T;
    Q_dfuds[3] = qdag_dfuds_rel_U;

	vector<uint256_t> results_partial_dfuds_back;

    uint8_t type_fun = argv[5] ? atoi(argv[5]) : 1;
    // size queue
    int64_t size_queue = argv[6] ? atoi(argv[6]) : 1000;


    high_resolution_clock::time_point start, stop;
    double total_time = 0.0;
    duration<double> time_span;

//    multiJoinPartialResultsDfudsBacktracking(Q_dfuds, grid_side, type_fun, size_queue, results_partial_dfuds_back);
    results_partial_dfuds_back.clear();

	uint256_t nodes_visited = 0;
    start = high_resolution_clock::now();

//    multiJoinPartialResultsDfudsBacktracking(Q_dfuds, grid_side, type_fun, size_queue, results_partial_dfuds_back, nodes_visited);
	multiJoinPartialResultsDfudsBacktracking(Q_dfuds, grid_side, type_fun, size_queue, results_partial_dfuds_back);

    stop = high_resolution_clock::now();
    time_span = duration_cast<microseconds>(stop - start);
    total_time = time_span.count();

    cout << /*"Multiway Join ended in " <<*/ total_time /*<< " seconds"*/ << endl;

    return 0;
}
