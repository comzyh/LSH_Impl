#include <cstdio>
#include <string>
#include <fstream>
#include <unordered_set>
#include "vendor/cmdline/cmdline.h"  // argparser
#include "lsh_impl/lsh_impl.hpp"

using namespace std;

int main(int argc, char **argv)
{
    // CommandLine parser
    cmdline::parser parser;
    parser.add<string>("base", 'b', "Base data path.", true);
    parser.add<string>("query", 'q', "Query data path.", true);
    parser.add<string>("ground", 'g', "Ground_Truth data path.", true);
    parser.add<int>("table_num", 't', "Number of table", false, 16);
    parser.add<int>("function_num", 'f', "Number of function", false, 3);
    parser.add<float>("W", 'w', "W", false, 100.0f);
    parser.add<int>("probe", 'p', "number of probe", false, 200);
    parser.add<int>("limit", 'l', "limit the query number", false, -1);
    parser.parse_check(argc, argv);

    //LoadData
    //http://corpus-texmex.irisa.fr/
    size_t filesize, dimension;

    lsh_impl::Matrix<float> base = lsh_impl::Matrix<float>::load_from_file(parser.get<string>("base").c_str());
    lsh_impl::Matrix<float> query = lsh_impl::Matrix<float>::load_from_file(parser.get<string>("query").c_str());
    lsh_impl::Matrix<int> ground_truth = lsh_impl::Matrix<int>::load_from_file(parser.get<string>("ground").c_str());

    printf("Base        : %8lu row, %4lu col\n", base.row, base.col);
    printf("Query       : %8lu row, %4lu col\n", query.row, query.col);
    printf("Ground truth: %8lu row, %4lu col\n", ground_truth.row, ground_truth.col);

    //LshIndexParams
    //size_t table_num, size_t function_num, float W, size_t probe_num
    lsh_impl::LshIndexParams index_params(
        parser.get<int>("table_num"),
        parser.get<int>("function_num"),
        parser.get<float>("W"),
        parser.get<int>("probe")
    );
    lsh_impl::LSH_Index<float> index(base, index_params);
    index.buildIndex();

    //find
    int nn = ground_truth.col;
    lsh_impl::Matrix<int> indices(query.row, nn);
    lsh_impl::Matrix<float> dists(query.row, nn);

    // lsh_impl::Matrix<float> fake_query;
    // fake_query.row = 100;
    // fake_query.col = base.col;
    // fake_query.data = base.data;


    //limit the count of query for debug
    int query_limit = parser.get<int>("limit");
    if (query_limit != -1) {
        query.row = query_limit;
    }

    index.knnSearch(query, indices, dists, nn, index_params);

    double recall = 0;
    for (size_t i = 0; i < query.row; i++) {
        // printf("query %4d:\n", i);
        lsh_impl::Vector<float> q(query.col, query[i]);
        unordered_set<int> truth;
        //GroundTruth
        for (size_t j = 0; j < nn; j++) {
            // printf("(%d, %.0f) ", ground_truth[i][j], q.square_dist(base[ground_truth[i][j]]));
            truth.insert(ground_truth[i][j]);
        }
        // printf("\n\n");
        //Mine
        for (size_t j = 0; j < nn; j++) {
            // printf("(%d, %.0f) ", indices[i][j], dists[i][j]);
            recall += truth.count(indices[i][j]);
        }
        // printf("\n"); 
    }
    recall /= query.row * nn;
    printf("Recall: %lf\n", recall);
    return 0;
}
