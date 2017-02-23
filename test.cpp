#include <cstdio>
#include <string>
#include <fstream>
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
    parser.add<int>("table_num", 't', "Number of table", false, 10);
    parser.add<int>("function_num", 'f', "Number of function", false, 1);
    parser.add<float>("W", 'w', "W", false, 10.0f);
    parser.add<int>("probe", 'p', "number of probe", false, 200);
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
    index.knnSearch(query, indices, dists, nn, index_params);

    for (size_t i = 0; i < 1; i++) {
        printf("query %4d:\n", i);
        for (size_t j = 0; j < nn; j++) {
            printf("%6d ", indices[i][j]);
        }
        printf("\n");
        for (size_t j = 0; j < nn; j++) {
            printf("%6d ", ground_truth[i][j]);
        }
        printf("\n");
    }
    return 0;
}
