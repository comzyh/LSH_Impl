#include <cstdio>
#include <string>
#include <fstream>
#include <unordered_set>
#include <chrono>   
#include "vendor/cmdline/cmdline.h"  // argparser
#include "lsh_impl/lsh_impl.hpp"
#include "vendor/flann/src/cpp/flann/flann.hpp"

using namespace std;
using namespace std::chrono;
using std::chrono::system_clock;


int main(int argc, char **argv)
{
    // CommandLine parser
    cmdline::parser parser;
    parser.add<string>("base", 'b', "base data path.", true);
    parser.add<string>("query", 'q', "query data path.", true);
    parser.add<string>("ground", 'g', "ground_Truth data path.", true);
    parser.add<int>("table_num", 't', "number of table", false, 16);
    parser.add<int>("function_num", 'f', "number of function", false, 3);
    parser.add<float>("W", 'w', "W", false, 100.0f);
    parser.add<int>("probe", 'p', "number of probe", false, 200);
    parser.add<int>("limit", 'l', "limit the query number", false, -1);
    parser.add<int>("kdtree", 'k', "using kdtree in flann", false, 0);
    parser.add<int>("checks", 'c', "checks using by kd-tree", false, 128);

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

     //find
    int nn = ground_truth.col;
    lsh_impl::Matrix<int> indices(query.row, nn);
    lsh_impl::Matrix<float> dists(query.row, nn);


    //limit the count of query for debug
    int query_limit = parser.get<int>("limit");
    if (query_limit != -1) {
        query.row = query_limit;
    }

    system_clock::time_point strat_time, build_index_finish_time, finish_time;
    strat_time = system_clock::now();

    if (parser.get<int>("kdtree")) {
        printf("Using KDTree, \n");
        flann::Matrix<float> flann_base(base.data, base.row, base.col);
        flann::Matrix<float> flann_query(query.data, query.row, query.col);
        flann::Matrix<int> flann_indices(indices.data, indices.row, indices.col);
        flann::Matrix<float> flann_dists(dists.data, dists.row, dists.col);


        //build Index
        if (parser.get<int>("kdtree") < 0) { //brute-force
            flann::Index<flann::L2<float> > index(flann_base, flann::LinearIndexParams()); 
            printf("Building index.\n");
            index.buildIndex();
            build_index_finish_time = system_clock::now();
            printf("Finsih indexing: %.4lf s\n", duration_cast<milliseconds>(build_index_finish_time - strat_time).count()/ 1000.0);

            //Search
            printf("Searching.\n");
            index.knnSearch(flann_query, flann_indices, flann_dists, nn, flann::SearchParams(parser.get<int>("checks")));
            finish_time = system_clock::now();
            printf("Finsih Searching: %.4lf s\n", duration_cast<milliseconds>(finish_time - build_index_finish_time).count()/ 1000.0);


        } else {
            flann::Index<flann::L2<float> > index(flann_base, flann::KDTreeIndexParams(parser.get<int>("kdtree"))); 
            printf("Building index.\n");
            index.buildIndex();
            build_index_finish_time = system_clock::now();
            printf("Finsih indexing: %.4lf s\n", duration_cast<milliseconds>(build_index_finish_time - strat_time).count()/ 1000.0);

            //Search
            printf("Searching.\n");
            index.knnSearch(flann_query, flann_indices, flann_dists, nn, flann::SearchParams(parser.get<int>("checks")));
            finish_time = system_clock::now();
            printf("Finsih Searching: %.4lf s\n", duration_cast<milliseconds>(finish_time - build_index_finish_time).count()/ 1000.0);

        }
    }  
    else {
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
        build_index_finish_time = system_clock::now();
        printf("Finsih indexing: %.4lf s\n", duration_cast<milliseconds>(build_index_finish_time - strat_time).count()/ 1000.0);

        index.knnSearch(query, indices, dists, nn, index_params);
        finish_time = system_clock::now();
        printf("Finsih Searching: %.4lf s\n", duration_cast<milliseconds>(finish_time - build_index_finish_time).count()/ 1000.0);
    }    
    double recall = 0;
    for (size_t i = 0; i < query.row; i++) {
        lsh_impl::Vector<float> q(query.col, query[i]);
        unordered_set<int> truth;
        //GroundTruth
        for (size_t j = 0; j < nn; j++) {
            truth.insert(ground_truth[i][j]);
        }
        //ANN
        for (size_t j = 0; j < nn; j++) {
            recall += truth.count(indices[i][j]);
        }
    }
    recall /= query.row * nn;
    printf("Recall: %lf\n", recall);
    return 0;
}
