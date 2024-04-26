#include <iostream>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/iterator.hpp>
#include <boost/iterator/function_output_iterator.hpp>
#include <vector>
#include <stdint.h>
#include <random>
#include "Flood.h"
#include "GridFile.h"
#include "query_generator.h"
#include "utils.h"
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<uint32_t, 6, bg::cs::cartesian> multiDimensionalPoint;
using Tree = bgi::rtree<multiDimensionalPoint, bgi::rstar<50>>; 
using namespace std;


int main()
{
    int dimensions = 6;
    int records = 1e8;
    int grid_cell_size = 65536;
    int seed = 55;
    int queries = 1000;

    /******************************************************************************/
    /* Dataset Generation */
    /******************************************************************************/
    vector<uint32_t> dataset(records*dimensions);
    using key_type = uint32_t;
    std::mt19937 gen(seed);
    uniform_int_distribution<key_type> key_distrib(0, 4e9);
    auto rand_uniform = [&gen, &key_distrib] { return key_distrib(gen); };
    std::generate(dataset.begin(), dataset.begin()+(records*dimensions/2), rand_uniform);

    double mean = 1e9;
    double stddev = 5e7;
    std::normal_distribution<double> distribution;
    distribution.param(std::normal_distribution<double>::param_type(mean, stddev));
    auto rand_normal = [&gen, &distribution] { return distribution(gen); };
	std::generate(dataset.begin()+(records*dimensions/2), dataset.begin()+(records*4), rand_normal);
    distribution.param(std::normal_distribution<double>::param_type(5e8, 3e7));
    auto rand_normal2 = [&gen, &distribution] { return distribution(gen); };
	std::generate(dataset.begin()+(records*4), dataset.begin()+(records*5), rand_normal2);
    distribution.param(std::normal_distribution<double>::param_type(2e9, 6e7));
    auto rand_normal3 = [&gen, &distribution] { return distribution(gen); };
	std::generate(dataset.begin()+(records*5), dataset.end(), rand_normal3);


    /******************************************************************************/
    /* Query Generation */
    /******************************************************************************/
    QueryGenerator qg;
    vector<double> col_selectivities(dimensions);
    vector<uint32_t> col_weights(dimensions);
    for (int i=0;i<col_selectivities.size();i++)
        col_selectivities[i] = 0.13; //0.13 6 filters.
    for (int i=0;i<col_weights.size();i++)
        col_weights[i] = 1;
    std::vector<vector<uint32_t>> queriesMD(queries*2);
    qg.generateQueries(dimensions, dataset, col_selectivities,col_weights,queriesMD,seed,dimensions,dimensions);
    vector<uint64_t> query_matches(queriesMD.size()/2);


    /******************************************************************************/
    /* R-Tree */
    /******************************************************************************/
    {
        Tree* rt;
    
        multiDimensionalPoint r1;
        multiDimensionalPoint r2;
        vector<multiDimensionalPoint> input(records);
        for (int i = 0; i < records; i++)
        {
            multiDimensionalPoint p;
            p.set<0>(dataset[0 * records + i]);
            p.set<1>(dataset[1 * records + i]);
            p.set<2>(dataset[2 * records + i]);
            p.set<3>(dataset[3 * records + i]);
            p.set<4>(dataset[4 * records + i]);
            p.set<5>(dataset[5 * records + i]);

            input[i] = p;
        }
        try{
            rt = new Tree(input);
        }
        catch(exception e)
        {
            cout << e.what()<<endl;
        }
        
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0, c=1; i < queriesMD.size(); i+=2, c++)
        {
            vector<multiDimensionalPoint> result2;
            r1.set<0>(queriesMD[i][0]);
            r1.set<1>(queriesMD[i][1]);
            r1.set<2>(queriesMD[i][2]);
            r1.set<3>(queriesMD[i][3]);
            r1.set<4>(queriesMD[i][4]);
            r1.set<5>(queriesMD[i][5]);

            r2.set<0>(queriesMD[i+1][0]);
            r2.set<1>(queriesMD[i+1][1]);
            r2.set<2>(queriesMD[i+1][2]);
            r2.set<3>(queriesMD[i+1][3]);
            r2.set<4>(queriesMD[i+1][4]);
            r2.set<5>(queriesMD[i+1][5]);

            bg::model::box<multiDimensionalPoint> mybox(r1, r2);
            rt->query(bgi::intersects(mybox), std::back_inserter(result2));
            query_matches[c-1] = result2.size();
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = end - start;
        uint64_t total_matches = std::accumulate(query_matches.begin(), query_matches.end(),0ULL);
        double avg_query_time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()/queries;
        double avg_query_matches = total_matches/queries;
        cout << "avg query matches = " << avg_query_matches <<endl;
        std::cout << "Average Range query time = " << avg_query_time << " millisecond" << std::endl;
        writeBenchmark("index-comparison.csv", "skewed_data_skewed_query", "R-Tree", avg_query_time, avg_query_matches);
    }

    /******************************************************************************/
    /* Grid File */
    /******************************************************************************/
    {
        vector<vector<uint32_t>> visitor(dimensions);
        for (int i = 0; i < visitor.size(); i++)
            visitor[i].resize(records);

        gf::GridFile gridfile(dimensions, records, grid_cell_size, dataset);

        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0, c=1; i < queriesMD.size(); i+=2, c++)
        {
            query_matches[c-1] = gridfile.range_query(dimensions, queriesMD[i], queriesMD[i + 1], visitor);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = end - start;
        uint64_t total_matches = std::accumulate(query_matches.begin(), query_matches.end(),0ULL);
        double avg_query_time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()/queries;
        double avg_query_matches = total_matches/queries;
        cout << "avg query matches = " << avg_query_matches <<endl;
        std::cout << "Average Range query time = " << avg_query_time << " millisecond" << std::endl;
        writeBenchmark("index-comparison.csv", "skewed_data_skewed_query", "GridFile", avg_query_time, avg_query_matches);
    }


    /******************************************************************************/
    /* Flood */
    /******************************************************************************/
    {
        vector<vector<uint32_t>> visitor(dimensions);
        for (int i = 0; i < visitor.size(); i++)
            visitor[i].resize(records);
            
        flood::GridFile flood(dimensions, records, grid_cell_size, dataset, queriesMD);

        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0, c=1; i < queriesMD.size(); i+=2, c++)
        {
            query_matches[c-1] = flood.range_query(dimensions, queriesMD[i], queriesMD[i + 1], visitor);

        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = end - start;
        uint64_t total_matches = std::accumulate(query_matches.begin(), query_matches.end(),0ULL);
        double avg_query_time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()/queries;
        double avg_query_matches = total_matches/queries;
        cout << "avg query matches = " << avg_query_matches <<endl;
        std::cout << "Average Range query time = " << avg_query_time << " millisecond" << std::endl;
        writeBenchmark("index-comparison.csv", "skewed_data_skewed_query", "Flood", avg_query_time, avg_query_matches);
    }
    
    return 0;
}