#include <iostream>
#include <vector>
#include <random>
#include "Flood.h"
#include <chrono>

int main(int argc, const char* argv[])
{
    int dimensions = 3;
    int records = 1e7;
    int grid_cell_size = 16384;
    int dataset_seed = 55;

    /* Printing dataset metadata */
    cout << "num dimensions = " << dimensions << endl;
    cout << "num records = " << records << endl;
    cout << "grid cell size (bytes) = " << grid_cell_size << endl;
    cout << endl;

    /* column-oriented data layout */
    vector<uint32_t> dataset;
    try {
        dataset.resize(dimensions * records);
    }
    catch (exception e)
    {
        cout << e.what();
    }

    /* Generating random data from uniform distribution. */
    using key_type = uint32_t;
    std::mt19937 gen(dataset_seed);
    uniform_int_distribution<key_type> key_distrib(0, 4e9);
    auto rand = [&gen, &key_distrib] { return key_distrib(gen); };
    std::generate(dataset.begin(), dataset.end(), rand);

    /* 4 Queries with different selectivities */
    vector<uint32_t> query_matches(4);
    std::vector<vector<uint32_t>> queriesMD;
    /* Query 1 */
    queriesMD.push_back({ 12341234,3213,1234142 });
    queriesMD.push_back({ 100000000,500000000,900000000 });
    /* Query 2 */
    queriesMD.push_back({ 100000,100000000,0 });
    queriesMD.push_back({ 3000000000,500000000,4000000000 });
    /* Query 3 */
    queriesMD.push_back({ 50000,400000000,0 });
    queriesMD.push_back({ 1000000000,900000000,4000000000 });
    /* Query 4 */
    queriesMD.push_back({ 70000000,98000000,0 });
    queriesMD.push_back({ 2000000000,120000000,4000000000 });

    /* visitor object result structure example */
    // attribute 1: 1, 2, 3, 394, 10
    // attribute 2: 8, 100, 23, 123, 20
    // attribute 3: 4, 403, 22, 99, 1, 4
    vector<vector<uint32_t>> visitor(dimensions);
    for (int i = 0; i < visitor.size(); i++)
        visitor[i].resize(records);
        
    flood::GridFile flood(dimensions, records, grid_cell_size, dataset, queriesMD, true, 0);

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0, c=1; i < queriesMD.size(); i+=2, c++)
    {
        query_matches[c-1] = flood.range_query_recursive_with_state_main(dimensions, queriesMD[i], queriesMD[i + 1], visitor);

    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;

    for (int i = 0, c = 1; i < queriesMD.size(); i += 2, c++)
    {
        cout << "query #" << c << " range query matches = " << query_matches[c-1] << endl;

    }
    std::cout << "Range query elapsed time = "
        << std::chrono::duration_cast<std::chrono::microseconds>(duration).count()
        << " microseconds" << std::endl;
}
