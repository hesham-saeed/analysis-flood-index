#include <iostream>
#include "../GridFile.h"
#include <vector>
#include <random>
#include <chrono>

int main(int argc, const char* argv[])
{
    int dimensions = 3;
    int records = 1e8;
    int grid_cell_size = 16384;
    
    try{
            dimensions = stoi(argv[1]);
            records = stoi(argv[2]);
            grid_cell_size = stoi(argv[3]);
    }
    catch (exception e){
        //cout << e.what() << endl;
        cout <<"Some args are missing or failed to be parsed. Defaults will be used instead.\n\n";
    }


    /* Printing dataset metadata */
    cout << "num dimensions = " << dimensions << endl;
    cout << "num records = " << records << endl;
    cout << "grid cell size (bytes) = " << grid_cell_size << endl;
    cout << endl;

    /* column oriented dataset. a new column begins after every set of (records). */
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
    std::mt19937 gen(55555);
    uniform_int_distribution<key_type> key_distrib(0, 4e9);
    auto rand = [&gen, &key_distrib] { return key_distrib(gen); };
    std::generate(dataset.begin(), dataset.end(), rand);


    /* Visitor object to store the range query results. */
    vector<vector<uint32_t>> visitor(dimensions);
    for (int i = 0; i < visitor.size(); i++)
    {
        /* visitor object result structure example */
        // attribute 1: 1, 2, 3, 394, 10
        // attribute 2: 8, 100, 23, 123, 20
        // attribute 3: 4, 403, 22, 99, 1, 4
        visitor[i].resize(records);
    }
    gf::GridFile gridfile(dimensions, records, grid_cell_size, dataset);
    dataset.clear();
    dataset.shrink_to_fit();

    auto start = std::chrono::high_resolution_clock::now();
    uint64_t res = gridfile.range_query(dimensions, { 100000000,100000000, 100000000 }, { 300000000,300000000, 300000000 }, visitor);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = end - start;


    for (int i=0;i<visitor.size();i++)
        visitor[i].resize(res);
    //cout <<"visitor size = " << visitor[0].size() << endl;
    cout << "\nRange query found = " << res  << " matches"<< endl << endl;
    std::cout << "Range query time elapsed = "
        << std::chrono::duration_cast<std::chrono::microseconds>(duration).count()
        << " microseconds" << std::endl;

}
