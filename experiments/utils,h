#pragma once 
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

template<typename ...T>
static bool writeBenchmark(string filename, string benchmark_name, string index_name, T&... args) {
    try
    {
        stringstream ss;
        ofstream file;
        file.open(filename, ios_base::app);
        ss << benchmark_name << "," << index_name;
        for (const auto& arg : { args... })
            ss << "," << arg;
        ss << endl;
        file << ss.str();
        file.close();
    }
    catch (std::ifstream::failure& readErr)
    {
        cout << readErr.what();
        return false;
    }
    return true;
}