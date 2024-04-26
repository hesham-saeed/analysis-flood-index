#!/bin/bash

EXPERIMENT='index_comparison'
OUTPUT_FILE='index-comparison.csv'

if ! test -f ${OUTPUT_FILE}; then
    echo "benchmark_name,index_name,avg_query_time,avg_query_matches" > ${OUTPUT_FILE}
fi

./build/experiments/index_comparison_skewed_data_skewed_query
./build/experiments/index_comparison_skewed_data_uniform_query
./build/experiments/index_comparison_uniform_data_skewed_query
./build/experiments/index_comparison_uniform_data_uniform_query