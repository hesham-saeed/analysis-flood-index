#include <vector>
#include <iostream>
#include <random>
#include <stdint.h>
using namespace std;

class QueryGenerator
{
public:
	bool has(vector<uint32_t> query_sample, int idx)
	{
		for (int i = 0; i < query_sample.size(); i++)
			if (query_sample[i] == idx)
				return true;

		return false;
	}
void generateQueries(
    int qdims,
	const vector<uint32_t>& dataset,
	const vector<double>& column_selectivities, //as in: the higher the priority the higher selectivity an attribute is 
	const vector<uint32_t>& column_weights, //as in: the higher the priority the higher selectivity an attribute is 
	std::vector<std::vector<uint32_t>>& output_queries,
	uint32_t seed, int min_filters, int max_filters)
{
	cout << "query generating "<< output_queries.size()/2 << " queries at seed = " << seed <<"...."<< endl;
	std::random_device rd;
	//std::mt19937 gen(rd());
	std::mt19937 gen(seed);
	std::discrete_distribution<uint32_t> discrete_distrib(column_weights.begin(), column_weights.end());

	/* Changing Dataset Layout from Row to Column Layout to sort each attribute. */
	std::vector<vector<uint32_t>> column_layout_dataset(qdims);
    for (int i=0;i<qdims;i++)
        column_layout_dataset[i].resize(dataset.size()/qdims);
	for (int i = 0; i < dataset.size()/qdims; i++) // first row
		for (int j = 0; j < qdims; j++) // first column
            column_layout_dataset[j][i] = (dataset[j * dataset.size()/qdims + i]);

    //cout << "sorting" << endl;
	/* Sorting each column */
	for (int i = 0; i < qdims; i++)
		std::sort(column_layout_dataset[i].begin(), column_layout_dataset[i].end());

    //cout << "sampling" << endl;
	vector<vector<uint32_t>> query_samples(output_queries.size());
	std::mt19937 gen1(seed);
	uniform_int_distribution<uint32_t> num_attributes_distrib(min_filters, max_filters);
	for (int i = 0; i < output_queries.size(); i++)
	{
		uint32_t number_of_attributes = num_attributes_distrib(gen1);
		//uint32_t number_of_attributes = qdims;
		vector<uint32_t> query_sample(number_of_attributes, 0);
		for (int j = 0; j < number_of_attributes; j++)
		{
			uint32_t attribute = discrete_distrib(gen);
			query_sample[j] = attribute;
		}
		query_samples[i] = query_sample;
	}

    //cout << "grouping" << endl;
	std::mt19937 gen2(seed);
	double largest_selectivity = *(std::max_element(column_selectivities.begin(), column_selectivities.end()));
	uniform_int_distribution<uint32_t> dataset_distrib(0, (dataset.size()/qdims - 1) * (1 - largest_selectivity));
	//auto rand = [&gen2, &key_distrib] { return key_distrib(gen2); };

	// producing queries
    //cout << "producing"<<endl;
	for (int i = 0; i < output_queries.size() - 1; i += 2)
	{
		vector<uint32_t> query_sample_start_interval(qdims);
		vector<uint32_t> query_sample_end_interval(qdims);
        uint32_t rand_start_index = dataset_distrib(gen2);
		for (int j = 0; j < query_sample_start_interval.size(); j++)
		{
			if (has(query_samples[i], j))
			{
				uint32_t random_value = column_layout_dataset[j][rand_start_index];
				query_sample_start_interval[j] = column_layout_dataset[j][rand_start_index];

				/* ending interval will depend on column selectivity */
				long long rand_end_index = rand_start_index + column_selectivities[j] * (dataset.size()/qdims);
				query_sample_end_interval[j] = column_layout_dataset[j][rand_end_index];
				//cout << rand_start_index << "\t" << rand_end_index << endl;
			}
			else {
				query_sample_start_interval[j] = column_layout_dataset[j][0];
				query_sample_end_interval[j] = column_layout_dataset[j][dataset.size()/qdims - 1];
			}
		}

		output_queries[i] = query_sample_start_interval;
		output_queries[i + 1] = query_sample_end_interval;
	}

	cout << "queries generated.\n";
}
};
