#pragma once

#include <iostream>
#include <vector>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <math.h>
#include <chrono>
using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

namespace flood {
class GridFile {
private:

	struct GridCell
	{
		//static uint32_t count;
		uint32_t** cellPoints2;
		uint32_t cell_number;
		uint32_t size;
		uint32_t idx;
		uint32_t dimensions;


		void setCellNumber(uint32_t cell_no) {
			cell_number = cell_no;
		}


		GridCell() {
			size = 16;
			cell_number = idx = 0;
			cellPoints2 = NULL;
			dimensions = 1;
		}

		~GridCell()
		{
			for (size_t i = 0; i < dimensions; ++i) {
				delete[] cellPoints2[i];
			}
			delete[] cellPoints2;
		}

		void initializeCells2(int dimensions, int initial_size) {
			size = initial_size;
			idx = 0;
			this->dimensions = dimensions;
			cellPoints2 = new uint32_t * [dimensions];
			for (int i = 0; i < dimensions; i++)
				cellPoints2[i] = new uint32_t[size];
		}

		void addPoint2(uint32_t* p, int dimensions) {
			if (idx == size)
			{
				size *= 2;
				uint32_t** tmp = new uint32_t * [dimensions];
				for (int i = 0; i < dimensions; i++)
				{
					tmp[i] = new uint32_t[size];
				}

				for (int i = 0; i < dimensions; i++)
				{
					for (int j = 0; j < idx; j++)
					{
						tmp[i][j] = cellPoints2[i][j];
					}
				}

				for (int i = 0; i < dimensions; ++i) {
					delete[] cellPoints2[i];
				}
				delete[]cellPoints2;
				cellPoints2 = tmp;
			}
			for (int i = 0; i < dimensions; i++)
			{
				cellPoints2[i][idx] = p[i];
			}
			idx++;

		}

		void addPoint(uint32_t* p) {
			if (idx == size)
			{
				size *= 2;
				uint32_t** tmp = new uint32_t * [size];
				for (int i = 0; i < idx; i++)
				{
					tmp[i] = cellPoints2[i];
				}
				delete[]cellPoints2;
				cellPoints2 = tmp;
			}
			cellPoints2[idx++] = p;
		}

		void shrink_cell(int& num_dimensions) {
			//uint32_t** tmp = new uint32_t * [idx];
			uint32_t** tmp = new uint32_t * [num_dimensions];
			for (int i = 0; i < num_dimensions; i++)
				tmp[i] = new uint32_t[idx];
			for (int i = 0; i < num_dimensions; i++)
			{
				for (int j = 0; j < idx; j++)
				{
					tmp[i][j] = cellPoints2[i][j];
				}
			}
			size = idx;

			for (size_t i = 0; i < num_dimensions; ++i) {
				delete[] cellPoints2[i];
			}
			delete[] cellPoints2;

			cellPoints2 = tmp;
		}



	};
private:
	int num_records;
	int insertedPoints;
	int num_dimensions;
	int min_num_grid_cells;
	int min_points_per_cell;
	int q;
	int gridCellsCount;
	vector<uint32_t> mulfactorArr;
	std::vector<uint32_t> rangeStartPos;
	std::vector<uint32_t> rangeEndPos;
	std::vector<uint32_t> dimsSize; //example : a grid of 2x6x3x8;
	vector<vector<pair<uint32_t, uint32_t>>> gridCellLimits; // size = number of created cells, to be freed as soon as the grid initialization is finalized.
	vector<vector<uint32_t>> attributeSlicers; // [0,10,20,30], [50,100], [100,200,300]
	vector<vector<uint32_t>> attributeSlicersPointers; // [0,10,20,30], [50,100], [100,200,300]
	GridCell* gridCells;
	vector<uint32_t> minValueForEachDimension;
	vector<uint32_t> maxValueForEachDimension;
	vector<uint32_t> orderDims;

	std::chrono::steady_clock::time_point total_scan_time;
	std::chrono::steady_clock::time_point total_refinement_time;
	uint64_t q_num_scanned_cells;
	uint64_t q_num_scanned_points;

public:
	/* Sorting grid cell */
	void sortGridCells(int dimension = 0) {
		for (int i = 0; i < gridCellsCount; i++) {
			vector<uint32_t> indexes(gridCells[i].size);
			for (int k = 0; k < gridCells[i].size; k++)
				indexes[k] = k;

			sort(indexes.begin(), indexes.end(), [&](uint32_t& a, uint32_t& b) {
				return gridCells[i].cellPoints2[0][a] < gridCells[i].cellPoints2[0][b];
				});

			uint32_t** cellPoints2Tmp = new uint32_t * [num_dimensions];
			for (int j = 0; j < num_dimensions; j++)
				cellPoints2Tmp[j] = new uint32_t[gridCells[i].size];

			for (int k = 0; k < gridCells[i].size; k++)
			{
				for (int j = 0; j < num_dimensions; j++)
				{
					cellPoints2Tmp[j][k] = gridCells[i].cellPoints2[j][indexes[k]];
				}
			}
			std::swap(gridCells[i].cellPoints2, cellPoints2Tmp);
			for (size_t i = 0; i < num_dimensions; ++i) {
				delete[] cellPoints2Tmp[i];
			}
			delete[] cellPoints2Tmp;
		}
		


	}

	void reOrderDataset(vector<uint32_t>& datasetMD)
	{
		vector<uint32_t> datasetMD2(datasetMD);
		for (int i = 0; i < num_records; i++)
		{
			//datasetMD2[i] = datasetMD2[orderDims[i]];
			for (int j = 0; j < num_dimensions; j++)
			{
				datasetMD2[j * num_records + i] = datasetMD[orderDims[j] * num_records + i];
			}
		}
		datasetMD = datasetMD2;
	}

	GridFile(uint32_t dimensions, int num_records, int CELL_SIZE,
		vector<uint32_t>& datasetMD,
		std::vector<std::vector<uint32_t>>& queriesMD, 
		bool sortedDimension = false, 
		int sortingDimension = 0) {

		this->num_records = num_records;
		this->num_dimensions = dimensions;
		this->insertedPoints = 0;
		cout << "*************************************" << endl;
		cout << "Constructing Flood Index..." << endl;
		cout << "*************************************" << endl;
		computeDimensionsLimits(datasetMD);
		dimsSize.resize(dimensions);
		orderDims.resize(dimensions);
		getDimSizes2(dimensions, CELL_SIZE, num_records, queriesMD);
		reOrderDataset(datasetMD);
		configureGrid(num_records);
		configureSlicers2(datasetMD);
		configureCellsLimits();
		assignGridNumber();

		cout << "Bulk loading the grid file with " << num_records << " records...." << endl;
		fillGrid(datasetMD);
		for (int i = 0; i < gridCellsCount; i++)
		{
			gridCells[i].shrink_cell(num_dimensions);
		}

		cout << "Sorting all grid cells..." << endl;
		if (sortedDimension)
			sortGridCells(sortingDimension);

		for (int i = 1; i <= dimensions; i++) {
			uint32_t mulfactor = 1;
			for (int j = 0; j < i - 1; j++) {
				mulfactor *= dimsSize[j];
			}
			mulfactorArr.push_back(mulfactor);
		}

		//uint32_t minPointsAtAllCells = 1000000;
		//uint32_t maxPointsAtAllCells = 0;
		//for (int i = 0; i < gridCellsCount; i++) {
		//	if (gridCells[i].size < minPointsAtAllCells)
		//		minPointsAtAllCells = gridCells[i].size;
		//	if (gridCells[i].size > maxPointsAtAllCells)
		//		maxPointsAtAllCells = gridCells[i].size;
		//}
		//cout << "Cell with lowest number of points has=" << minPointsAtAllCells << endl;
		//cout << "Cell with highest number of points has=" << maxPointsAtAllCells << endl;
		rangeStartPos.resize(num_dimensions);
		rangeEndPos.resize(num_dimensions);
		cout << "****************************************\n";
	}

	~GridFile()
	{
		delete[] gridCells;
	}
private:

	void computeDimensionsLimits(std::vector<uint32_t>& dataset) {
		minValueForEachDimension.resize(num_dimensions);
		maxValueForEachDimension.resize(num_dimensions);
		for (int i = 0; i < num_dimensions; i++)
		{
			minValueForEachDimension[i] = UINT32_MAX;
			maxValueForEachDimension[i] = 0;
		}

		for (int i = 0; i < num_records; i++)
		{
			for (int j = 0; j < num_dimensions; j++)
			{
				if (minValueForEachDimension[j] > dataset[j * num_records + i])
					minValueForEachDimension[j] = dataset[j * num_records + i];
				if (maxValueForEachDimension[j] < dataset[j * num_records + i])
					maxValueForEachDimension[j] = dataset[j * num_records + i];
			}
		}
	}

	std::vector<pair<double, uint32_t>> getFrequencyRatios(vector<vector<uint32_t>>& queriesMD)
	{
		std::vector<pair<double, uint32_t>> frequencies(num_dimensions);
		for (int i = 0; i < frequencies.size(); i++)
			frequencies[i].first = 0;
		for (int i = 0; i < num_dimensions; i++)
			frequencies[i].second = i;
		for (int i = 0; i < queriesMD.size(); i += 2)
		{
			for (int j = 0; j < num_dimensions; j++)
			{
				if (queriesMD[i][j] <= minValueForEachDimension[j] && queriesMD[i + 1][j] >= maxValueForEachDimension[j]) {
					//do not count as this attribute is not used
				}
				else
				{
					frequencies[j].first++;
				}
			}
		}
		return frequencies;
	}


	std::vector<pair<double, uint32_t>> getDimSizes2(
		uint32_t dimensions, 
		uint32_t PAGE_SIZE, 
		uint32_t dataset_size, 
		vector<vector<uint32_t>>& queriesMD) 
	{
		uint32_t maxElementPerPage = std::round((double)PAGE_SIZE / (dimensions * 4));
		if (maxElementPerPage == 0)
			maxElementPerPage++;
		//cout << "maxElementsPerPage = " << maxElementPerPage << endl;
		dimsSize.resize(dimensions);  //-1 because of excluding the sort dimension

		uint32_t expectedNumberOfGridCells = dataset_size / maxElementPerPage; //the actual number of grid cells will depend on partitionsPerDimensions
		//uint32_t expectedNumberOfGridCells = dataset_size / 4096; //the actual number of grid cells will depend on partitionsPerDimensions
		if (dataset_size % maxElementPerPage != 0)
			expectedNumberOfGridCells++;

		//int dimensions_prioritized_count = 0;
		//for (int i = 0; i < dimensions_priority.size() - 1; i++) { //-1 because of excluding the sort dimension
		//	if (dimensions_priority[i] != 0)
		//		dimensions_prioritized_count++;
		//}

		auto frequencies = getFrequencyRatios(queriesMD);
		for (int i = 0; i < frequencies.size(); i++) {
			if (frequencies[i].first == 0)
				frequencies[i].first = 1;
		}

		// Treat all dimensions equally.
		//for (int i = 0; i < frequencies.size(); i++) {
		//	frequencies[i].first = 1;
		//}
		bool asc = false;

		////Desccending Order
		//asc = false;
		//std::sort(frequencies.begin(), frequencies.end(), [](auto& left, auto& right) {
		//	return left.first > right.first;
		//	});
		//int min_num = frequencies[frequencies.size() - 1].first;

		////Ascending Order
		std::sort(frequencies.begin(), frequencies.end(), [](auto& left, auto& right) {
			return left.first < right.first;
			});
		int min_num = frequencies[0].first;
		asc = true;

		//double partitionsPerDimension = std::pow(expectedNumberOfGridCells, 1.0 / dimensions);
		//double partitionsPerDimension = std::pow(expectedNumberOfGridCells, 1.0 / dimensions_prioritized_count);
		//partitionsPerDimension = ceil(partitionsPerDimension);

		//find the higest frequency to make it a sort dim


		//for (int i = 1; i < frequencies.size(); i++)
		//	if (frequencies[i].first < min_num)
		//		min_num = frequencies[i].first;

		cout << "Attribute frequencies = ";
		for (int i = 0; i < frequencies.size(); i++)
			cout << frequencies[i].first << " ";
		cout << endl;

		//divide by the minimum frequency
		for (int i = 0; i < frequencies.size(); i++)
			frequencies[i].first = std::round(frequencies[i].first / min_num);

		//frequencies[frequencies.size() - 1].first
		//std::sort(frequencies.begin(), frequencies.end(), [](auto& left, auto& right) {
		//	return left.first < right.first;
		//});

		//for (int i = 0; i < frequencies.size(); i++)
		//{
		//	cout << frequencies[i].first << " ";
		//}

		double product_of_elems = 1;
		for (int i = 0; i < frequencies.size(); i++)
			product_of_elems *= frequencies[i].first;
		
		double partitionsPerDimension = std::pow(expectedNumberOfGridCells / product_of_elems, 1.0 / (num_dimensions - 1));
		cout << "parittionPerDimension = " << partitionsPerDimension << endl;
		dimsSize[0] = 1;

		if (asc)
			for (int i = 1; i < num_dimensions; i++)
				dimsSize[i] = std::round(partitionsPerDimension * frequencies[i].first);
		else
			for (int i = 1; i < num_dimensions; i++)
				dimsSize[i] = std::round(partitionsPerDimension * frequencies[i - 1].first);

		for (int i = 0; i < num_dimensions; i++) {
			orderDims[i] = frequencies[i].second;
		}

		cout << "Number of splits per dimension = ";
		for (int i = 0; i < num_dimensions; i++)
			cout << dimsSize[i] << " ";

		//cout << "\nOrder of dims = ";
		//for (int i = 0; i < num_dimensions; i++)
		//	cout << orderDims[i] << " ";

		cout << "\nExpected number of grid cells = " << expectedNumberOfGridCells << endl;
		min_num_grid_cells = expectedNumberOfGridCells;
		min_points_per_cell = maxElementPerPage;
		return frequencies;
	}


	string getDimensions() {
		string s = "Dimensional splits = ";
		for (int i = 0; i < dimsSize.size(); i++) {
			string tmp = to_string(dimsSize[i]) + " ";
			s += tmp;
		}
		return s;
	}

	void configureGrid(size_t dataset_size) {
		gridCellsCount = 1;
		for (int i = 0; i < dimsSize.size(); i++)
		{
			gridCellsCount *= dimsSize[i];
		}

		for (int i = 0; i < dimsSize.size(); i++)
		{
			if ((gridCellsCount / dimsSize[i] * (dimsSize[i] - 1)) > min_num_grid_cells)
			{
				gridCellsCount = (gridCellsCount / dimsSize[i] * (dimsSize[i] - 1));
				dimsSize[i] -= 1;
			}
			else
				break;
		}
		cout << "Actual number of grid cells = " << gridCellsCount << endl;
		//cout << getDimensions() << endl;

		gridCells = new GridCell[gridCellsCount];
		for (int i = 1; i <= gridCellsCount; i++)
		{
			gridCells[i - 1].setCellNumber(i);
			gridCells[i - 1].initializeCells2(this->num_dimensions, min_points_per_cell);
		}
		gridCellLimits.resize(gridCellsCount);
		size_t initialCellSize = dataset_size / gridCellsCount + ((dataset_size % gridCellsCount == 0) ? 0 : 1);


	}

	void configureSlicers2(vector<uint32_t>& dataset)
	{
		attributeSlicers.resize(num_dimensions); //-1 for the sorting dimension
		attributeSlicersPointers.resize(num_dimensions); //-1 for the sorting dimension

		for (int j = 0, k = 0; j < dimsSize.size(); j++) {
			if (dimsSize[j] == 0)
				continue;
			vector<uint32_t> tmp(num_records);
			for (int i = 0; i < num_records; i++) {
				//tmp[i] = dataset[i].p[j];
				tmp[i] = dataset[j * num_records + i];
			}
			/*cout << "Printing tmp in configureSlicers()" << endl;
			for (int i = 0; i < tmp.size(); i++) {
				cout << tmp[i] << endl;
			}*/
			sort(tmp.begin(), tmp.end());
			//cout << "min max = " << tmp[0] << " , " << tmp[tmp.size() - 1] << endl;
			uint32_t base = tmp.size() / dimsSize[j];
			//slicers[j]
			attributeSlicers[k].push_back(tmp[0]);
			//slicersPointer[j].push_back(0);
			for (uint32_t i = base, l = 1; l < dimsSize[j] && i < tmp.size(); i += base, l++) {
				if (attributeSlicers[k][attributeSlicers[k].size() - 1] >= tmp[i])
				{
					//handling duplicates
					int x = i + 1;
					while ((x < tmp.size()) && attributeSlicers[k][attributeSlicers[k].size() - 1] >= tmp[x])
						x++;
					if (x < tmp.size())
					{
						attributeSlicers[k].push_back(tmp[x]);
						attributeSlicersPointers[k].push_back(0);
					}
					else
					{
						attributeSlicers[k].push_back(tmp[x - 1] + 1);
						attributeSlicersPointers[k].push_back(0);
					}
					//i = x;
				}
				else
				{
					attributeSlicers[k].push_back(tmp[i]);
					attributeSlicersPointers[k].push_back(0);
				}
			}
			attributeSlicers[k].push_back(tmp[tmp.size() - 1]);
			attributeSlicersPointers[k].push_back(0);
			k++;
		}
	}

	void configureCellsLimits() {
		//TODO divide domain and set limits
		uint32_t mulfactor = 1;
		for (int j = 0, l = 0; j < dimsSize.size(); j++) {
			if (j > 0)
				mulfactor *= dimsSize[j - 1];
			uint32_t k = 0;
			for (int i = 0; i < gridCellsCount; i++) {

				if (i % mulfactor == 0)
					k++;

				if (k == attributeSlicers[l].size())
					k = 1;

				if (k == 1)
					gridCellLimits[i].push_back(make_pair(attributeSlicers[j][k - 1], attributeSlicers[j][k]));
				else
					gridCellLimits[i].push_back(make_pair(attributeSlicers[j][k - 1] + 1, attributeSlicers[j][k]));

				//if (k == mulfactor)
				//	k++;
				//if (k == slicers[j].size())
				//	k = 1;

			}
			l++;
		}
		//cout << "\n----- cell limits ---------\n";
		//for (int i = 0; i < gridCellsCount; i++)
		//{
		//	cout << "cell " << i;
		//	for (int j = 0; j < num_dimensions; j++)
		//	{
		//		cout << "= {" << gridCellLimits[i][j].first << ", " << gridCellLimits[i][j].second << "}\t";
		//	}
		//	cout << endl;
		//}
	}

	void assignGridNumber() {
		for (int i = 0; i < gridCellsCount; i++) {

			for (int j = 0; j < attributeSlicers.size(); j++) //iterate over slicer for every dimension
			{
				for (int k = 0; k < attributeSlicers[j].size(); k++) {
					if (k == 0) {
						if ((attributeSlicers[j][k] == (gridCellLimits[i][j].first)) && attributeSlicersPointers[j][k] < 1) {
							//if ((slicers[j][k] == (gridCells[i].cellLimits[j].first)) && slicersPointer[j][k] < 1) {
							attributeSlicersPointers[j][k] = gridCells[i].cell_number;
						}
					}
					else
					{
						if (((attributeSlicers[j][k] + 1) == (gridCellLimits[i][j].first)) && attributeSlicersPointers[j][k] < 1) {
							//if (((slicers[j][k] + 1) == (gridCells[i].cellLimits[j].first)) && slicersPointer[j][k] < 1) {
							attributeSlicersPointers[j][k] = gridCells[i].cell_number;
						}
					}

				}

			}

		}

	}

	void fillGrid(vector<uint32_t>& dataset) {
		int total_recs = 0;

		uint32_t* p = new uint32_t[num_dimensions]; //placeholder
		//uint32_t* ordered_point = new uint32_t[num_dimensions];
		for (int i = 0; i < num_records; i++)
		{
			total_recs++;

			//cout << "Inserting ";
			for (int j = 0; j < num_dimensions; j++)
			{
				p[j] = dataset[j * num_records + i];
			}
			//for (int j = 0; j < num_dimensions; j++)
			//{
			//	ordered_point[j] = p[orderDims[j]];
			//}
			insertPoint2(dimsSize.size(), p, 0);
		}
		//delete[]p; //freeing placeholder
		//delete[]ordered_point;
		cout << "DONE! Bulk loaded points = " << insertedPoints << endl;
	}

	void insertPoint2(uint32_t dim, uint32_t* point, uint32_t start) {
		uint32_t mulfactor = 1;
		for (int i = 0; i < dim - 1; i++) {
			mulfactor *= dimsSize[i];
		}

		auto ptr = std::lower_bound(attributeSlicers[dim - 1].begin(), attributeSlicers[dim - 1].end(), point[dim - 1]);
		if (ptr == attributeSlicers[dim - 1].end())
		{
			ptr--;
		}
		uint32_t tmpIndex = ptr - attributeSlicers[dim - 1].begin();

		if (tmpIndex != 0)
			tmpIndex--;
		uint32_t tmp = attributeSlicersPointers[dim - 1][tmpIndex];

		start += tmp - 1;
		if (dim == 1) {
			bool inLimits = true;
			for (int i = 0; i < gridCellLimits[start].size(); i++) {
				if (!(point[i] >= gridCellLimits[start][i].first && point[i] <= gridCellLimits[start][i].second)) {
					inLimits = false;
					break;
				}
			}
			if (inLimits) {
				gridCells[start].addPoint2(point, num_dimensions);
				++insertedPoints;
			}
			else
			{

			}
			return;
		}
		insertPoint2(dim - 1, point, start);
	}

	int visitorIdx = 0;

	unsigned long long resIndex = 0;

	uint32_t probe_sorted_grid_cell_for_matches(uint32_t page, vector<uint32_t> range_start,
		vector<uint32_t> range_end, vector<vector<uint32_t>>& visitor) {

		auto t0 = high_resolution_clock::now();
		uint32_t sum = 0;


		uint32_t cellStartIdx = std::lower_bound(gridCells[page].cellPoints2[0], gridCells[page].cellPoints2[0] + gridCells[page].size, range_start[0]) - gridCells[page].cellPoints2[0];

		if (cellStartIdx > 0)
			cellStartIdx--;
		uint32_t cellEndIdx = std::upper_bound(gridCells[page].cellPoints2[0], gridCells[page].cellPoints2[0] + gridCells[page].size, range_end[0]) - gridCells[page].cellPoints2[0];

		if (cellEndIdx > 0 && cellEndIdx == gridCells[page].size)
			cellEndIdx--;
		auto t1 = high_resolution_clock::now();
		total_refinement_time += t1 - t0;

		//auto t0 = high_resolution_clock::now();
		for (int i = cellStartIdx; i <= cellEndIdx; i++)
		{
			q_num_scanned_points++;
			//q++;
			if (checkMatchingPoint(gridCells[page].cellPoints2, i, visitor, range_start, range_end)) {
				//visitor[resIndex] = gridCells[page].cellPoints2[i];
				//cout << "resIndex = " << resIndex << endl;
				sum++;
			}
			else {
				/*cout << "point mismatch " << "{";
				cout << gridCells[page].cellPoints2[0][i] << "," << gridCells[page].cellPoints2[1][i] << "," << gridCells[page].cellPoints2[2][i];
				cout << "}" << ", cell =  " << page << "\tcellStartIdx = " << i << endl;*/
					
			}

		}
		//auto t1 = high_resolution_clock::now();
		//total_scan_time += t1 - t0;

		//cout << page << ": " << sum << endl;

		return sum;
	}

	bool checkMatchingPoint(uint32_t** point, int idx, vector<vector<uint32_t>>& visitor, const std::vector<uint32_t>& range_start, const std::vector<uint32_t>& range_end)
	{
		for (int i = 0; i < num_dimensions; i++) {
			visitor[i][visitorIdx] = point[i][idx];
			if (!(point[i][idx] >= range_start[i] && point[i][idx] <= range_end[i])) {
				return false;

			}
		}
		//++visitorIdx;
		return true;
	}


public:

	uint32_t range_query_recursive_with_state_main(uint32_t dim, vector<uint32_t> range_start,
		vector<uint32_t> range_end, vector<vector<uint32_t>>& visitor, uint32_t start = 0, uint32_t end = 0, int state = 0)
	{
		uint32_t sum = 0;
		vector<uint32_t> r1(dim);
		vector<uint32_t> r2(dim);
		for (int i = 0; i < num_dimensions; ++i) {
			r1[i] = range_start[orderDims[i]];
			r2[i] = range_end[orderDims[i]];
		}

		for (int i = dim; i > 1; --i) {
			uint32_t tmpIndex = std::lower_bound(attributeSlicers[i - 1].begin(), attributeSlicers[i - 1].end(), r1[i - 1]) - attributeSlicers[i - 1].begin();
			uint32_t tmpIndex2 = std::lower_bound(attributeSlicers[i - 1].begin(), attributeSlicers[i - 1].end(), r2[i - 1]) - attributeSlicers[i - 1].begin();
			if (tmpIndex == attributeSlicers[i - 1].size())
				return 0;
			if (tmpIndex2 == attributeSlicers[i - 1].size())
				tmpIndex2 -= 2; // end() + 1 more element at slicer than slicer_position
			else if (tmpIndex2 != 0)
				tmpIndex2--;
			if (tmpIndex != 0)
				tmpIndex--;

			rangeStartPos[i - 1] = attributeSlicersPointers[i - 1][tmpIndex];
			rangeEndPos[i - 1] = attributeSlicersPointers[i - 1][tmpIndex2];
		}

		return range_query_recursive_with_state_cut(dim, r1, r2, visitor, start, end, state);

	}
	

	uint32_t range_query_recursive_with_state_cut(uint32_t dim, vector<uint32_t> range_start,
		vector<uint32_t> range_end, vector<vector<uint32_t>>& visitor, uint32_t start = 0, uint32_t end = 0, int state = 0)
	{

		uint32_t sum = 0;
		start += rangeStartPos[dim - 1];
		end += rangeEndPos[dim - 1];
		if (start > 0)
			start--;
		if (end > 0)
			end--;

		//stopping condition
		if (dim == 2) {
			//auto t0 = high_resolution_clock::now();
			uint32_t count = 0;
			while (start <= end) {
				//q_num_scanned_cells++;
				//cout << start << endl;
				count += probe_sorted_grid_cell_for_matches(start, range_start, range_end, visitor);
				start += mulfactorArr[dim - 1];
			}
			/*auto t1 = high_resolution_clock::now();
			total_time += t1 - t0;*/
			return count;
		}

		if (dim > 2)
		{
			while (start <= end) {
				sum += range_query_recursive_with_state_cut(dim - 1, range_start, range_end, visitor, start, start, 1); //normal
				start += mulfactorArr[dim - 1];
			}
		}

		return sum;
	}

};
}