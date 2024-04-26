#pragma once

#include <iostream>
#include <vector>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <math.h>

using namespace std;

namespace gf {
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


		GridCell(){
			size = 16;
			cell_number = idx = 0;
		}

		~GridCell()
		{
			//cout << "~GridCell Destructor" << endl;
			for (size_t i = 0; i < dimensions; ++i) {
				delete[] cellPoints2[i];
			}
			delete[] cellPoints2;
		}

		void initializeCells2(int dimensions, int initial_size) {
			//cell_number = ++count;
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

public:
	GridFile() {
		cout << "Constructing GridFile()" << endl;
	}

public:
	GridFile(int dimensions, int num_records, int CELL_SIZE, vector<uint32_t>& datasetMD)
	{
		this->num_records = num_records;
		this->num_dimensions = dimensions;
		this->insertedPoints = 0;
		cout << "*************************************" << endl;
		cout << "Constructing GridFile..." << endl;
		cout << "*************************************" << endl;
		//gridCells = new GridCell[total_cells];
		//for (int i = 0; i < total_cells; i++)
		//	gridCells->initializeCells2(dimensions, initial_size);
		dimsSize.resize(dimensions);
		getDimSizes2(dimensions, CELL_SIZE, num_records);
		configureGrid(num_records);
		configureSlicers(datasetMD);
		configureCellsLimits();
		assignGridNumber();
		cout << "Bulk loading the grid file with " << num_records << " records...." << endl;
		fillGrid(datasetMD);

		for (int i = 1; i <= dimensions; i++) {
			uint32_t mulfactor = 1;
			for (int j = 0; j < i - 1; j++) {
				mulfactor *= dimsSize[j];
			}
			mulfactorArr.push_back(mulfactor);
		}

		for (int i = 0; i < gridCellsCount; i++)
		{
			gridCells[i].shrink_cell(num_dimensions);
		}

		uint32_t minPointsAtAllCells = 1000000;
		uint32_t maxPointsAtAllCells = 0;
		for (int i = 0; i < gridCellsCount; i++) {
			if (gridCells[i].size < minPointsAtAllCells)
				minPointsAtAllCells = gridCells[i].size;
			if (gridCells[i].size > maxPointsAtAllCells)
				maxPointsAtAllCells = gridCells[i].size;
		}
		//cout << "Cell with lowest number of points has=" << minPointsAtAllCells << endl;
		//cout << "Cell with highest number of points has=" << maxPointsAtAllCells << endl;

		rangeStartPos.resize(num_dimensions);
		rangeEndPos.resize(num_dimensions);


	}

	~GridFile()
	{
		//cout << "~GridFile Destructor" << endl;
		delete [] gridCells;
	}
private:
	std::vector<uint32_t> getDimSizes2(uint32_t dimensions, uint32_t PAGE_SIZE, uint32_t dataset_size) {
		uint32_t maxElementPerPage = PAGE_SIZE / (dimensions * 4);
		if (maxElementPerPage == 0)
			maxElementPerPage++;

		//total theoritical needed number of cells. data distribution will define the actual created number of cells.
		uint32_t minRequiredNumberOfGridCells = dataset_size / maxElementPerPage; 

		if (dataset_size % maxElementPerPage != 0)
			minRequiredNumberOfGridCells++;

		double partitionsPerDimension = std::pow(minRequiredNumberOfGridCells, 1.0 / dimensions);
		partitionsPerDimension = ceil(partitionsPerDimension);

		for (int i = 0; i < dimensions; i++)
			dimsSize[i] = partitionsPerDimension;

		cout << "Number of splits per dimension = " << dimsSize[0] << endl;
		//cout << "Minimum required number of grid cells = " << minRequiredNumberOfGridCells << endl;
		min_num_grid_cells = minRequiredNumberOfGridCells;
		min_points_per_cell = maxElementPerPage;

		return dimsSize;
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
		cout << "Total number of grid cells = " << gridCellsCount << endl;
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

	void configureSlicers(uint32_t** dataset) {
		attributeSlicers.resize(dimsSize.size());
		attributeSlicersPointers.resize(dimsSize.size());

		for (int j = 0; j < dimsSize.size(); j++) {
			vector<uint32_t> tmp(num_records);
			for (int i = 0; i < num_records; i++) {
				tmp[i] = dataset[j][i];
			}
			uint32_t minPerDimension = *(std::min_element(tmp.begin(), tmp.end()));
			uint32_t maxPerDimension = *(std::max_element(tmp.begin(), tmp.end()));
			uint32_t splittingValue = (maxPerDimension - minPerDimension) / dimsSize[j];

			//sort(tmp.begin(), tmp.end());
			uint32_t base = tmp.size() / dimsSize[j];

			attributeSlicers[j].push_back(tmp[0]);

			for (uint32_t i = base, l = 1; l < dimsSize[j]; i += base, l++) {
				attributeSlicers[j].push_back(attributeSlicers[j][attributeSlicers[j].size() - 1] + splittingValue);
				attributeSlicersPointers[j].push_back(0);
			}
			attributeSlicers[j].push_back(tmp[tmp.size() - 1]);
			attributeSlicersPointers[j].push_back(0);
		}
	}

	void configureSlicers(vector<uint32_t>& dataset) {
		attributeSlicers.resize(dimsSize.size());
		attributeSlicersPointers.resize(dimsSize.size());

		for (int j = 0; j < dimsSize.size(); j++) {
			vector<uint32_t> tmp(num_records);
			for (int i = 0; i < num_records; i++) {
				tmp[i] = dataset[j*num_records+i];
			}
			
			uint32_t minPerDimension = *(std::min_element(tmp.begin(), tmp.end()));
			uint32_t maxPerDimension = *(std::max_element(tmp.begin(), tmp.end()));
			uint32_t splittingValue = (maxPerDimension - minPerDimension) / dimsSize[j];

			//sort(tmp.begin(), tmp.end());
			uint32_t base = tmp.size() / dimsSize[j];

			attributeSlicers[j].push_back(minPerDimension);

			for (uint32_t i = base, l = 1; l < dimsSize[j]; i += base, l++) {
				attributeSlicers[j].push_back(attributeSlicers[j][attributeSlicers[j].size() - 1] + splittingValue);
				attributeSlicersPointers[j].push_back(0);
			}
			attributeSlicers[j].push_back(maxPerDimension);
			attributeSlicersPointers[j].push_back(0);
		}
	}

	void configureCellsLimits() {
		uint32_t mulfactor = 1;
		for (int j = 0; j < dimsSize.size(); j++) {
			if (j > 0)
				mulfactor *= dimsSize[j - 1];
			uint32_t k = 0;
			for (int i = 0; i < gridCellsCount; i++) {

				if (i % mulfactor == 0)
					k++;

				if (k == attributeSlicers[j].size())
					k = 1;

				if (k == 1)
					gridCellLimits[i].push_back(make_pair(attributeSlicers[j][k - 1], attributeSlicers[j][k]));
				else
					gridCellLimits[i].push_back(make_pair(attributeSlicers[j][k - 1] + 1, attributeSlicers[j][k]));

			}
		}
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
			return;
		}
		insertPoint2(dim - 1, point, start);
	}

	void fillGrid(vector<uint32_t>& dataset) {
		int total_recs = 0;

		uint32_t* p = new uint32_t[num_dimensions]; //placeholder

		for (int i = 0; i < num_records; i++)
		{
			total_recs++;

			//cout << "Inserting ";
			for (int j = 0; j < num_dimensions; j++)
			{
				p[j] = dataset[j * num_records + i];
			}
			insertPoint2(dimsSize.size(), p, 0);
		}
		delete []p; //freeing placeholder
		cout << "Bulk loaded points = " << insertedPoints << endl;
	}

	int visitorIdx = 0;

	bool checkMatchingPoint(uint32_t** point, int idx, vector<vector<uint32_t>>& visitor, const std::vector<uint32_t>& range_start, const std::vector<uint32_t>& range_end)
	{
		for (int i = 0; i < num_dimensions; i++) {
			visitor[i][visitorIdx] = point[i][idx];
			if (!(point[i][idx] >= range_start[i] && point[i][idx] <= range_end[i])) {
				return false;

			}
		}
		++visitorIdx;
		return true;
	}

uint32_t range_query_recurse(uint32_t dim, vector<uint32_t>& range_start,
		vector<uint32_t>& range_end, vector<vector<uint32_t>>& visitor, uint32_t start = 0, uint32_t end = 0)
	{
		uint32_t sum = 0;

		start += rangeStartPos[dim - 1];
		end += rangeEndPos[dim - 1];
		if (start > 0)
			start--;
		if (end > 0)
			end--;

		if (dim == 1) {
			uint32_t count = 0;
			while (start <= end) {
				for (int i = 0; i < gridCells[start].size; i++) {
					//for (int i = 0; i < gridCells[start].cellPoints.size(); i++) {
					q++;
					if (checkMatchingPoint(gridCells[start].cellPoints2,i,visitor, range_start, range_end)) {
						count++;
						//visitor[resIndex] = gridCells[start].cellPoints2[i];
					}
				}
				start++;
			}
			return count;
		}

		while (start <= end && dim > 1) {

			sum += range_query_recurse(dim - 1, range_start, range_end,visitor, start, start); //optimized

			start += mulfactorArr[dim - 1];
		}


		return sum;
	}

public:
	uint32_t range_query(uint32_t dim, vector<uint32_t> range_start,
		vector<uint32_t> range_end, vector<vector<uint32_t>>& visitor, uint32_t start = 0, uint32_t end = 0)
	{
		visitorIdx = 0;
		auto r1 = range_start;
		auto r2 = range_end;

		for (int i = dim; i > 0; --i) {
			uint32_t tmpIndex = std::lower_bound(attributeSlicers[i - 1].begin(), attributeSlicers[i - 1].end(), range_start[i - 1]) - attributeSlicers[i - 1].begin();
			uint32_t tmpIndex2 = std::lower_bound(attributeSlicers[i - 1].begin(), attributeSlicers[i - 1].end(), range_end[i - 1]) - attributeSlicers[i - 1].begin();
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

		return range_query_recurse(dim, range_start, range_end, visitor, start, end);

	}

	
};
}