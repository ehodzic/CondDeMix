#ifndef __DATAMATRIX_H_INCLUDED__
#define __DATAMATRIX_H_INCLUDED__

#include <vector>
#include <map>
#include <set>
#include <random>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

const long double EPS = 1e-6;
const long double fractionEPS = 1e-4;

struct CellTypeMarkerEntry {
	std::string cellTypeName;
	std::string geneName;
	long double distance;
	int numCellTypesGeneIsExpressedIn;

	bool operator<(const CellTypeMarkerEntry & other) const {
		if (cellTypeName != other.cellTypeName) return cellTypeName < other.cellTypeName;
		if (distance != other.distance) return distance > other.distance;
		else return geneName < other.geneName;
	}
};

struct DataMatrix {
	std::vector< std::string > rowNames;
	std::vector< std::string > colNames;
	std::map< std::string, int > rowNameToIdx;
	std::map< std::string, int > colNameToIdx;
	std::vector< std::vector<long double> > val;

	DataMatrix() {}
	DataMatrix(const DataMatrix & Q);
	DataMatrix(const std::vector< std::string > & _rowNames, const std::vector< std::string > & _colNames, const std::map< std::string, int > & _rowNameToIdx, const std::map< std::string, int > & _colNameToIdx);

	int numRows() const { return rowNameToIdx.size(); }
	int numCols() const { return colNameToIdx.size(); }

	int getRowIdx(const std::string & rowName) const;
	int getColIdx(const std::string & colName) const;
	long double getVal(const std::string & rowName, const std::string & colName) const;
	void setVal(const std::string & rowName, const std::string & colName, const long double x);

	void intersectRows(const DataMatrix & Q);
	void intersectRows(const std::map< std::string, int > & M);
	void intersectColumns(const DataMatrix & Q);
	void intersectColumns(const std::map< std::string, int > & M);
	void generateSyntheticDataForTheGammaMatrix(const DataMatrix & B, const DataMatrix & S, const DataMatrix & F, const std::vector< CellTypeMarkerEntry > & markers, const std::string & outputFolder);
	void generateRandomValuesForTheGammaMatrix(const DataMatrix & B);
	void generateRandomValuesForTheFractionMatrix(const long double correlationNoseFactor);
	void generateRandomValuesForTheSignatureMatrix();
	void calculateBulkExpressionWithNoise(const DataMatrix & B, const DataMatrix & S, const std::string & sigName, const DataMatrix & F, const DataMatrix & Gamma, const DataMatrix & Alpha, const long double noiseFraction);
	void reorderColumns(const std::vector< std::string > & newOrder);
	void reorderRows(const std::vector< std::string > & new_rowNames);

	void printColumnCorrelationStats(const char * msg) const;
	long double correlationPearson(const int j1, const int j2) const;
	long double correlationPearson(const std::vector<long double> & A, const std::vector<long double> & B) const;
};

#endif