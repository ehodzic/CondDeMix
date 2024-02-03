#ifndef __IO_H_INCLUDED__
#define __IO_H_INCLUDED__

#include <map>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cctype>
#include "DataMatrix.h"

namespace IO {
	extern std::map<char, std::string> consoleParameters;
	extern std::string outputFolder;

	void printUsage();
	void printHeader(const std::string & text, const int width = 0, const char fillh = '*', const char fillv = '*', const char fillc = '*');
	void printRow(const std::string & text, const int width = 0, const char fillv = '*', const int padding = 0);
	void printBreak(const int width = 0, const char fillh = '*', const char fillc = '*');
	void printMatrix(const std::string & filename, const DataMatrix & M, const std::string & txt = "NAME", const DataMatrix * const B = nullptr);
	void printSyntheticInput(const DataMatrix & B, const DataMatrix & F, const DataMatrix & M, const DataMatrix & gamma, const DataMatrix & alpha, const DataMatrix & S);
	void printSolution(const int iter, const DataMatrix & B, const DataMatrix & F, const DataMatrix & M, const DataMatrix & gamma, const DataMatrix & alpha, const DataMatrix & S, const std::string & sigName, const bool onlyBasic = false);
	void printDone();
	void printStats(const int iter, const long double lastL2Error, const long double objVal);
	void printHeatmapInput(const int iter, const DataMatrix & F, const DataMatrix & S);

	void readInputParameters(const int argc, char * argv[]);
	void readMatrix(const std::string & filename, DataMatrix & M);
	void readStringSet(const std::string & filename, std::set<std::string> & S);
	void readCellTypeMarkerEntries(const std::string & filename, std::vector< CellTypeMarkerEntry > & V);
}

#endif