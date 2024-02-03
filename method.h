#ifndef __METHOD_H_INCLUDED__
#define __METHOD_H_INCLUDED__

#include "IO.h"
#include "DataMatrix.h"
#include "gurobi_c++.h"
#include <chrono>
#include <cmath>

extern std::map< int, std::string > GUROBI_OPT_STATUS_CODES;
const long double MIN_VAL_CUTOFF = 1e-6;
const long double MIN_NONZERO_FRACTION = MIN_VAL_CUTOFF;
const long double OBJ_COEFF_ALPHA = 0.0000001;
const int MAX_NUM_JUMPS = 10;

void solveLP(const DataMatrix & M, const DataMatrix & B, const DataMatrix & S, const std::string & sigName, const std::set< std::string > & markers, DataMatrix & F, DataMatrix & gamma, DataMatrix & alpha, long double & objVal, const int maxNonZeroGammas, const bool extrapolateFromMarkers = false);

void runMethod(const DataMatrix & M, const DataMatrix & B, const DataMatrix & S, const std::string & sigName, const std::vector< CellTypeMarkerEntry > & markers, DataMatrix & F, DataMatrix & gamma, DataMatrix & alpha);

long double solveStageTwo(const std::vector< std::string > & samples, const std::vector< std::string > cellTypes, const DataMatrix & B, const DataMatrix & S, const DataMatrix & M, const DataMatrix & solF, const int sigIdx, const std::string & geneName, const std::vector<int> & expressedInCellTypes, GRBEnv * env, DataMatrix & solAlpha, DataMatrix & solGamma, const int nonZeroIdx);

long double getSquareError(const std::vector< std::string > & samples, const std::vector< std::string > cellTypes, const DataMatrix & B, const DataMatrix & S, const DataMatrix & M, const DataMatrix & solF, const int sigIdx, const std::string & geneName, const std::vector<int> & expressedInCellTypes, const DataMatrix & solAlpha, const DataMatrix & solGamma);

#endif