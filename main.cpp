#include <iostream>
#include <cassert>
#include "IO.h"
#include "DataMatrix.h"
#include "method.h"

int main(int argc, char * argv[]) {
	bool genSyntheticData = false;
	bool redrawSamples = false;
	bool initializeRandom = false;
	DataMatrix M, B, S;
	std::string sigName;
	std::vector< CellTypeMarkerEntry > markers;

	IO::readInputParameters(argc, argv);
	genSyntheticData = (IO::consoleParameters['t'] == "1");
	redrawSamples = (IO::consoleParameters['l'] != "0");
	initializeRandom = (IO::consoleParameters['r'] == "1");

	IO::readMatrix(IO::consoleParameters['b'], B);

	if (genSyntheticData) {	// We are just generating synthetic data, nothing else
		const int n = 100;	// the number of samples
		sigName = "sigX";
		// B_{m x c}, F_{c x n}, M_{m x n}
		std::vector< std::string > sampleNames(n), sigNames(1);
		std::map< std::string, int > sampleIndices, sigIndices;
		char temp[100];
		for (int i = 0; i < n; i++) {
			sprintf(temp, "sample_%d", i);
			sampleNames[i] = temp;
			sampleIndices[temp] = i;
		}
		sigNames[0] = sigName;
		sigIndices[sigName] = 0;

		// Structure of F_{c x n}, M_{m x n}, S_{n x 1}
		// generated F
		DataMatrix F(B.colNames, sampleNames, B.colNameToIdx, sampleIndices);
		M = DataMatrix(B.rowNames, sampleNames, B.rowNameToIdx, sampleIndices);
		S = DataMatrix(sampleNames, sigNames, sampleIndices, sigIndices);
		DataMatrix gamma(B);
		DataMatrix alpha(B);
		S.generateRandomValuesForTheSignatureMatrix();

		long double bulkNoiseFactor;
		long double fractionDeviationFactor;
		sscanf(IO::consoleParameters['x'].c_str(), "%Lf", &bulkNoiseFactor);
		sscanf(IO::consoleParameters['y'].c_str(), "%Lf", &fractionDeviationFactor);
		F.generateRandomValuesForTheFractionMatrix(fractionDeviationFactor);
		
		gamma.generateSyntheticDataForTheGammaMatrix(B, S, F, markers, IO::outputFolder);

		M.calculateBulkExpressionWithNoise(B, S, sigName, F, gamma, alpha, bulkNoiseFactor);
		// Just output the generated data
		IO::printSyntheticInput(B, F, M, gamma, alpha, S);
	}
	else {
		IO::readMatrix(IO::consoleParameters['b'], B);
		IO::readMatrix(IO::consoleParameters['m'], M);
		IO::readMatrix(IO::consoleParameters['s'], S);
		if (IO::consoleParameters.count('k')) IO::readCellTypeMarkerEntries(IO::consoleParameters['k'], markers);
		std::cout << "Intersecting the genes of M and B, and samples of M and S...\n";
		M.intersectRows(B);
		B.intersectRows(M);
		M.intersectColumns(S.rowNameToIdx);
		S.intersectRows(M.colNameToIdx);
		std::cout << "M (" << M.numRows() << " x " << M.numCols() <<")\n";
		std::cout << "B (" << B.numRows() << " x " << B.numCols() <<")\n";
		std::cout << "S (" << S.numRows() << " x " << S.numCols() <<")\n";
		assert(M.rowNames == B.rowNames);
		assert(M.colNames == S.rowNames);
		DataMatrix F(B.colNames, M.colNames, B.colNameToIdx, M.colNameToIdx), gamma(B), alpha(B);
		std::cout << "F (" << F.numRows() << " x " << F.numCols() <<")\n";
		std::cout << "Gamma (" << gamma.numRows() << " x " << gamma.numCols() <<")\n";
		std::cout << "Alpha (" << alpha.numRows() << " x " << alpha.numCols() <<")\n";
		if (S.numCols() > 1) {
			if (!IO::consoleParameters.count('c')) {
				std::cerr << "\n<ERROR> The signature matrix has multiple columns, but the parameter -c is missing. Program does not know which signature column to use in the model. Exiting.\n";
				exit(0);
			}
			else if (IO::consoleParameters.count('c') && !S.colNameToIdx.count(IO::consoleParameters['c'])) {
				std::cerr << "\n<ERROR> The signature matrix has multiple columns, but no column is named as provided in the -c parameter. Exiting.\n";
				exit(0);
			}
			sigName = IO::consoleParameters['c'];
		}
		else {
			sigName = S.colNames[0];
		}
		std::cout << "Signature " << sigName << " used in the model.\n";
		
		if (!S.colNameToIdx.count(sigName)) {
			std::cerr << "\n\n<ERROR> Signature matrix has no column named '" << sigName << "'. Exiting.\n";
			exit(0);
		}
		if (initializeRandom) { gamma.generateRandomValuesForTheGammaMatrix(B); }
		runMethod(M, B, S, sigName, markers, F, gamma, alpha);
	}
	return 0;
}