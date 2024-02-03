#include "DataMatrix.h"

DataMatrix::DataMatrix(const DataMatrix & Q) {
	rowNames = Q.rowNames;
	colNames = Q.colNames;
	rowNameToIdx = Q.rowNameToIdx;
	colNameToIdx = Q.colNameToIdx;
	val.resize(numRows());
	for (int i = 0; i < numRows(); i++) {
		val[i].resize(numCols(), 0);
	}
}

DataMatrix::DataMatrix(const std::vector< std::string > & _rowNames, const std::vector< std::string > & _colNames, const std::map< std::string, int > & _rowNameToIdx, const std::map< std::string, int > & _colNameToIdx) {
	rowNames = _rowNames;
	colNames = _colNames;
	rowNameToIdx = _rowNameToIdx;
	colNameToIdx = _colNameToIdx;
	val.resize(numRows());
	for (int i = 0; i < numRows(); i++) {
		val[i].resize(numCols(), 0);
	}
}

int DataMatrix::getRowIdx(const std::string & rowName) const {
	auto it = rowNameToIdx.find(rowName);
	if (it != rowNameToIdx.end()) return it->second;
	else {
		std::cerr << "\n\n<ERROR in DataMatrix::rowIdx> No row found with the name '" << rowName << "'. Exiting.\n";
		exit(0);
	}
}

int DataMatrix::getColIdx(const std::string & colName) const {
	auto it = colNameToIdx.find(colName);
	if (it != colNameToIdx.end()) return it->second;
	else {
		std::cerr << "\n\n<ERROR in DataMatrix::rowIdx> No row found with the name '" << colName << "'. Exiting.\n";
		exit(0);
	}
}

long double DataMatrix::getVal(const std::string & rowName, const std::string & colName) const {
	const int rowIdx = getRowIdx(rowName);
	const int colIdx = getColIdx(colName);
	return val[rowIdx][colIdx];
}

void DataMatrix::setVal(const std::string & rowName, const std::string & colName, const long double x) {
	const int rowIdx = getRowIdx(rowName);
	const int colIdx = getColIdx(colName);
	val[rowIdx][colIdx] = x;
}

// Intersects the rows of the current DataMatrix (based on their names) with the rows of the DataMatrix Q, discarding any that are not present in Q. This function puts the kept rows in the lexicographic order (by their names). The columns of the two DataMatrices need not match.
void DataMatrix::intersectRows(const DataMatrix & Q) {
	std::vector< std::string > new_rowNames;
	std::map< std::string, int > new_rowNameToIdx;
	std::vector< std::vector<long double> > new_val;
	for (auto & it : rowNameToIdx) {
		if (Q.rowNameToIdx.count(it.first)) {
			int idx = new_rowNameToIdx.size();
			new_rowNameToIdx[it.first] = idx;
			new_rowNames.push_back(it.first);
		}
	}
	new_val.resize(new_rowNames.size());
	for (int i = 0; i < new_rowNames.size(); i++) {
		const std::string & rowName = new_rowNames[i];
		int oldIdx = rowNameToIdx[rowName];
		new_val[i] = val[oldIdx];
	}
	rowNames = new_rowNames;
	rowNameToIdx = new_rowNameToIdx;
	val = new_val;
}

// The same as the above, but it intersects the row names with those from an std::map, rather than the rows of another DatMatrix.
void DataMatrix::intersectRows(const std::map< std::string, int > & M) {
	std::vector< std::string > new_rowNames;
	std::map< std::string, int > new_rowNameToIdx;
	std::vector< std::vector<long double> > new_val;
	for (auto & it : rowNameToIdx) {
		if (M.count(it.first)) {
			int idx = new_rowNameToIdx.size();
			new_rowNameToIdx[it.first] = idx;
			new_rowNames.push_back(it.first);
		}
	}
	new_val.resize(new_rowNames.size());
	for (int i = 0; i < new_rowNames.size(); i++) {
		const std::string & rowName = new_rowNames[i];
		int oldIdx = rowNameToIdx[rowName];
		new_val[i] = val[oldIdx];
	}
	rowNames = new_rowNames;
	rowNameToIdx = new_rowNameToIdx;
	val = new_val;
}

// Intersects the columns of the current DataMatrix (based on their names) with the columns of the DataMatrix Q, discarding any that are not present in Q. This function puts the kept columns in a lexicographic order (by their names). The rows of the two DataMatrices need not match.
void DataMatrix::intersectColumns(const DataMatrix & Q) {
	std::vector< std::string > new_colNames;
	std::map< std::string, int > new_colNameToIdx;
	std::vector< std::vector<long double> > new_val;

	for (auto & it : colNameToIdx) {
		if (Q.colNameToIdx.count(it.first)) {
			int idx = new_colNameToIdx.size();
			new_colNameToIdx[it.first] = idx;
			new_colNames.push_back(it.first);
		}
	}
	new_val.resize(numRows());
	std::vector<long double> newRow(new_colNameToIdx.size());
	for (int i = 0; i < numRows(); i++) {
		for (int j = 0; j < new_colNameToIdx.size(); j++) {
			const std::string & colName = new_colNames[j];
			const int colIdx = colNameToIdx[colName];
			newRow[j] = val[i][colIdx];
		}
		val[i] = newRow;
	}

	colNames = new_colNames;
	colNameToIdx = new_colNameToIdx;
	val = new_val;
}
// The same as the above, but it intersects the column names with those from an std::map, rather than the columns of another DatMatrix.
void DataMatrix::intersectColumns(const std::map< std::string, int > & M) {
	std::vector< std::string > new_colNames;
	std::map< std::string, int > new_colNameToIdx;
	std::vector< std::vector<long double> > new_val;

	for (auto & it : colNameToIdx) {
		if (M.count(it.first)) {
			int idx = new_colNameToIdx.size();
			new_colNameToIdx[it.first] = idx;
			new_colNames.push_back(it.first);
		}
	}
	new_val.resize(numRows());
	std::vector<long double> newRow(new_colNameToIdx.size());
	for (int i = 0; i < numRows(); i++) {
		for (int j = 0; j < new_colNameToIdx.size(); j++) {
			const std::string & colName = new_colNames[j];
			const int colIdx = colNameToIdx[colName];
			newRow[j] = val[i][colIdx];
		}
		val[i] = newRow;
	}

	colNames = new_colNames;
	colNameToIdx = new_colNameToIdx;
	val = new_val;
}

void DataMatrix::generateSyntheticDataForTheGammaMatrix(const DataMatrix & B, const DataMatrix & S, const DataMatrix & F, const std::vector< CellTypeMarkerEntry > & markers, const std::string & outputFolder) {
	{
		std::string dim = "";
		if (this->numRows() != B.numRows()) dim = "rows";
		else if (this->numCols() != B.numCols()) dim = "columns";
		if (dim.length() > 0) {
			std::cerr << "\n\n<ERROR in DataMatrix::generateRandomValuesForTheGammaMatrix> The number of " << dim << " between Gamma and B do not match. Exiting.\n";
			exit(0);
		}
	}
	std::random_device rd;
	std::mt19937 RNG(rd());

	long double smallestNonzeroSig = -1;
	for (int j = 0; j < S.numRows(); j++) {
		if (S.val[j][0] > 0) {
			if ((smallestNonzeroSig == -1) or (smallestNonzeroSig > S.val[j][0])) {
				smallestNonzeroSig = S.val[j][0];
			}
		}
	}

	std::vector<long double> statsMin(F.numRows());
	std::vector<long double> statsMax(F.numRows());
	std::vector<long double> statsMedian(F.numRows());
	std::vector<long double> statsAvg(F.numRows());
	for (int k = 0; k < F.numRows(); k++) {
		std::vector<long double> rowCopy(F.val[k]);
		std::sort(rowCopy.begin(), rowCopy.end());
		statsMin[k] = rowCopy[0];
		statsMax[k] = rowCopy[F.numCols() - 1];
		statsMedian[k] = rowCopy[F.numCols()/2];
		statsAvg[k] = 0;
		for (const long double & val : rowCopy) {
			statsAvg[k] += val;
		}
		statsAvg[k] /= rowCopy.size();
	}

	std::vector<int> geneIndices(numRows());
	std::ofstream fout(outputFolder + "/statsGamma.tsv");
	if (!fout) {
		std::cerr << "\n<ERROR in DataMatrix::generateSyntheticDataForTheGammaMatrix> Cannot open [" << (outputFolder + "/statsGamma.tsv") << "] for writing.\n";
		exit(0);
	}
	// Mode where a gene can have at most a single non-zero gamma assigned to it
	std::vector<bool> alreadyDone(numRows(), false);
	std::vector<int> colOrder(numCols());
	for (int k = 0; k < numCols(); k++) colOrder[k] = k;
	std::shuffle(colOrder.begin(), colOrder.end(), RNG);
	for (const int k : colOrder) {
		fout << colNames[k] << '\t';
		// std::cout << F.rowNames[k];
		// if (statsMax[k] < fractionEPS) { std::cerr << "Cell type'" << F.rowNames[k] << "' has all fractions below " << fractionEPS << ". Not assigning gamma to it.\n"; continue; }
		if (statsMedian[k] < fractionEPS) { std::cerr << "Cell type'" << F.rowNames[k] << "' has the median fractions below " << fractionEPS << ". Not assigning gamma to it.\n"; }
		// if (statsAvg[k] < fractionEPS) { std::cerr << "Cell type'" << F.rowNames[k] << "' has the average fractions below " << fractionEPS << ". Not assigning gamma to it.\n"; continue; }
		else {
			int numGenesExpressed = 0;
			for (int i = 0; i < numRows(); i++) {
				if (B.val[i][k] > 1 && !alreadyDone[i]) {
					geneIndices[numGenesExpressed] = i;
					numGenesExpressed++;
				}
			}
			if (numGenesExpressed == 0) continue;

			// Sort by reference count, in descending order
			std::sort(geneIndices.begin(), geneIndices.begin() + numGenesExpressed, [k, &B] (const int i1, const int i2) { return B.val[i1][k] > B.val[i2][k]; } );

			// For each reference count, find the minimum possible gamma value that would still result in the nonnegative reference count for every B + gamma * S_j. Select only those that allow value to be lower than or equal to 5
			std::vector<long double> minPossibleGamma(numGenesExpressed);
			int numToAssignNegativeGamma = 0;
			for (int i = 0; i < numGenesExpressed; i++) {
				const int geneIdx = geneIndices[i];
				minPossibleGamma[i] = -B.val[geneIdx][k];
				for (int j = 0; j < S.numRows(); j++) {
					if (S.val[j][0] > 0) {
						const long double temp = -B.val[geneIdx][k] / S.val[j][0];
						if (minPossibleGamma[i] < temp) minPossibleGamma[i] = temp;
					}
				}
				if (minPossibleGamma[i] <= -60) minPossibleGamma[i] = -60;
				if (minPossibleGamma[i] <= -5.0) numToAssignNegativeGamma++;
				else break;
			}
			int numToAssignTotal = numGenesExpressed * 0.1;
			if (numToAssignTotal == 0) continue;
			if (numToAssignNegativeGamma > numToAssignTotal/2) numToAssignNegativeGamma = numToAssignTotal/2;
			const int numToAssignPositiveGamma = numToAssignNegativeGamma;
			numToAssignTotal = numToAssignPositiveGamma + numToAssignNegativeGamma;

			fout << "(all) " << numGenesExpressed << "\t(+) " << numToAssignPositiveGamma << '\n';
			fout << numToAssignNegativeGamma << "(-)";

			for (int i = 0; i < numToAssignNegativeGamma; i++) {
				// fprintf(stderr, "%Lf\t", minPossibleGamma[i]);
				fout << '\t' << minPossibleGamma[i];
				const int geneIdx = geneIndices[i];
				std::uniform_real_distribution<> uniReal(minPossibleGamma[i], -5.0);
				val[geneIdx][k] = uniReal(RNG);
				alreadyDone[geneIdx] = true;
			}
			// fprintf(stderr, "\n");
			for (int i = numToAssignNegativeGamma; i < numToAssignTotal; i++) {
				const int geneIdx = geneIndices[i];
				std::uniform_real_distribution<> uniReal(5.0, 60.0);
				val[geneIdx][k] = uniReal(RNG);
				alreadyDone[geneIdx] = true;
			}
		}
		fout << '\n';
	}
	fout.close();
}

// Generates random values to simulate a Gamma matrix. (gene, cell type) pairs that are not expressed in the single cell reference matrix B are given a value of 0.
void DataMatrix::generateRandomValuesForTheGammaMatrix(const DataMatrix & B) {
	{
		std::string dim = "";
		if (this->numRows() != B.numRows()) dim = "rows";
		else if (this->numCols() != B.numCols()) dim = "columns";
		if (dim.length() > 0) {
			std::cerr << "\n\n<ERROR in DataMatrix::generateRandomValuesForTheGammaMatrix> The number of " << dim << " between Gamma and B do not match. Exiting.\n";
			exit(0);
		}
	}
	const int RND_MAXABSVAL = 200;
	// const int RND_MINABSVAL = 10;
	const int RND_MAXDECIMALPRECISION = 2;
	const int RND_TOTALVAL = RND_MAXABSVAL * pow(10, RND_MAXDECIMALPRECISION);
	std::random_device rd;
	std::mt19937 * RNG = new std::mt19937(rd());
	std::uniform_int_distribution<int> uni(0, RND_TOTALVAL);
	// std::uniform_int_distribution<int> uniZero(0, 100);
	for (int i = 0; i < numRows(); i++) {
		for (int j = 0; j < numCols(); j++) {
			if (B.val[i][j] == 1) val[i][j] = 0;
			else {
				// if (uniZero(*RNG) <= 80) val[i][j] = 0;
				// else {
					val[i][j] = ((long double)(uni(*RNG) - RND_TOTALVAL/2)) / pow(10, RND_MAXDECIMALPRECISION);
					// if (val[i][j] < 0) val[i][j] -= RND_MINABSVAL;
					// else if (val[i][j] > 0) val[i][j] += RND_MINABSVAL;
				// }
			}
		}
	}
	delete RNG;
}

// Generates random values to simulate a cell-type fraction matrix F. 
void DataMatrix::generateRandomValuesForTheFractionMatrix(const long double correlationNoseFactor) {
	// const int RND_MAXVAL = 100000;
	std::random_device rd;
	std::mt19937 RNG(rd());
	// std::mt19937 RNG(29); // FIXED SEED
	// std::uniform_int_distribution<int> uni(0, RND_MAXVAL);
	std::exponential_distribution<> dist(1);
	// std::uniform_int_distribution<int> uniZero(0, 100);
	for (int k = 0; k < numRows(); k++) {
		// for (int j = 0; j < numCols(); j++) {
		// 	if (uniZero(RNG) <= 25) val[k][j] = 0;
		// 	else val[k][j] = 3 * dist(RNG);
		// }
		// if (uniZero(RNG) <= 20) val[k][0] = 0;
		// else 
		val[k][0] = 2 * dist(RNG);
	}
	{	// Turn the first column into a vector of percentages summing up to 1
		const int j = 0;
		long double totalSum = 0;
		for (int k = 0; k < numRows(); k++) {
			totalSum += val[k][j];
		}
		for (int k = 0; k < numRows(); k++) {
			val[k][j] = val[k][j] / totalSum;
		}
	}
	for (int j = 1; j < numCols(); j++) {
		for (int k = 0; k < numRows(); k++) {
			// if (val[k][0] > fractionEPS) {
				std::uniform_real_distribution<> uni(-val[k][0] * correlationNoseFactor, val[k][0] * correlationNoseFactor);
				val[k][j] = val[k][0] + uni(RNG);
			// }
			// else {
			// 	std::uniform_real_distribution<> uni(0, smallestNonZeroEntry * correlationNoseFactor);
			// 	val[k][j] = uni(RNG);
			// }
		}
		long double totalSum = 0;
		for (int k = 0; k < numRows(); k++) {
			totalSum += val[k][j];
		}
		for (int k = 0; k < numRows(); k++) {
			val[k][j] = val[k][j] / totalSum;
		}
	}
	printColumnCorrelationStats("Correlated vectors");
}

// Generates random values to simulate a mutational signature strength matrix S. 
void DataMatrix::generateRandomValuesForTheSignatureMatrix() {
	const int RND_MAXABSVAL = 300;
	// const int RND_MAXDECIMALPRECISION = 4;
	// const int RND_TOTALVAL = RND_MAXABSVAL * pow(10, RND_MAXDECIMALPRECISION);
	std::random_device rd;
	std::mt19937 * RNG = new std::mt19937(rd());
	// std::uniform_int_distribution<int> uni(0, RND_TOTALVAL);
	std::uniform_real_distribution<> uni(1, RND_MAXABSVAL);
	std::uniform_int_distribution<int> uniZero(0, 100);
	for (int i = 0; i < numRows(); i++) {
		for (int j = 0; j < numCols(); j++) {
			if (uniZero(*RNG) <= 20) val[i][j] = 0;
			// else val[i][j] = uni(*RNG) / pow(10, RND_MAXDECIMALPRECISION);
			else val[i][j] = uni(*RNG);
		}
	}
	delete RNG;
}

// Calculates the bulk expression based on the model, and adds noise to within noiseFraction % of the calculated model value. 
void DataMatrix::calculateBulkExpressionWithNoise(const DataMatrix & B, const DataMatrix & S, const std::string & sigName, const DataMatrix & F, const DataMatrix & Gamma, const DataMatrix & Alpha, const long double noiseFraction) {
	// int n = numCols();	// the number of samples
	// int m = numRows();	// the number of genes
	// int c = B.numCols();	// the number of cell types
	const int sigIdx = S.colNameToIdx.find(sigName)->second;
	std::random_device rd;
	std::mt19937 RNG(rd());
	for (int i = 0; i < numRows(); i++) {
		for (int j = 0; j < numCols(); j++) {
			const std::string & sampleName = colNames[j];
			if (S.rowNameToIdx.count(sampleName)) {
				const int idx = S.rowNameToIdx.find(sampleName)->second;
				const long double Sj = S.val[idx][sigIdx];
				val[i][j] = 0;
				for (int k = 0; k < B.numCols(); k++) {
					val[i][j] += F.val[k][j] * (B.val[i][k] + Sj * Gamma.val[i][k] + Alpha.val[i][k]);
					if (noiseFraction > 0) {
						std::uniform_real_distribution<> uni(-val[i][j] * noiseFraction, val[i][j] * noiseFraction);
						val[i][j] += uni(RNG);
					}
				}
			}
			else {
				std::cerr << "<ERROR in DataMatrix::calculateBulkExpression> Could not find sample name " << sampleName << " in the signature matrix. Exiting.\n";
				exit(0);
			}
		}
	}
}

// Reorders the columns based on the order provided in 'new_colNames'
void DataMatrix::reorderColumns(const std::vector< std::string > & new_colNames) {
	if (new_colNames.size() != numCols()) {
		std::cerr << "<ERROR in DataMatrix::reorderColumns> the column sizes do not match. Exiting.\n";
		exit(0);
	}
	std::map< std::string, int > new_colNameToIdx;
	std::vector<int> indices(numCols(), 0);
	for (int j = 0; j < numCols(); j++) {
		const std::string & colName = new_colNames[j];
		const std::string & colName_sanitized = colName.substr(0, colName.find("|"));
		if (!colNameToIdx.count(colName_sanitized)) {
			std::cerr << "<ERROR in DataMatrix::reorderColumns> Cannot find column named " << colName_sanitized << " in the matrix. Exiting.\n";
			exit(0);
		}
		indices[j] = colNameToIdx[colName_sanitized];
		new_colNameToIdx[colName] = j;
	}
	std::vector<long double> newRow(numCols());
	for (int i = 0; i < numRows(); i++) {
		for (int j = 0; j < numCols(); j++) {
			newRow[j] = val[i][ indices[j] ];
		}
		val[i] = newRow;
	}

	colNames = new_colNames;
	colNameToIdx = new_colNameToIdx;
}

// Reorders the columns based on the order provided in 'new_colNames'
void DataMatrix::reorderRows(const std::vector< std::string > & new_rowNames) {
	if (new_rowNames.size() != numRows()) {
		std::cerr << "<ERROR in DataMatrix::reorderRows> the row sizes do not match. Exiting.\n";
		exit(0);
	}
	std::map< std::string, int > new_rowNameToIdx;
	std::vector<int> indices(numRows(), 0);
	for (int i = 0; i < numRows(); i++) {
		const std::string & rowName = new_rowNames[i];
		const std::string & rowName_sanitized = rowName.substr(0, rowName.find("|"));
		if (!rowNameToIdx.count(rowName_sanitized)) {
			std::cerr << "<ERROR in DataMatrix::reorderRows> Cannot find row named " << rowName_sanitized << " in the matrix. Exiting.\n";
			exit(0);
		}
		indices[i] = rowNameToIdx[rowName_sanitized];
		new_rowNameToIdx[rowName] = i;
	}
	std::vector< std::vector<long double> > new_val(numRows());
	for (int i = 0; i < numRows(); i++) {
		new_val[i] = val[ indices[i] ];
	}

	val = new_val;
	rowNames = new_rowNames;
	rowNameToIdx = new_rowNameToIdx;
}

long double DataMatrix::correlationPearson(const int j1, const int j2) const {
	long double r = 0;
	long double avgA = 0, avgB = 0;
	long double stdA = 0, stdB = 0;
	for (int i = 0; i < numRows(); i++) {
		avgA += val[i][j1];
		avgB += val[i][j2];
	}
	avgA /= numRows();
	avgB /= numRows();
	for (int i = 0; i < numRows(); i++) {
		r += (val[i][j1] - avgA) * (val[i][j2] - avgB);
		stdA += (val[i][j1] - avgA) * (val[i][j1] - avgA);
		stdB += (val[i][j2] - avgB) * (val[i][j2] - avgB);
	}
	if (stdA == 0) {
		std::cerr << "<ERROR in DataMatrix::correlationPearson> Centered column j1 is a zero vector. Exiting.\n";
		exit(0);
	}
	if (stdB == 0) {
		std::cerr << "<ERROR in DataMatrix::correlationPearson> Centered column j2 is a zero vector. Exiting.\n";
		exit(0);
	}
	r = r / (sqrt(stdA * stdB));
	return r;
}

long double DataMatrix::correlationPearson(const std::vector<long double> & A, const std::vector<long double> & B) const {
	if (A.size() != B.size()) {
		std::cerr << "<ERROR in DataMatrix::correlationPearson> the sizes of A (" << A.size() << ") and B (" << B.size() << ") do not match. Exiting.\n";
		exit(0);
	}
	const int n = A.size();
	long double r = 0;
	long double avgA = 0, avgB = 0;
	long double stdA = 0, stdB = 0;
	for (int i = 0; i < n; i++) {
		avgA += A[i];
		avgB += B[i];
	}
	avgA /= n;
	avgB /= n;
	for (int i = 0; i < n; i++) {
		r += (A[i] - avgA) * (B[i] - avgB);
		stdA += (A[i] - avgA) * (A[i] - avgA);
		stdB += (B[i] - avgB) * (B[i] - avgB);
	}
	if (stdA == 0) {
		std::cerr << "<ERROR in DataMatrix::correlationPearson> Centered A is a zero vector. Exiting.\n";
		exit(0);
	}
	if (stdB == 0) {
		std::cerr << "<ERROR in DataMatrix::correlationPearson> Centered B is a zero vector. Exiting.\n";
		exit(0);
	}
	r = r / (sqrt(stdA * stdB));
	return r;
}

void DataMatrix::printColumnCorrelationStats(const char * msg) const {
	std::vector<long double> corrMin(numCols()), corrMax(numCols()), corrAvg(numCols());
	for (int j1 = 0; j1 < numCols(); j1++) {
		corrMin[j1] = 5;
		corrMax[j1] = -5;
		corrAvg[j1] = 0;
		for (int j2 = 0; j2 < numCols(); j2++) {
			if (j1 != j2) {
				long double r = correlationPearson(j1, j2);
				if (r < corrMin[j1]) corrMin[j1] = r;
				if (r > corrMax[j1]) corrMax[j1] = r;
				corrAvg[j1] += r;
			}
		}
		corrAvg[j1] /= (numCols() - 1);
	}
	fprintf(stderr, "%s\n", msg);
	fprintf(stderr, "%12s\t%9s\t%9s\t%9s\n", "Sample", "corrMin", "corrMax", "corrAvg");
	for (int j = 0; j < numCols(); j++) {
		fprintf(stderr, "%8d\t%9Lf\t%9Lf\t%9Lf\n", j, corrMin[j], corrMax[j], corrAvg[j]);
	}
}