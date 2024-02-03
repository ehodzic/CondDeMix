#include "IO.h"

namespace IO {
	std::map<char, std::string> consoleParameters;
	std::string outputFolder;

	void printHeader(const std::string & text, const int width /*= 0*/, const char fillh /*= '*'*/, const char fillv /*= '*'*/, const char fillc /*= '*'*/) {
		using std::cout;
		using std::endl;
		using std::setfill;
		using std::setw;
		int n = text.length();
		if (width > n) n = width;
		int before = n/2 - text.length()/2;
		int after = (n - n/2) - (text.length() - text.length()/2);
		cout << std::right << endl;
		cout << fillc << setfill(fillh) << setw(n + 3) << fillc << "\n";
		cout << fillv << setfill(' ') << setw(before + text.length() + 1) << text << setw(after + 2) << fillv << endl;
		cout << fillc << setfill(fillh) << setw(n + 3) << fillc << "\n";
	}

	void printRow(const std::string & text, const int width /*= 0*/, const char fillv /*= '*'*/, const int padding /*= 0*/) {
		if (padding >= width) {
			std::cerr <<  "<ERROR in IO::printRow> Padding is wider than the row width. Exiting.\n";
			exit(0);
		}
		for (int i = 0; i < text.length(); ) {
			while (i < text.length() && text[i] == ' ') i++;
			int j = std::min(i + width - padding, int(text.length()));
			if (j < text.length() && !isspace(text[j])) {
				while (j > i && text[j - 1] != ' ') j--;
			}
			std::cout << fillv << " " << std::left << std::setfill(' ') << std::setw(padding) << "" << std::setw(width - padding) << text.substr(i, j - i) << " " << fillv << std::endl;
			i = j;
		}
	}

	void printBreak(const int width /*= 0*/, const char fillh /*= '*'*/, const char fillc /*= '*'*/) {
		std::cout << fillc << std::right << std::setfill(fillh) << std::setw(width + 3) << fillc << "\n";
	}

	void printUsage() {
		const int LINEWIDTH = 60;
		const int DESCPAD = 6;
		IO::printHeader("HOW TO RUN", LINEWIDTH, '-', '|', '+');
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("./CondDeMix -m bulkfile -b referencefile -s sigfile -c sigName", LINEWIDTH, ' ', 1);
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("e.g.: ./CondDeMix -m input/M.tsv -b input/B.tsv -s input/S.tsv", LINEWIDTH, ' ', 1);
		IO::printHeader("List of parameters", LINEWIDTH, '-', '|', '+');
		IO::printRow("-m (required)", LINEWIDTH, ' ');
		IO::printRow("The bulk expression file, in form of a tab-separated data matrix (genes x samples). Data is organized so that the first row is the header row. The header row contains a name for the first column, which is followed by all the sample names. Every subsequent row begins with the name of a gene, then followed by the expression of that gene in all the samples listed in the header.", LINEWIDTH, ' ', DESCPAD);
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("-b (required)", LINEWIDTH, ' ');
		IO::printRow("The single-cell reference file, in form of a tab-separated data matrix (genes x cell types). Data is organized so that the first row is the header row. The header row contains a name for the first column, which is followed by all the cell-type names. Every subsequent row begins with the name of a gene, then followed by the counts of that gene in all the cell types listed in the header.", LINEWIDTH, ' ', DESCPAD);
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("-s (required)", LINEWIDTH, ' ');
		IO::printRow("The mutational signature file, in form of a tab-separated data matrix (samples x signature) with two columns. Data is organized so that the first row is the header row. The header row contains a name for the first column, which is followed by the mutational signature name. Every subsequent row begins with the name of a sample, then followed by the signature strengths in that sample for the mutational signatures listed in the header row.", LINEWIDTH, ' ', DESCPAD);
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("-c (optional)", LINEWIDTH, ' ');
		IO::printRow("A string indicating the name of the signature column (matching the name in the header of the signature matrix -s) to be used in the model.", LINEWIDTH, ' ', DESCPAD);
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("-f (optional; default value is 'output/testRun')", LINEWIDTH, ' ');
		IO::printRow("The path of the output folder, in which all the output files are placed.", LINEWIDTH, ' ', DESCPAD);
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("-v (optional; default value is 0)", LINEWIDTH, ' ');
		IO::printRow("A binary flag (1 or 0) for the verbose output mode. Verbose output mode prints all the output files after each iteration, instead of just the final result. Enabling this flag may use a lot of space to store all the resulting matrices.", LINEWIDTH, ' ', DESCPAD);
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("-t (optional; default value is 0)", LINEWIDTH, ' ');
		IO::printRow("A binary flag (1 or 0) for the synthetic-data simulation mode. In synthetic-data simulation mode, random input data is generated for the cell type fractions matrix, the gamma matrix, and the mutational signature strengths, while the original single-cell reference matrix is kept. The bulk expression is then calculated according to our model, and the accuracy of the recovery of the original fractions and gamma values is measured.", LINEWIDTH, ' ', DESCPAD);
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("-r (optional; default value is 0)", LINEWIDTH, ' ');
		IO::printRow("A binary flag (1 or 0) for the random initialization mode. In random initialization mode, the Gamma matrix is initialized to random values instead of all-zeroes.", LINEWIDTH, ' ', DESCPAD);
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("-l (optional; default value is 0)", LINEWIDTH, ' ');
		IO::printRow("A parameter that denotes the percentage of samples to randomly redraw with repetitions.", LINEWIDTH, ' ', DESCPAD);
		IO::printRow(" ", LINEWIDTH, ' ');
		IO::printRow("-k (optional; by default not used)", LINEWIDTH, ' ');
		IO::printRow("Name of the file containing a list of marker genes, each on a separate line.", LINEWIDTH, ' ', DESCPAD);
		IO::printBreak(LINEWIDTH, '-', '+');
	}

	void readInputParameters(const int argc, char * argv[]) {
		if (argc <= 1) {
			IO::printUsage();
			exit(0);
		}
		const char consoleFlags[] = {'m', 's', 'b', 'f', 't', 'v', 'r', 'c', 'l', 'k', 'p', 'x', 'y'};
		const int numFlags = sizeof(consoleFlags) / sizeof(consoleFlags[0]);
		bool optional[200] = {};
		optional['f'] = true;
		optional['v'] = true;
		optional['t'] = true;
		optional['r'] = true;
		optional['l'] = true;
		optional['k'] = true;
		optional['p'] = true;
		optional['c'] = true;
		optional['x'] = true;
		optional['y'] = true;
		for (int i = 1; i < argc; i++) {
			if ( argv[i][0] == '-' && argv[i][1] ) {
				IO::consoleParameters[ argv[i][1] ] = std::string( argv[i + 1] );
				i++;
			}
		}
		for (int i = 0; i < numFlags; i++) {
			const char flag = consoleFlags[i];
			if ( !optional[flag] && !IO::consoleParameters.count(flag) ) {
				std::cerr << "\n<ERROR in IO::readInputParameters> Missing value for the required parameter '" << flag << "'. Exiting the program.\n";
				exit(0);
			}
		}

		if (!IO::consoleParameters.count('f')) {
			IO::consoleParameters['f'] = "output/testRun";
			std::cout << "Value for parameter 'f' not specified. Defaulting to writing all output to the \"./output/testRun\" folder.\n";
		}

		if (IO::consoleParameters.count('v') && IO::consoleParameters['v'] == "1") {
			std::cout << "<PARAM v=1> Running the program in verbose output mode. (WARNING: Can be hard drive space intensive)\n";
		}
		else {
			IO::consoleParameters['v'] = "0";
			std::cout << "Value for parameter 'v' not specified or invalid. Defaulting to not using the verbose output mode (only the matrices of the final solution will be saved, rather than every step of the iterative algortihm).\n";
		}

		if (IO::consoleParameters.count('t') && IO::consoleParameters['t'] == "1") {
			std::cout << "<PARAM t=1> Running in synthetic-data simulation mode.\n";
			if (!IO::consoleParameters.count('x')) {
				std::cerr << "\n<ERROR in IO::readInputParameters> Parameter 'x' not supplied in synthetic data generation mode. Exiting.\n";
				exit(0);
			}
			if (!IO::consoleParameters.count('y')) {
				std::cerr << "\n<ERROR in IO::readInputParameters> Parameter 'y' not supplied in synthetic data generation mode. Exiting.\n";
				exit(0);
			}
		}
		else {
			IO::consoleParameters['t'] = "0";
			std::cout << "Value for parameter 't' not specified or invalid. Defaulting to not running synthetic-data simulations.\n";
		}

		if (IO::consoleParameters.count('r') && IO::consoleParameters['r'] == "1") {
			std::cout << "<PARAM r=1> Running in randomized gamma initialization mode.\n";
		}
		else {
			IO::consoleParameters['r'] = "0";
			std::cout << "Value for parameter 'r' not specified or invalid. Defaulting to initializing gamma and alpha matrices with zeroes.\n";
		}

		if (IO::consoleParameters.count('l')) {
			long double frac;
			std::stringstream ss(IO::consoleParameters['l']);
			ss >> frac;
			if (frac > 0 && frac <= 1) {
				std::cout << "<PARAM l=" << frac << "> Running in redrawn samples mode.\n";
			}
			else {
				IO::consoleParameters['l'] = "0";
				std::cout << "Value for parameter 'l' not specified or invalid. Defaulting to no redrawing of samples.\n";
			}
		}
		else {
			IO::consoleParameters['l'] = "0";
			std::cout << "Value for parameter 'l' not specified or invalid. Defaulting to no redrawing of samples.\n";
		}

		IO::outputFolder = IO::consoleParameters['f'] + "/";
		std::cout << "Using the output folder '" << IO::outputFolder << "'. All existing data in that folder will be erased.\n";
		system(("rm -rI " + IO::outputFolder).c_str());
		system(("mkdir -pv " + IO::outputFolder).c_str());
	}

	void readMatrix(const std::string & filename, DataMatrix & M) {
		M = DataMatrix();
		int numCols = 0;
		int numRows = 0;
		char cc;
		std::string temp;
		std::ifstream fin(filename.c_str());
		if (!fin) {
			std::cerr << "\n<ERROR in IO::readMatrix> Cannot open [" << filename << "] for reading.\n";
			exit(0);
		}
		// Determining the number of columns
		std::getline(fin, temp);
		for (int i = 0; i < temp.length(); i++) {
			if (temp[i] == '\t') { numCols++; }
		}
		M.colNames.resize(numCols);
		// Read the column names
		{
			bool isFirst = true;
			std::string colName;
			std::stringstream ss(temp);
			int j = 0;
			while (std::getline(ss, colName, '\t')) {
				if (isFirst) {	// discarding the first cell of the header
					isFirst = false;
					continue;
				}
				if (M.colNameToIdx.count(colName)) {
					std::cerr << "\n<ERROR in IO::readMatrix> [" << colName << "] is repeated multiple times in the header row of [" << filename << "]. I am assuming that this is not intended and I am stopping the program execution.\n";
					exit(0);
				}
				M.colNameToIdx[colName] = j;
				M.colNames[j] = colName;
				j++;
			}
		}

		// Determine the number of rows
		while (fin >> temp) {
			numRows++;
			while (!fin.eof() && fin.get() != '\n');
		}
		// std::cerr << numRows << " rows.\n";
		fin.clear();
		fin.seekg(0, fin.beg);
		// Read their names
		while (!fin.eof() && fin.get() != '\n'); // Skip the header row
		for (int i = 0; i < numRows; i++) {
			std::getline(fin, temp);
			// std::cerr << temp << '\n';
			for (int j = 0; j < temp.length(); j++) {
				if (temp[j] == '\t') {
					temp = temp.substr(0, j);
					break;
				}
			}
			// std::cerr << temp << " ";
			if (M.rowNameToIdx.count(temp)) {
				std::cerr << "\n<ERROR in IO::readMatrix> [" << temp << "] is repeated multiple times in the first column of [" << filename << "]. I am assuming that this is not intended and I am stopping the program execution.\n";
				exit(0);
			}
			M.rowNameToIdx[temp] = i;
			M.rowNames.push_back(temp);
			// while (!fin.eof() && fin.get() != '\n');
		}
		// Allocate the space for the expression matrix
		M.val.resize(M.numRows());
		for (int i = 0; i < M.numRows(); i++) {
			M.val[i].resize(M.numCols());
		}
		// Read the expression matrix.
		fin.clear();
		fin.seekg(0, fin.beg);
		std::getline(fin, temp); // Skip the header row
		for (int i = 0; i < M.numRows(); i++) {
			// fin >> temp; // Skip the header column
			while (!fin.eof() && fin.get() != '\t'); // Skip the header column
			// std::cerr << temp;
			for (int j = 0; j < M.numCols(); j++) {
				fin >> M.val[i][j];
				// std::cerr << "\t" << M.val[i][j];
			}
			// std::cerr << "\n";
		}
		fin.close();
		std::cout << "Read the data matrix in " << M.numRows() << " rows and " << M.numCols() << " columns from [" << filename << "].\n";
	}

	void readStringSet(const std::string & filename, std::set<std::string> & S) {
		std::ifstream fin(filename.c_str());
		if (!fin) {
			std::cerr << "\n<ERROR in IO::readMatrix> Cannot open [" << filename << "] for reading.\n";
			exit(0);
		}
		else {
			for (std::string line; std::getline(fin, line);) {
				S.insert(line);
			}
			fin.close();
			std::cout << "Read " << S.size() << " marker names from [" << filename << "].\n";
		}
	}

	void readCellTypeMarkerEntries(const std::string & filename, std::vector< CellTypeMarkerEntry > & V) {
		std::ifstream fin(filename.c_str());
		if (!fin) {
			std::cerr << "\n<ERROR in IO::readMatrix> Cannot open [" << filename << "] for reading.\n";
			exit(0);
		}
		else {
			// struct CellTypeMarkerEntry {
			// 	std::string cellTypeName;
			// 	std::string geneName;
			// 	long double distance;
			// 	int numCellTypesGeneIsExpressedIn;
			int tabIdx[3];
			{ std::string line; std::getline(fin, line); }	// get rid of the first line
			for (std::string line; std::getline(fin, line);) {
				CellTypeMarkerEntry temp;
				const int len = (int)line.length();
				for (int i = 0, idx = 0; i < len; i++) {
					if (line[i] == '\t') {
						tabIdx[idx] = i;
						idx++;
					}
				}
				temp.cellTypeName = line.substr(0, tabIdx[0]);
				temp.geneName = line.substr(tabIdx[0] + 1, tabIdx[1] - tabIdx[0] - 1);
				sscanf(line.c_str() + tabIdx[1] + 1, "%Lf", &temp.distance);
				sscanf(line.c_str() + tabIdx[2] + 1, "%d", &temp.numCellTypesGeneIsExpressedIn);
				V.push_back(temp);
				// fprintf(stderr, "%s\t%s\t%Lf\t%d\n", temp.cellTypeName.c_str(), temp.geneName.c_str(), temp.distance, temp.numCellTypesGeneIsExpressedIn);
			}
			fin.close();
			std::cout << "Read " << V.size() << " markers from [" << filename << "].\n";
		}
	}

	void printMatrix(const std::string & filename, const DataMatrix & M, const std::string & txt /*= "NAME"*/, const DataMatrix * const B /*= nullptr*/) {
		std::ofstream fout(filename.c_str());
		if (!fout) {
			std::cerr << "\n<ERROR in IO::printMatrix> Cannot open [" << filename << "] for writing.\n";
			exit(0);
		}
		fout << txt;
		for (int j = 0; j < M.numCols(); j++) {
			fout << "\t" << M.colNames[j];
		}
		fout << "\n";
		for (int i = 0; i < M.numRows(); i++) {
			fout << M.rowNames[i];
			for (int j = 0; j < M.numCols(); j++) {
				if (B != nullptr && B->val[i][j] <= 1)
					fout << "\tNaN";
				else
					fout << "\t" << M.val[i][j];
			}
			fout << "\n";
		}
		fout.close();
		std::cout << filename << " written.\n";
	}

	void printSyntheticInput(const DataMatrix & B, const DataMatrix & F, const DataMatrix & M, const DataMatrix & gamma, const DataMatrix & alpha, const DataMatrix & S) {
		std::string outFolder;
		{
			std::stringstream ss;
			ss << IO::outputFolder << "syntheticData/";
			system(("rm -rf " + ss.str()).c_str());
			system(("mkdir -pv " + ss.str()).c_str());
			ss >> outFolder;
		}
		{
			std::stringstream ss;
			ss << outFolder << "F.tsv";
			printMatrix(ss.str(), F, "cellType\\sample");
		}
		{
			std::stringstream ss;
			ss << outFolder << "gamma.tsv";
			printMatrix(ss.str(), gamma, "gene\\cellType", &B);
		}
		{
			std::stringstream ss;
			ss << outFolder << "alpha.tsv";
			printMatrix(ss.str(), alpha, "gene\\cellType", &B);
		}
		{
			std::stringstream ss;
			ss << outFolder << "M.tsv";
			printMatrix(ss.str(), M, "gene\\sample");
		}
		{
			std::stringstream ss;
			ss << outFolder << "S.tsv";
			printMatrix(ss.str(), S, "sample\\sig");
		}
		std::cout << "Generated matrices M, S, F, gamma, and alpha have been written to the '" << outFolder << "' folder.\n";
	}

	void printSolution(const int iter, const DataMatrix & B, const DataMatrix & F, const DataMatrix & M, const DataMatrix & gamma, const DataMatrix & alpha, const DataMatrix & S, const std::string & sigName, const bool onlyBasic /* = false*/) {
		const int sigIdx = S.colNameToIdx.find(sigName)->second;
		std::string itFolder;
		{
			std::stringstream ss;
			ss << IO::outputFolder << "it" << std::setfill('0') << std::setw(4) << std::right << iter << "/";
			system(("rm -rf " + ss.str()).c_str());
			system(("mkdir -pv " + ss.str()).c_str());
			ss >> itFolder;
		}
		{
			std::stringstream ss;
			ss << itFolder << "F.tsv";
			printMatrix(ss.str(), F, "cellType\\sample");
		}
		{
			std::stringstream ss;
			ss << itFolder << "gamma.tsv";
			printMatrix(ss.str(), gamma, "gene\\cellType", &B);
		}
		{
			std::stringstream ss;
			ss << itFolder << "alpha.tsv";
			printMatrix(ss.str(), alpha, "gene\\cellType", &B);
		}

		if (onlyBasic) return;
		
		std::ofstream fout, fout2, fout3;
		{
			std::stringstream ss;
			ss << itFolder << "modelExpression.tsv";
			fout.open(ss.str().c_str());
			if (!fout) {
				std::cerr << "\n<ERROR in IO::printSolution> Cannot open [" << ss.str() << "] for writing.\n";
				exit(0);
			}
		}
		{
			std::stringstream ss;
			ss << itFolder << "bulkVsModelExpression.tsv";
			fout2.open(ss.str().c_str());
			if (!fout) {
				std::cerr << "\n<ERROR in IO::printSolution> Cannot open [" << ss.str() << "] for writing.\n";
				exit(0);
			}
		}
		{
			std::stringstream ss;
			ss << itFolder << "error.tsv";
			fout3.open(ss.str().c_str());
			if (!fout) {
				std::cerr << "\n<ERROR in IO::printSolution> Cannot open [" << ss.str() << "] for writing.\n";
				exit(0);
			}
		}
		fout << "gene\\sample";
		fout2 << "gene\\sample";
		fout3 << "gene\ttotalError_L2\taverageBulk\taverageModel\n";
		for (int j = 0; j < M.numCols(); j++) {
			fout << "\t" << M.colNames[j];
			fout2 << "\t" << M.colNames[j] << " bulk\t" << M.colNames[j] << " model";
		}
		fout << "\n";
		fout2 << "\n";
		long double sumOfSquaresError = 0;
		for (int i = 0; i < M.numRows(); i++) {
			fout << M.rowNames[i];
			fout2 << M.rowNames[i];
			fout3 << M.rowNames[i];
			long double totalBulkVal = 0;
			long double totalModelVal = 0;
			long double geneSumOfSquares = 0;
			for (int j = 0; j < M.numCols(); j++) {
				const std::string & sampleName = M.colNames[j];
				if (S.rowNameToIdx.count(sampleName)) {
					const int idx = S.rowNameToIdx.find(sampleName)->second;
					const long double Sj = S.val[idx][sigIdx];
					long double modelVal = 0;
					for (int k = 0; k < B.numCols(); k++) {
						modelVal += F.val[k][j] * (B.val[i][k] + Sj * gamma.val[i][k] + alpha.val[i][k]);
					}
					fout << "\t" << modelVal;
					long double err = M.val[i][j] - modelVal;
					fout2 << "\t" << M.val[i][j] << "\t" << modelVal;
					sumOfSquaresError += err * err;
					geneSumOfSquares += err * err;
					totalModelVal += modelVal;
					totalBulkVal += M.val[i][j];
				}
				else {
					std::cerr << "<ERROR in IO::printSolution> Could not find sample name " << sampleName << " in the signature matrix. Exiting.\n";
					exit(0);
				}
			}
			fout << "\n";
			fout2 << "\n";
			fout3 << "\t" << geneSumOfSquares << "\t" << (totalBulkVal / M.numCols()) << "\t" << (totalModelVal / M.numCols()) << "\n";
		}
		fout.close();
		fout2.close();
		fout3.close();

		{
			std::stringstream ss;
			ss << itFolder << "totalError.txt";
			std::ofstream fout4(ss.str().c_str()); 
			if (!fout4) {
				std::cerr << "\n<ERROR in IO::printSolution> Cannot open [" << ss.str() << "] for writing.\n";
				exit(0);
			}
			else {
				fout4 << "Sum of squares error: " << sumOfSquaresError << "\n";
				std::cout << "Sum of squares error: " << sumOfSquaresError << "\n";
				fout4.close();
			}
		}
	}

	void printDone() {
		std::stringstream ss;
		ss << IO::outputFolder << "done.txt";
		std::ofstream fout(ss.str().c_str()); 
		if (!fout) {
			std::cerr << "\n<ERROR in IO::printDone> Cannot open [" << ss.str() << "] for writing.\n";
			exit(0);
		}
		else {
			fout << "Program has ended successfully.\n\n";
			fout.close();
		}
	}

	void printStats(const int iter, const long double lastL2Error, const long double objVal) {
		const int LINEWIDTH = 40;
		const int LEFTWIDTH = 16;
		const int RIGHTWIDTH = 24;
		// IO::printHeader("Iteration " + std::to_string(iter) + " solved", LINEWIDTH, '-', '|', '+');
		IO::printBreak(LINEWIDTH, '-', '+');
		{
			std::stringstream ss;
			ss << std::right << std::setfill(' ') << std::setw(LEFTWIDTH) << "Last error:" << std::setw(RIGHTWIDTH) << lastL2Error;
			std::cout << "| " << ss.str() << " |\n";
		}
		{
			std::stringstream ss;
			ss << std::right << std::setfill(' ') << std::setw(LEFTWIDTH) << "Objective val:" << std::setw(RIGHTWIDTH) << objVal;
			std::cout << "| " << ss.str() << " |\n";
		}
		{
			double long ratio;
			if (lastL2Error == -1) ratio = 0;
			else ratio = objVal / lastL2Error;
			std::stringstream ss;
			ss << std::right << std::setfill(' ') << std::setw(LEFTWIDTH) << "Ratio:" << std::setw(RIGHTWIDTH) << ratio;
			std::cout << "| " << ss.str() << " |\n";
		}
		IO::printBreak(LINEWIDTH, '-', '+');
	}

	void printHeatmapInput(const int iter, const DataMatrix & F, const DataMatrix & S) {
		const int numSignatures = S.numCols();
		const int numCellTypes = F.numRows();
		std::ofstream fout;
		std::string itFolder;
		{
			std::stringstream ss;
			ss << IO::outputFolder << "it" << std::setfill('0') << std::setw(4) << std::right << iter << "/";
			ss >> itFolder;
		}

		// all cell types
		{
			std::stringstream ss;
			ss << itFolder << "F_forHeatmap_allTypes.txt";
			fout.open(ss.str().c_str());
			if (!fout) {
				std::cerr << "\n<ERROR in IO::printHeatmapInput> Cannot open [" << ss.str() << "] for writing.\n";
				exit(0);
			}
		}
		// header (cell types and signatures)
		fout << "cellType\\sample";
		for (int j = 0; j < numCellTypes; j++) {
			fout << "\t" << F.rowNames[j];
		}
		for (int j = 0 ; j < numSignatures; j++) {
			fout << "\t" << S.colNames[j];
		}
		fout << "\n";
		// data
		for (int i = 0; i < F.numCols(); i++) {
			const std::string & sampleName = F.colNames[i];
			// const std::string & sample_sanitized = sample.substr(0, sample.find("|"));
			fout << sampleName;
			for (int j = 0; j < numCellTypes; j++) {
				const std::string & cellTypeName = F.rowNames[j];
				if (!F.rowNameToIdx.count(cellTypeName)) {
					std::cerr << "<ERROR in IO::printHeatmapInput> Cannot find cell type " << cellTypeName << " in the fractions matrix. Exiting. \n";
					exit(0);
				}
				fout << "\t" << F.getVal(cellTypeName, sampleName);	// rescaled fractions
			}
			if (!S.rowNameToIdx.count(sampleName)) {
				std::cerr << "<ERROR in IO::printHeatmapInput> Cannot find sample " << sampleName << " in the Signature matrix. Exiting. \n";
				exit(0);
			}
			const int sampleIdxInS = S.getRowIdx(sampleName);
			for (int j = 0 ; j < numSignatures; j++) {
				fout << "\t" << S.val[sampleIdxInS][j];
			}
			fout << "\n";
		}
		fout.close();
		
		{
			std::stringstream ss;
			ss << itFolder << "makePlot_allTypes.py";
			fout.open(ss.str().c_str());
			if (!fout) {
				std::cerr << "\n<ERROR in IO::printHeatmapInput> Cannot open [" << ss.str() << "] for writing.\n";
				exit(0);
			}
		}
		fout << "import pandas as pd\n";
		fout << "import seaborn as sns\n";
		fout << "import numpy as np\n";
		fout << "from matplotlib import pyplot as plt\n";
		fout << "sns.set(font_scale=1.5)\n";
		fout << "figTitle = \"Spearman correlation (rescaled fractions)\"\n";
		fout << "inputFilename = \"F_forHeatmap_allTypes.txt\"\n";
		fout << "outputFilename = \"F_forHeatmap_allTypes\"\n";
		fout << "df = pd.read_csv(inputFilename, sep='\\t', header = 0, index_col = 0)\n";
		fout << "newDF = df.corr(method = \"spearman\")\n";
		fout << "newDF = newDF.iloc[:-" << numSignatures << ",-" << numSignatures << ":].iloc[:, ::-1]\n";
		fout << "newDF.to_csv(outputFilename + \"_correlation.tsv\", sep = '\t')\n";
		fout << "sns_plot = sns.clustermap(data = newDF, row_cluster=False, col_cluster=False, cmap = \"vlag\", vmin = -0.6, vmax = 0.6, figsize=(18, 12), annot=True)\n";
		fout << "sns_plot.fig.suptitle(figTitle, x = 0.46, y = 0.9, fontsize = 36)\n";
		fout << "sns_plot.savefig(outputFilename, dpi=150)\n";
		fout.close();
	}
}