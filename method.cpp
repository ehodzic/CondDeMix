#include "method.h"

std::map< int, std::string > GUROBI_OPT_STATUS_CODES = {
	{1, "Model is loaded, but no solution information is available."},
	{2, "Model was solved to optimality (subject to tolerances), and an optimal solution is available."},
	{3, "Model was proven to be infeasible."},
	{4, "Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize."},
	{5, "Model was proven to be unbounded. Important note: an unbounded status indicates the presence of an unbounded ray that allows the objective to improve without limit. It says nothing about whether the model has a feasible solution. If you require information on feasibility, you should set the objective to zero and reoptimize."},
	{6, "Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. No solution information is available."},
	{7, "Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter."},
	{8, "Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter."},
	{9, "Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter."},
	{10, "Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter."},
	{11, "Optimization was terminated by the user."},
	{12, "Optimization was terminated due to unrecoverable numerical difficulties."},
	{13, "Unable to satisfy optimality tolerances; a sub-optimal solution is available."},
	{14, "An asynchronous optimization call was made, but the associated optimization run is not yet complete."},
	{15, "User specified an objective limit (a bound on either the best objective or the best bound), and that limit has been reached."},
	{16, "Optimization terminated because the work expended exceeded the value specified in the WorkLimit parameter."},
};

void runMethod(const DataMatrix & M, const DataMatrix & B, const DataMatrix & S, const std::string & sigName, const std::vector< CellTypeMarkerEntry > & markers, DataMatrix & F, DataMatrix & gamma, DataMatrix & alpha) {
	bool verboseOutput = (IO::consoleParameters['v'] == "1");
	const long double STOPPING_RATIO = 0.97;
	const bool printOnlyBasicOutput = true;
	const int maxNonZeroGammas = 1;
	long double lastL2Error = -1;
	int iter = 0;
	int numJumps = 0;
	DataMatrix bestF, bestAlpha, bestGamma;
	long double bestObj = -1;

	std::set< std::string > markerGenes;
	for (const auto & entry : markers) {
		markerGenes.insert(entry.geneName);
	}
	while (true) {
		iter++;
		long double objVal = 0;
		IO::printHeader("Iteration " + std::to_string(iter), 40, '-', '|', '+');
		solveLP(M, B, S, sigName, markerGenes, F, gamma, alpha, objVal, maxNonZeroGammas, false);
		IO::printStats(iter, lastL2Error, objVal);
		if (iter == 1 || verboseOutput == true) {
			IO::printSolution(iter, B, F, M, gamma, alpha, S, sigName, printOnlyBasicOutput);
			IO::printHeatmapInput(iter, F, S);
		}
		if (bestObj == -1 || objVal < bestObj) {
			bestF = F;
			bestAlpha = alpha;
			bestGamma = gamma;
			bestObj = objVal;
			numJumps = 0;
		}
		if (lastL2Error != -1 && lastL2Error < objVal) numJumps++;
		if (	// Ending optimization if
			(objVal < MIN_VAL_CUTOFF)	// found an optimal solution
			|| (numJumps == MAX_NUM_JUMPS)	// exceeded the allowed number of successive objVal increases
			|| (iter > 4 && lastL2Error >= objVal && objVal / lastL2Error >= STOPPING_RATIO)	// though still decreasing, convergence slowed down too much
		) {
			if (!verboseOutput) {
				IO::printSolution(iter, B, bestF, M, bestGamma, bestAlpha, S, sigName, printOnlyBasicOutput);
				IO::printHeatmapInput(iter, bestF, S);
			}
			break;
		}
		lastL2Error = objVal;
		// break;
	}
	if (!markerGenes.empty()) {
		std::cout << "Cell type proportions have been computed using only marker genes.\n";
		std::cout << "Now extrapolating alphas and gammas for the rest of the genes from those fractions.\n";
		iter++;
		long double objVal = 0;
		IO::printHeader("Iteration " + std::to_string(iter), 40, '-', '|', '+');
		solveLP(M, B, S, sigName, markerGenes, F, gamma, alpha, objVal, maxNonZeroGammas, true);
		IO::printSolution(iter, B, F, M, gamma, alpha, S, sigName, printOnlyBasicOutput);
		IO::printHeatmapInput(iter, F, S);
	}
	IO::printDone();
}

void solveLP(const DataMatrix & M, const DataMatrix & B, const DataMatrix & S, const std::string & sigName, const std::set< std::string > & markerGenes, DataMatrix & solF, DataMatrix & solGamma, DataMatrix & solAlpha, long double & objVal, const int maxNonZeroGammas, const bool extrapolateFromMarkers /*= false*/) {
	const std::vector< std::string > & samples(M.colNames);
	const std::vector< std::string > & genes(M.rowNames);
	const std::vector< std::string > & cellTypes(B.colNames);
	const int numSamples = (int)samples.size();
	const int numGenes = (int)genes.size();
	const int numCellTypes = (int)cellTypes.size();

	// Check that each data matrix has the samples, genes and cell types that we want to us
	for (const std::string & sampleName : samples) {
		std::string evilMatrix = "";
		if (!M.colNameToIdx.count(sampleName)) evilMatrix = "bulk";
		if (!solF.colNameToIdx.count(sampleName)) evilMatrix = "proportion";
		if (!S.rowNameToIdx.count(sampleName)) evilMatrix = "signature";
		if (evilMatrix.length() > 0) {
			std::cerr << "\n\n<ERROR in solveLP> Sample '" << sampleName << "' not found in the " << evilMatrix << " matrix. Exiting.\n";
			exit(0);
		}
	}
	for (const std::string & geneName : genes) {
		std::string evilMatrix = "";
		if (!M.rowNameToIdx.count(geneName)) evilMatrix = "bulk";
		if (!B.rowNameToIdx.count(geneName)) evilMatrix = "reference";
		if (evilMatrix.length() > 0) {
			std::cerr << "\n\n<ERROR in solveLP> Gene '" << geneName << "' not found in the " << evilMatrix << " matrix. Exiting.\n";
			exit(0);
		}
	}
	for (const std::string & cellTypeName : cellTypes) {
		std::string evilMatrix = "";
		if (!solF.rowNameToIdx.count(cellTypeName)) evilMatrix = "proportion";
		if (!B.colNameToIdx.count(cellTypeName)) evilMatrix = "reference";
		if (evilMatrix.length() > 0) {
			std::cerr << "\n\n<ERROR in solveLP> Cell type '" << cellTypeName << "' not found in the " << evilMatrix << " matrix. Exiting.\n";
			exit(0);
		}
	}

	// B_{numGenes x numCellTypes}, F_{numCellTypes x numSamples}, M_{numGenes x numSamples}
	const int sigIdx = S.colNameToIdx.find(sigName)->second;
	if (!extrapolateFromMarkers) {
	try {
		GRBEnv * env = new GRBEnv(true);
		env->set(GRB_IntParam_OutputFlag, 0);
		// env->set(GRB_IntParam_TSPort, replaceThisIfNeeded);
		// env->set(GRB_StringParam_UserName, "replaceThisIfNeeded");
		// env->set(GRB_StringParam_TokenServer, "replaceThisIfNeeded");
		// env->set(GRB_IntParam_Method, 1);	// dual simplex (barrier by default in gurobi 10.0)
		env->set(GRB_IntParam_DualReductions, 0);
		env->start();

		IO::printRow("Stage #1: Solving for F", 40, '|', 0);
		std::chrono::steady_clock::time_point timeBegin = std::chrono::steady_clock::now();

		int nextPrint = 0;
		const int nextPrintStep = numSamples * 0.05;
		for (int j = 0; j < numSamples; j++) {
			if (j == nextPrint) {
				std::cout << "\r" << j << "/" << numSamples;
				nextPrint += nextPrintStep;
			}
			
			GRBModel * model = new GRBModel(*env);
			GRBQuadExpr objective = 0;
			GRBVar * Y = model->addVars(numGenes, GRB_CONTINUOUS);
			GRBVar * F = model->addVars(numCellTypes, GRB_CONTINUOUS);
			for (int i = 0; i < numGenes; i++) {
				Y[i].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
			}

			const std::string & sampleName = samples[j];
			const int sampleIdxInS = S.getRowIdx(sampleName);
			const long double Sj = (S.val[sampleIdxInS][sigIdx] < MIN_VAL_CUTOFF) ? 0 : S.val[sampleIdxInS][sigIdx];

			for (int i = 0; i < numGenes; i++) {
				objective += Y[i] * Y[i];
			}
			model->setObjective(objective, GRB_MINIMIZE);

			// error constraints
			for (int i = 0; i < numGenes; i++) {
				const std::string & geneName = genes[i];
				if (markerGenes.empty() || markerGenes.count(geneName)) {
					GRBQuadExpr rhs = 0;
					for (int k = 0; k < numCellTypes; k++) {
						const std::string & cellTypeName = cellTypes[k];
						// rhs += F[k] * (B.val[i][k] + Sj * solGamma.val[i][k] + solAlpha.val[i][k]);
						rhs += F[k] * (B.getVal(geneName, cellTypeName) + Sj * solGamma.getVal(geneName, cellTypeName) + solAlpha.getVal(geneName, cellTypeName));
					}
					// rhs -= M.val[i][j];
					rhs -= M.getVal(geneName, sampleName);
					model->addConstr(Y[i] == rhs);
				}
				else {	// (!markerGenes.empty() && !markerGenes.count(geneName))
					model->addConstr(Y[i] == 0);
				}
			}

			// Sum of fractions
			{
				GRBLinExpr lhs = 0;
				for (int k = 0; k < numCellTypes; k++) {
					lhs += F[k];
					model->addConstr(F[k] <= 1);
					model->addConstr(F[k] >= 0);
				}
				model->addConstr(lhs == 1);
			}

			model->optimize();
			// model->write("model.sol");
			int optimizationStatusCode = model->get(GRB_IntAttr_Status);
			if (optimizationStatusCode != 2) {
				std::cerr << "\nSolution not optimal. Returned status code:\n" << optimizationStatusCode << ": " << GUROBI_OPT_STATUS_CODES[optimizationStatusCode] << "\n\n";
			}

			for (int k = 0; k < numCellTypes; k++) {
				const std::string & cellTypeName = cellTypes[k];
				long double val_f = F[k].get(GRB_DoubleAttr_X);
				if (val_f < MIN_NONZERO_FRACTION) {
					val_f = 0;
				}
				// solF.val[k][j] = val_f;
				solF.setVal(cellTypeName, sampleName, val_f);
			}
			delete [] Y;
			delete [] F;
			delete model;
		}
		delete env;
		std::cout << "\r";

		std::chrono::steady_clock::time_point timeEnd = std::chrono::steady_clock::now();
		IO::printRow("Optimization time duration = " + std::to_string(std::chrono::duration_cast<std::chrono::seconds>(timeEnd - timeBegin).count()) + "[s]", 40, '|', 0);
	}
	catch (GRBException e){
		std::cerr << "Error code = " << e.getErrorCode() << std::endl;
		std::cerr << e.getMessage() << std::endl;
		exit(0);
	}
	catch (...) {
		std::cerr << "Exception during optimization" << std::endl;
		exit(0);
	}
	}
	
	std::vector<int> expressedInCellTypes(numCellTypes);
	try {
		GRBEnv * env = new GRBEnv(true);
		env->set(GRB_IntParam_OutputFlag, 0);
		// env->set(GRB_IntParam_TSPort, replaceIfNeeded);
		// env->set(GRB_StringParam_UserName, "replaceIfNeeded");
		// env->set(GRB_StringParam_TokenServer, "replaceIfNeeded");
		env->set(GRB_IntParam_Method, 0);	// primal simplex
		env->set(GRB_IntParam_DualReductions, 0);
		env->start();

		IO::printRow("Stage #2: Solving for gamma and alpha", 40, '|', 0);
		std::chrono::steady_clock::time_point timeBegin = std::chrono::steady_clock::now();

		int nextPrint = 0;
		const int nextPrintStep = (numGenes/100 >= 1) ? (numGenes/100) : 1;
		for (int i = 0; i < numGenes; i++) {
			if (i == nextPrint) {
				std::cerr << "\r" << i << "/" << numGenes;
				nextPrint += nextPrintStep;
			}

			//////////////////
			// Second stage //
			//////////////////

			const std::string & geneName = genes[i];
			// Not solving gammas/alphas for this gene if:
			// (i):
			//		There is a set of markers to begin with,
			//		&& We are to extrapolate from markers,
			//		&& this gene is a marker
			// or (ii):
			//		There is a set of markers to begin with,
			//		&& We are NOT to extrapolate from markers
			//		&& this gene is NOT a marker
			if (
				(!markerGenes.empty() && extrapolateFromMarkers && markerGenes.count(geneName))
				|| (!markerGenes.empty() && !extrapolateFromMarkers && !markerGenes.count(geneName))
			) {
				continue;
			}

			if (maxNonZeroGammas == 0) {
				GRBModel * model = new GRBModel(*env);
				GRBQuadExpr objective = 0;
				GRBVar * gamma = model->addVars(numCellTypes, GRB_CONTINUOUS);
				GRBVar * alpha = model->addVars(numCellTypes, GRB_CONTINUOUS);
				GRBVar * Y = model->addVars(numSamples, GRB_CONTINUOUS);
				for (int k = 0; k < numCellTypes; k++) {
					const std::string & cellTypeName = cellTypes[k];
					if (B.getVal(geneName, cellTypeName) >= 1) {
						gamma[k].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
						alpha[k].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
					}
				}
				for (int j = 0; j < numSamples; j++) {
					Y[j].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
				}

				objective = 0;
				for (int j = 0; j < numSamples; j++) {
					objective += Y[j] * Y[j];
				}
				for (int k = 0; k < numCellTypes; k++) {
					const std::string & cellTypeName = cellTypes[k];
					if (B.getVal(geneName, cellTypeName) >= 1) {
						objective += OBJ_COEFF_ALPHA * alpha[k] * alpha[k];
					}
				}

				model->setObjective(objective, GRB_MINIMIZE);

				for (int k = 0; k < numCellTypes; k++) {
					const std::string & cellTypeName = cellTypes[k];
					if (B.getVal(geneName, cellTypeName) >= 1) {
						bool nukeIt = true;
						for (int j = 0; j < numSamples; j++) {
							const std::string & sampleName = samples[j];
							const int sampleIdxInS = S.getRowIdx(sampleName);
							const long double Sj = S.val[sampleIdxInS][sigIdx] < MIN_VAL_CUTOFF ? 0 : S.val[sampleIdxInS][sigIdx];
							if (Sj != 0 && solF.getVal(cellTypeName, sampleName) != 0) {
								nukeIt = false;
								break;
							}
						}
						if (nukeIt) {
							model->addConstr(gamma[k] == 0);
							model->addConstr(alpha[k] == 0);
						}
						else {
						}
					}
				}

				// error
				for (int j = 0; j < numSamples; j++) {
					const std::string & sampleName = samples[j];
					const int sampleIdxInS = S.getRowIdx(sampleName);
					const long double Sj = S.val[sampleIdxInS][sigIdx] < MIN_VAL_CUTOFF ? 0 : S.val[sampleIdxInS][sigIdx];
					GRBLinExpr e = 0;
					for (int k = 0; k < numCellTypes; k++) {
						const std::string & cellTypeName = cellTypes[k];
						if (B.getVal(geneName, cellTypeName) >= 1) {
							e += solF.getVal(cellTypeName, sampleName) * (B.getVal(geneName, cellTypeName) + Sj * gamma[k] + alpha[k]);
							model->addConstr(B.getVal(geneName, cellTypeName) + Sj * gamma[k] + alpha[k] >= 0);
						}
					}
					model->addConstr(Y[j] == M.getVal(geneName, sampleName) - e);
				}

				model->optimize();
				int optimizationStatusCode = model->get(GRB_IntAttr_Status);
				if (optimizationStatusCode != 2) {
					std::cerr << "\nSolution not optimal. Returned status code:\n" << optimizationStatusCode << ": " << GUROBI_OPT_STATUS_CODES[optimizationStatusCode] << "\n\n";
				}

				for (int k = 0; k < numCellTypes; k++) {
					const std::string & cellTypeName = cellTypes[k];
					if (B.getVal(geneName, cellTypeName) >= 1) {
						long double val = gamma[k].get(GRB_DoubleAttr_X);
						if (fabs(val) < MIN_VAL_CUTOFF) {
							val = 0;
						}
						// solGamma.val[i][k] = val;
						solGamma.setVal(geneName, cellTypeName, val);

						val = alpha[k].get(GRB_DoubleAttr_X);
						if (fabs(val) < MIN_VAL_CUTOFF) {
							val = 0;
						}
						// solAlpha.val[i][k] = val;
						solAlpha.setVal(geneName, cellTypeName, val);
					}
				}

				delete [] Y;
				delete [] gamma;
				delete [] alpha;
				delete model;
			}
			else 
			{
				expressedInCellTypes.clear();
				for (int k = 0; k < numCellTypes; k++) {
					const std::string & cellTypeName = cellTypes[k];
					if (B.getVal(geneName, cellTypeName) >= 1) {
						expressedInCellTypes.push_back(k);
					}
					solGamma.setVal(geneName, cellTypeName, 0);
					solAlpha.setVal(geneName, cellTypeName, 0);
				}
				// The code below supports only maxNonZeroGammas = 1 so far
				long double bestSol_error = -1;
				int bestSol_nonZeroIdx = -1;
				long double bestSol_alpha = -1;
				long double bestSol_gamma = -1;				
				for (int idxNonZero = 0; idxNonZero < (int)expressedInCellTypes.size(); idxNonZero++) {
					long double err = solveStageTwo(samples, cellTypes, B, S, M, solF, sigIdx, geneName, expressedInCellTypes, env, solAlpha, solGamma, idxNonZero);
					const int k = expressedInCellTypes[idxNonZero];
					const std::string & cellTypeName = cellTypes[k];
					if (bestSol_error == -1 || err < bestSol_error) {
						bestSol_error = err;
						bestSol_nonZeroIdx = idxNonZero;
						bestSol_gamma = solGamma.getVal(geneName, cellTypeName);
						bestSol_alpha = solAlpha.getVal(geneName, cellTypeName);
					}
					solGamma.setVal(geneName, cellTypeName, 0);
					solAlpha.setVal(geneName, cellTypeName, 0);
				}
				// Set the best solution
				{
					const int k = expressedInCellTypes[bestSol_nonZeroIdx];
					const std::string & cellTypeName = cellTypes[k];
					solAlpha.setVal(geneName, cellTypeName, bestSol_alpha);
					solGamma.setVal(geneName, cellTypeName, bestSol_gamma);
				}
			}
		}/**/
		delete env;
		std::cerr << "\r";


		std::chrono::steady_clock::time_point timeEnd = std::chrono::steady_clock::now();
		IO::printRow("Optimization time duration = " + std::to_string(std::chrono::duration_cast<std::chrono::seconds>(timeEnd - timeBegin).count()) + "[s]", 40, '|', 0);
	}
	catch (GRBException e) {
		std::cerr << "Error code = " << e.getErrorCode() << std::endl;
		std::cerr << e.getMessage() << std::endl;
		exit(0);
	}
	catch (...) {
		std::cerr << "Exception during optimization" << std::endl;
		exit(0);
	}

	// Recompute the sum of squares because of rounding down small values to 0
	objVal = 0;
	for (const std::string & geneName : genes) {
		if (!markerGenes.empty() && !extrapolateFromMarkers && !markerGenes.count(geneName)) continue;

		for (const std::string & sampleName : samples) {
			const int sampleIdxInS = S.getRowIdx(sampleName);
			const long double Sj = (S.val[sampleIdxInS][sigIdx] < MIN_VAL_CUTOFF) ? 0 : S.val[sampleIdxInS][sigIdx];
			
			long double modelVal = 0;
			for (const std::string & cellTypeName : cellTypes ) {
				modelVal += solF.getVal(cellTypeName, sampleName) * (B.getVal(geneName, cellTypeName) + Sj * solGamma.getVal(geneName, cellTypeName) + solAlpha.getVal(geneName, cellTypeName));
			}
			long double err = M.getVal(geneName, sampleName) - modelVal;
			objVal += err * err;
		}

		for (const std::string & cellTypeName : cellTypes ) {
			long double alphaVal = solAlpha.getVal(geneName, cellTypeName);
			objVal += OBJ_COEFF_ALPHA * alphaVal * alphaVal;
		}
	}
}

long double solveStageTwo(const std::vector< std::string > & samples, const std::vector< std::string > cellTypes, const DataMatrix & B, const DataMatrix & S, const DataMatrix & M, const DataMatrix & solF, const int sigIdx, const std::string & geneName, const std::vector<int> & expressedInCellTypes, GRBEnv * env, DataMatrix & solAlpha, DataMatrix & solGamma, const int nonZeroIdx) {
	const int numSamples = (int)samples.size();
	const int numCellTypes = (int)expressedInCellTypes.size();

	const long double EPS_ZERO = 1e-15;
	const int k0 = expressedInCellTypes[nonZeroIdx];
	const std::string & cellTypeName_0 = cellTypes[k0];
	solGamma.setVal(geneName, cellTypeName_0, 0);
	solAlpha.setVal(geneName, cellTypeName_0, 0);
	bool nukeIt = true;
	for (const std::string & sampleName : samples) {
		const int sampleIdxInS = S.getRowIdx(sampleName);
		const long double Sj = (S.val[sampleIdxInS][sigIdx] < MIN_VAL_CUTOFF) ? 0 : S.val[sampleIdxInS][sigIdx];
		if (Sj != 0 && solF.getVal(cellTypeName_0, sampleName) > EPS_ZERO) {
			nukeIt = false;
			break;
		}
	}
	if (nukeIt == false) {
		GRBModel model(*env);
		GRBVar gamma = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
		GRBVar alpha = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
		GRBVar * Y = model.addVars(numSamples, GRB_CONTINUOUS);
		for (int j = 0; j < numSamples; j++) {
			Y[j].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
		}

		GRBQuadExpr objective = 0;
		for (int j = 0; j < numSamples; j++) {
			objective += Y[j] * Y[j];
		}
		objective += OBJ_COEFF_ALPHA * alpha * alpha;
		model.setObjective(objective, GRB_MINIMIZE);

		// error
		for (int j = 0; j < numSamples; j++) {
			const std::string & sampleName = samples[j];
			const int sampleIdxInS = S.getRowIdx(sampleName);
			const long double Sj = (S.val[sampleIdxInS][sigIdx] < MIN_VAL_CUTOFF) ? 0 : S.val[sampleIdxInS][sigIdx];
			// Constr: non-negative modified reference counts
			model.addConstr(B.getVal(geneName, cellTypeName_0) + Sj * gamma + alpha >= 0);
			// Constr: error variable
			long double sumRefFrac = 0;
			for (const int k : expressedInCellTypes) {
				const std::string & cellTypeName = cellTypes[k];
				const long double bk = B.getVal(geneName, cellTypeName);
				const long double fkj = solF.getVal(cellTypeName, sampleName);
				sumRefFrac += bk * fkj;
			}
			const long double fk0j = solF.getVal(cellTypeName_0, sampleName);
			model.addConstr(Y[j] == M.getVal(geneName, sampleName) - sumRefFrac - fk0j*Sj*gamma - fk0j*alpha);
		}

		model.optimize();
		int optimizationStatusCode = model.get(GRB_IntAttr_Status);
		if (optimizationStatusCode != 2) {
			std::cerr << "\nSolution not optimal. Returned status code:\n" << optimizationStatusCode << ": " << GUROBI_OPT_STATUS_CODES[optimizationStatusCode] << "\n\n";
		}

		long double val = gamma.get(GRB_DoubleAttr_X);
		if (fabs(val) < MIN_VAL_CUTOFF) {
			val = 0;
		}
		solGamma.setVal(geneName, cellTypeName_0, val);

		val = alpha.get(GRB_DoubleAttr_X);
		if (fabs(val) < MIN_VAL_CUTOFF) {
			val = 0;
		}
		solAlpha.setVal(geneName, cellTypeName_0, val);

		delete [] Y;
	}

	// Compute the error for this solution
	return getSquareError(samples, cellTypes, B, S, M, solF, sigIdx, geneName, expressedInCellTypes, solAlpha, solGamma);
}

long double getSquareError(const std::vector< std::string > & samples, const std::vector< std::string > cellTypes, const DataMatrix & B, const DataMatrix & S, const DataMatrix & M, const DataMatrix & solF, const int sigIdx, const std::string & geneName, const std::vector<int> & expressedInCellTypes, const DataMatrix & solAlpha, const DataMatrix & solGamma) {
	const int numSamples = (int)samples.size();
	const int numCellTypes = (int)expressedInCellTypes.size();

	long double totalError = 0;
	for (const std::string & sampleName : samples) {
		const int sampleIdxInS = S.getRowIdx(sampleName);
		const long double Sj = (S.val[sampleIdxInS][sigIdx] < MIN_VAL_CUTOFF) ? 0 : S.val[sampleIdxInS][sigIdx];
		
		long double modelVal = 0;
		for (const int k : expressedInCellTypes) {
			const std::string & cellTypeName = cellTypes[k];
			modelVal += solF.getVal(cellTypeName, sampleName) * (B.getVal(geneName, cellTypeName) + Sj * solGamma.getVal(geneName, cellTypeName) + solAlpha.getVal(geneName, cellTypeName));
		}
		long double err = M.getVal(geneName, sampleName) - modelVal;
		totalError += err * err;
	}

	for (const int k : expressedInCellTypes) {
		const std::string & cellTypeName = cellTypes[k];
		
		long double alphaVal = solAlpha.getVal(geneName, cellTypeName);
		totalError += OBJ_COEFF_ALPHA * alphaVal * alphaVal;
	}
	return totalError;
}
