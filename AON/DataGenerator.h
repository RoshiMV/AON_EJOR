#ifndef DATAGENERATOR_H
#define DATAGENERATOR_H

#include <vector>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h> 

using namespace std;

class ProblemInstance;

class DataGenerator
{
public:
	vector<ProblemInstance*>* Instances;

	int numSuppliers;
	int numScenarios;
	int numSets;
	int numScenariosForUB; // number of scenarios used in calculating the upper bound
	int hCost;
	int rc;
	int ecLB;
	int ecUB;
	double demandMean;
	double demandSD;
	double spotMarketMean;
	double spotMarketSD;
	int ciLB;
	int ciUB;
	double piLB;
	double piUB;
	int numInstances;
	vector<vector<int>> matrixOfScenarios;
	//vector<double> probabilityOfScenarios;

	DataGenerator() { Instances = new vector<ProblemInstance*>; 	srand(time(NULL)); }
	~DataGenerator() { Instances->clear(); }

	double Norm_Standard_CDF(double z_value_In);
	double Norm_Standard_PDF(double z_value_In);
	double Norm_CDF(double Q_value_In, double mean, double standard_dev);
	double Norm_Standard_CDF_INV(double probability_In);
	double Norm_CDF_INV(double probability, double mean, double standard_dev);
	double Poly_value(int n, double a[], double x_In);


	void GenerateRandomInstances(string InputFileName);
	void WriteRandomInstances(string FolderName);
	void EstimateStatisticalUpperBound(string FolderNameIn, string FileName);
	void EstimateStatisticalUpperBoundForAll(string FolderName);

	void WriteExactModelsInGamsWithKRIPOPT(string FolderName);
	void WriteExactModelsInGamsWithKRConopt(string FolderName);
	void WriteExactModelsInGamsWithKRLindoglobal(string FolderName);

	void WriteExactModelsInGamsNoKRIPOPT(string FolderName);
	void WriteExactModelsInGamsNoKRConopt(string FolderName);
	void WriteExactModelsInGamsNoKRLindoglobal(string FolderName);

	void WriteSAAModelsInGams(string FolderName);
	void WriteExactGAMSBatchFile(string FolderName);
	void WriteExactGAMSBatchFileNoKR(string FolderName);
	void WriteSaaGAMSBatchFile(string FolderName);
	void WriteSAAModelsInGamsNoKR(string FolderName);
	void WriteApproxHeuristicGAMSBatchFile(string FolderName);
	void WriteNoKRHeuristicGAMSBatchFile(string FolderName);
	void WriteSaaGAMSBatchBatchFile(string FolderName);

	void WriteHeuristicModelsGamsWithKR(string FolderName);
	void WriteHeuristicModelsGamsWithKR_ApproximateIterative(string FolderName);
	void WriteHeuristicModelsGamsWithKR_ExactIterative(string FolderName);
	void WriteHeuristicModelsGamsNoKR(string FolderName);
	void GenerateAllPossibleScenarios();

	vector<double> CalcScenarioProbabilities(double* instanceProbabilities);
	double CalcScenarioCostOfGivenSolution(ProblemInstance* pIn, double scenariodemand, double scenarioSpotPrice, double QtotalIn, double KR);
private:

	void CreateProblemInstance(int insNo);

};



#endif
