#ifndef TESTFUNCTIONS_H
#define TESTFUNCTIONS_H
#include "DataGenerator.h"
#include "ProblemInstance.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <math.h>  
using namespace std;


class TestFunctions {
public:

	void CallDataGenerator() {

		DataGenerator DG;
		cout << "Generating random problem instances . . . " << endl;
		DG.GenerateRandomInstances("Input.csv");
		cout << "Writing problem instance files . . . " << endl;
		DG.WriteRandomInstances("Instances");

	}

	void GenerateAndTestProblemsWithKR() {

		DataGenerator DG;

		// Generating Random Instances
		cout << "Generating random problem instances . . . " << endl;
		DG.GenerateRandomInstances("Input.csv");
		cout << "Writing problem instance files . . . " << endl;
		DG.WriteRandomInstances("InstancesSAA");

		// Generating Exact Models with KR
		DG.WriteExactModelsInGamsWithKRIPOPT("InstancesExactGamsWithKRIPOPT");
		DG.WriteExactModelsInGamsWithKRConopt("InstancesExactGamsWithKRConopt");
		DG.WriteExactModelsInGamsWithKRLindoglobal("InstancesExactGamsWithKRLindoglobal");
		DG.WriteExactGAMSBatchFile("InstancesExactGamsWithKRIPOPT");

		DG.WriteSAAModelsInGams("InstancesSAA");
		DG.WriteSaaGAMSBatchFile("InstancesSAA");
		DG.WriteSaaGAMSBatchBatchFile("InstancesSAA");

		DG.WriteHeuristicModelsGamsWithKR_ApproximateIterative("InstancesHeuristic");
		DG.WriteApproxHeuristicGAMSBatchFile("InstancesHeuristic");
		//DG.WriteHeuristicModelsGamsWithKR_ExactIterative("Instances");
	}

	void UpperBoundEstimateForSAAtest() {

		DataGenerator DG;
		DG.GenerateRandomInstances("Input.csv");
		DG.EstimateStatisticalUpperBoundForAll("Instances");

	}

	void GenerateAndTestProblemsNoKR() {

		DataGenerator DG;

		// Generating Random Instances
		cout << "Generating random problem instances . . . " << endl;
		DG.GenerateRandomInstances("Input.csv");
		cout << "Writing problem instance files . . . " << endl;
		DG.WriteRandomInstances("Instances");

		DG.WriteExactModelsInGamsNoKRIPOPT("InstancesExactGamsNoKRIPOPT");
		//DG.WriteExactModelsInGamsNoKRConopt("InstancesExactGamsNoKRConopt");
		//DG.WriteExactModelsInGamsNoKRLindoglobal("InstancesExactGamsNoKRLindoglobal");

		DG.WriteExactGAMSBatchFileNoKR("InstancesExactGamsNoKRIPOPT");

		DG.WriteSAAModelsInGamsNoKR("InstancesSAA");
		DG.WriteSaaGAMSBatchFile("InstancesSAA");
		DG.WriteSaaGAMSBatchBatchFile("InstancesSAA");

		DG.WriteHeuristicModelsGamsNoKR("InstancesHeuristic");
		DG.WriteNoKRHeuristicGAMSBatchFile("InstancesHeuristic");

	}

	void TestHeuristicNoKR() {

		DataGenerator DG;
		DG.GenerateRandomInstances("Input.csv");
		//DG.WriteExactModelsInGamsNoKR("InstancesExactGamsNoKR");
		DG.WriteHeuristicModelsGamsNoKR("InstancesHeuristicNoKR");

	}

	void SenesitivityAnalysisWithSAA() {

		DataGenerator DG;

		// Generating Random Instances
		cout << "Generating random problem instances . . . " << endl;
		DG.GenerateRandomInstances("Input.csv");
		cout << "Writing problem instance files . . . " << endl;
		DG.WriteRandomInstances("Instances");

		// Generating Exact Models with KR
		DG.WriteExactModelsInGamsWithKRIPOPT("InstancesExactGamsWithKRIPOPT");

		DG.WriteSAAModelsInGams("InstancesSAA");
		DG.WriteSaaGAMSBatchFile("InstancesSAA");
		DG.WriteSaaGAMSBatchBatchFile("InstancesSAA");

		
	}
};


#endif