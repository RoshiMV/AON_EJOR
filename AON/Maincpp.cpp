#include <iostream>
#include <vector>
#include <algorithm>
#include <time.h>
#include <chrono>
#include "DataGenerator.h"
#include "TestFunctions.h"


using namespace std;

int main()
{

	typedef std::chrono::high_resolution_clock Time;
	typedef std::chrono::milliseconds ms;
	typedef std::chrono::duration<float> fsec;

	TestFunctions testDataGenerator;
	
	//testDataGenerator.GenerateAndTestProblemsWithKR();
	//testDataGenerator.GenerateAndTestProblemsNoKR();
	testDataGenerator.SenesitivityAnalysisWithSAA();
	//testDataGenerator.UpperBoundEstimateForSAAtest();
	//testDataGenerator.TestHeuristicNoKR();

	system("pause");
	return 0;
}