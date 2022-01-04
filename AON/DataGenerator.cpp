#include "DataGenerator.h"
#include "ProblemInstance.h"
#include "Parameters.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <math.h>  
#include <chrono>
using namespace std;

// Distribution functions

double DataGenerator::Norm_Standard_CDF(double z_value_In) {


	//  http://people.sc.fsu.edu/~jburkardt/cpp_src/asa066/asa066.html
	//   Norm_Standard_CDF computes the cumulative density of the standard normal distribution.
	//
	//    Input:
	//		double z_value_In
	//
	//    Output:
	//		double cdf_prob_out
	//		the integrals of the standard normal distribution over the interval ( - Infinity, z_value]

	double a0 = 0.5;
	double a1 = 0.398942280444;
	double a2 = 0.399903438504;
	double a3 = 5.75885480458;
	double a4 = 29.8213557808;
	double a5 = 2.62433121679;
	double a6 = 48.6959930692;
	double a7 = 5.92885724438;
	double b0 = 0.398942280385;
	double b1 = 0.000000038052;
	double b2 = 1.00000615302;
	double b3 = 0.000398064794;
	double b4 = 1.98615381364;
	double b5 = 0.151679116635;
	double b6 = 5.29330324926;
	double b7 = 4.8385912808;
	double b8 = 15.1508972451;
	double b9 = 0.742380924027;
	double b10 = 30.789933034;
	double b11 = 3.99019417011;
	double y;
	double zabs;
	double pdf;
	double q;
	double cdf_prob_out;

	zabs = fabs(z_value_In);
	//  |Z| between 0 and 1.28
	if (zabs <= 1.28)
	{
		y = a0 * z_value_In * z_value_In;
		pdf = exp(-y) * b0;
		q = a0 - zabs * (a1 - a2 * y
			/ (y + a3 - a4
				/ (y + a5 + a6
					/ (y + a7))));
	}
	//  |Z| between 1.28 and 12.7
	else if (zabs <= 12.7)
	{
		y = a0 * z_value_In * z_value_In;
		pdf = exp(-y) * b0;
		q = pdf
			/ (zabs - b1 + b2
				/ (zabs + b3 + b4
					/ (zabs - b5 + b6
						/ (zabs + b7 - b8
							/ (zabs + b9 + b10
								/ (zabs + b11))))));
	}
	//  Z far out in tail.
	else {
		q = 0.0;
		pdf = 0.0;
	}

	if (z_value_In < 0.0) {
		cdf_prob_out = q;
		q = 1.0 - cdf_prob_out;
	}
	else {
		cdf_prob_out = 1.0 - q;
	}

	return cdf_prob_out;
}

double DataGenerator::Norm_Standard_PDF(double z_value_In) {

	//   Norm_Standard_CDF computes the cumulative density of the standard normal distribution.
	//
	//    Input:
	//		double z_value_In
	//
	//    Output: 
	//		double pdf_out, the value of the standard normal distribution at z_value_In
	//
	double a0 = 0.5;
	double b0 = 0.398942280385;
	double y;
	double zabs;
	double pdf_out;

	zabs = fabs(z_value_In);
	//  |Z| between 0 and 1.28
	if (zabs <= 1.28)
	{
		y = a0 * z_value_In * z_value_In;
		pdf_out = exp(-y) * b0;
	}
	//  |Z| between 1.28 and 12.7
	else if (zabs <= 12.7)
	{
		y = a0 * z_value_In * z_value_In;
		pdf_out = exp(-y) * b0;
	}
	//  Z far out in tail.
	else {
		pdf_out = 0.0;
	}

	return pdf_out;
}

double DataGenerator::Norm_CDF(double value_In, double mean, double standard_dev) {
	double temp = (value_In - mean) / standard_dev;
	return Norm_Standard_CDF(temp);
}

double DataGenerator::Norm_Standard_CDF_INV(double probability_In) {

	//    Input:
	//		double probability_In, the value of the cumulative probability densitity function.
	//    0 < P < 1.  If P is outside this range, an "infinite" result is returned.
	//
	//    Output,:
	//		double CDF_INV_val

	Parameters pr;
	static double a[4] = { 3.3871327179, 50.434271938, 159.29113202, 59.109374720 };
	static double b[4] = { 1.0, 17.895169469, 78.757757664, 67.187563600 };
	static double c[4] = { 1.4234372777, 2.7568153900, 1.3067284816, 0.17023821103 };
	static double const1 = 0.180625;
	static double const2 = 1.6;
	static double d[3] = { 1.0, 0.73700164250, 0.12021132975 };
	static double e[4] = { 6.6579051150, 3.0812263860, 0.42868294337, 0.017337203997 };
	static double f[3] = { 1.0, 0.24197894225, 0.012258202635 };
	double q;
	double r;
	static double split1 = 0.425;
	static double split2 = 5.0;
	double value;

	if (probability_In <= 0.0) {
		value = -pr.infinity;
		return value;
	}

	if (1.0 <= probability_In) {
		value = pr.infinity;
		return value;
	}

	q = probability_In - 0.5;

	if (fabs(q) <= split1) {
		r = const1 - q * q;
		value = q * Poly_value(4, a, r) / Poly_value(4, b, r);
	}
	else {
		if (q < 0.0) {
			r = probability_In;
		}
		else
		{
			r = 1.0 - probability_In;
		}

		if (r <= 0.0) {
			value = -1.0;
			exit(1);
		}

		r = sqrt(-log(r));

		if (r <= split2) {
			r = r - const2;
			value = Poly_value(4, c, r) / Poly_value(3, d, r);
		}
		else {
			r = r - split2;
			value = Poly_value(4, e, r) / Poly_value(3, f, r);
		}

		if (q < 0.0) {
			value = -value;
		}

	}

	return value;
}

double DataGenerator::Norm_CDF_INV(double probability, double mean, double standard_dev) {
	double result_tmp = DataGenerator::Norm_Standard_CDF_INV(probability);
	result_tmp = standard_dev* result_tmp + mean;
	return result_tmp;
}

double DataGenerator::Poly_value(int n, double a[], double x_In) {

	//    POLY_VALUE evaluates a real polynomial.
	//		https://people.sc.fsu.edu/~jburkardt/cpp_src/asa241/asa241.cpp
	//
	//      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)
	//
	//    Input: 
	//		int N, the order of the polynomial.
	//      double A[N], the coefficients of the polynomial.
	//		A[0] is the constant term.
	//		double X, the point at which the polynomial is to be evaluated.
	//
	//    Output:
	//		double POLY_VALUE
	//
	int i;
	double value_tmp;
	value_tmp = 0.0;

	for (i = n - 1; 0 <= i; i--)
	{
		value_tmp = value_tmp * x_In + a[i];
	}

	return value_tmp;
}



void DataGenerator::GenerateRandomInstances(string InputFileName)
{
	ifstream input_file(InputFileName.c_str());

	if (!input_file.is_open())
	{
		cout << "File " << InputFileName << " not found!" << endl;
		system("pause");
		return;
	}
	string x;
	getline(input_file, x, '\n');

	while (!input_file.eof())
	{
		getline(input_file, x, '\n');
		if (x.size()<1) continue;

		string token;
		stringstream s(x);
		getline(s, token, ',');
		numSuppliers = atoi(token.c_str());
		getline(s, token, ',');
		numScenarios = atoi(token.c_str());
		getline(s, token, ',');
		numSets = atoi(token.c_str());
		getline(s, token, ',');
		numScenariosForUB = atoi(token.c_str());
		getline(s, token, ',');
		hCost = atoi(token.c_str());
		getline(s, token, ',');
		rc = atof(token.c_str());
		getline(s, token, ',');
		ecLB = atoi(token.c_str());
		getline(s, token, ',');
		ecUB = atoi(token.c_str());
		getline(s, token, ',');
		demandMean = atof(token.c_str());
		getline(s, token, ',');
		demandSD = atof(token.c_str());
		getline(s, token, ',');
		spotMarketMean = atof(token.c_str());
		getline(s, token, ',');
		spotMarketSD = atof(token.c_str());
		getline(s, token, ',');
		ciLB = atof(token.c_str());
		getline(s, token, ',');
		ciUB = atof(token.c_str());
		getline(s, token, ',');
		piLB = atof(token.c_str());
		getline(s, token, ',');
		piUB = atof(token.c_str());
		getline(s, token, ',');
		numInstances = atoi(token.c_str());

		for (int i = 0; i < numInstances; i++)
		{
			CreateProblemInstance(i);
		}

	}

	input_file.close();
}

void DataGenerator::GenerateAllPossibleScenarios(){
	int n = numSuppliers;
	int combinations = pow(2, n);
	int total = 0;
	vector<vector<int>> matrix;
	for (int j = 0; j < combinations; j++) {
		int num = j;
		vector<int> temp_row;
		for (int k = 0; k < n; k++)
		{
			total = num % 2;
			temp_row.insert(temp_row.begin(), total);
			num /= 2;
		}
		matrix.push_back(temp_row);
	}
	matrixOfScenarios = matrix;
}

vector<double> DataGenerator::CalcScenarioProbabilities(double* instanceProbabilitiesIn) {
	int combinations = pow(2, numSuppliers);
	vector<double> result;
	for (int combination_tmp = 0; combination_tmp < combinations; combination_tmp++) {
		double prob_tmp = 1.0;
		for (int k = 0; k < numSuppliers; k++) {
			int test = matrixOfScenarios[combination_tmp][k];
			if (test == 1) {
				prob_tmp *= instanceProbabilitiesIn[k];
			}
			else {
				prob_tmp *= (1.0 - instanceProbabilitiesIn[k]);
			}
		}
		result.push_back(prob_tmp);
	}
	return result;
}

void DataGenerator::CreateProblemInstance(int insnoIn)
{
	ProblemInstance* p = new ProblemInstance(insnoIn, numSuppliers);
	double ranno = (double)rand() / RAND_MAX;
	int temp1 = (ecUB - ecLB)*ranno + ecLB;
	p->AssignExecutionCost(temp1);
	ranno = (double)rand() / RAND_MAX;
	temp1 = hCost; //(50 - hCost)*ranno + hCost;
	p->Assignh(temp1);

	int* cost_sort_tmp = new int[numSuppliers];

	for (int i = 0; i < numSuppliers; i++)
	{
		double ranno2 = (double)rand() / RAND_MAX;
		int temp2 = floor((double)((ciUB - ciLB))*ranno2 + (double)ciLB);	
		cost_sort_tmp[i]= temp2;
	}

	// Sort the cost array in the nonincreasing order of critical fractiles
	double compare_tmp;
	for (int i = 0; i < numSuppliers; i++)
	{
		for (int j = i + 1; j < numSuppliers; j++)
		{
			/*double fraction_c_i = (spotMarketMean - ((double)cost_sort_tmp[i])) / (spotMarketMean + ((double)hCost));
			double fraction_c_j = (spotMarketMean - ((double)cost_sort_tmp[j])) / (spotMarketMean + ((double)hCost));*/
			if (cost_sort_tmp[j] < cost_sort_tmp[i])
			{
				compare_tmp = cost_sort_tmp[i];
				cost_sort_tmp[i] = cost_sort_tmp[j];
				cost_sort_tmp[j] = compare_tmp;
			}
		}
	}
	p->AssignCostArray(cost_sort_tmp);

	//// This part is temporary for sorting the pi values
	double* prob_sort_tmp = new double[numSuppliers];
	for (int i = 0; i < numSuppliers; i++)
	{
		double ranno3 = (double)rand() / RAND_MAX;
		double temp3 = (double)((piUB - piLB)) * ranno3 + (double)piLB;
		prob_sort_tmp[i] = temp3;
	}

	for (int i = 0; i < numSuppliers; i++)
	{
		for (int j = i + 1; j < numSuppliers; j++)
		{
			if (prob_sort_tmp[j] < prob_sort_tmp[i])
			{
				compare_tmp = prob_sort_tmp[i];
				prob_sort_tmp[i] = prob_sort_tmp[j];
				prob_sort_tmp[j] = compare_tmp;
			}
		}
	}
	p->AssignProbArray(prob_sort_tmp);
	//// to this point
	/*for (int i = 0; i < numSuppliers; i++)
	{
		double ranno3 = (double)rand() / RAND_MAX;
		double temp3 = (double)((piUB - piLB))*ranno3 + (double)piLB;
		p->AssignProb(i, temp3);
	}*/

	Instances->push_back(p);
}

void DataGenerator::WriteRandomInstances(string FolderName)
{
	vector<ProblemInstance*>::iterator it;
	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine demandGenerator(seed1);
	std::normal_distribution<double> distribution1(demandMean, demandSD);
	std::default_random_engine priceGenerator(seed2);
	std::normal_distribution<double> distribution2(spotMarketMean, spotMarketSD);

	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		for (int j = 0; j < numSets; j++) {
			stringstream convert_t, convert_ins, convert_set;
			convert_t << numSuppliers;
			convert_ins << (*it)->instanceNo;
			convert_set << j;
			// Write the first inc file for the demand information
			string out_file = FolderName + "/demand_" + convert_t.str() + "_" + convert_ins.str() + "_set_" + convert_set.str() + ".inc";
			// name of each sxcel file is numberOfSuppliers_probInstanceNum_setNum(set of scenarios)
			ofstream output;
			output.open(out_file.c_str());
			output << "parameter data(m, j)/" << endl;
			output << "$ondelim" << endl;
			for (int i = 0; i < numScenarios; i++)
			{
				double ranno_demand = round(distribution1(demandGenerator)*1000.0) / 1000.0;
				output << i + 1 << '\t' << "Demand" << '\t' << ranno_demand << endl;
			}
			output << "$offdelim" << endl;
			output.close();
			// Write the second inc file for the SP information
			out_file = FolderName + "/SP_" + convert_t.str() + "_" + convert_ins.str() + "_set_" + convert_set.str() + ".inc";
			output.open(out_file.c_str());
			output << "$ondelim" << endl;
			for (int i = 0; i < numScenarios; i++)
			{
				double ranno_price = round(distribution2(priceGenerator)*1000.0) / 1000.0;
				output << i + 1 << '\t' << "SP" << '\t' << ranno_price << endl;
			}
			output << "$offdelim" << endl;
			output << "/;" << endl;
			output.close();
			// Write the third inc file for the R data
			out_file = FolderName + "/Rdata_" + convert_t.str() + "_" + convert_ins.str() + "_set_" + convert_set.str() + ".inc";
			output.open(out_file.c_str());
			output << "parameter R(m, i)/" << endl;
			output << "$ondelim" << endl;
			/*for (int k = 0; k < numSuppliers; k++) {
			output << k + 1 << '\t';
			}
			output << endl;*/
			for (int k = 0; k < numSuppliers; k++)
			{
				for (int i = 0; i < numScenarios; i++) {
					output << i + 1 << '\t';
					output << k + 1;
					int availableInScenario = 0;
					double ranno = round(((double)rand() / RAND_MAX)*1000.0) / 1000.0;
					if (ranno < (*it)->prob_tmp[k]) {
						availableInScenario = 1;
					}
					output << '\t' << availableInScenario << endl;
				}
			}
			output << "$offdelim" << endl;
			output << "/;" << endl;
			output.close();

		}

	}

	//for (it = Instances->begin(); it < Instances->end(); it++)
	//{
	//	for (int j = 0; j < numSets; j++) {
	//		stringstream convert_t, convert_ins, convert_set;
	//		convert_t << numSuppliers;
	//		convert_ins << (*it)->instanceNo;
	//		convert_set << j;
	//		string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_" + convert_set.str() + ".csv";
	//		// name of each sxcel file is numberOfSuppliers_probInstanceNum_setNum(set of scenarios)
	//		ofstream output;
	//		output.open(out_file.c_str());
	//		output << " ,Demand,SP,";
	//		for (int k = 0; k < numSuppliers; k++) {
	//			output << "," << k + 1;
	//		}
	//		output << endl;
	//		for (int i = 0; i < numScenarios; i++)
	//		{				
	//			double ranno_demand = round(distribution1(demandGenerator)*1000.0)/ 1000.0;
	//			double ranno_price = round(distribution2(priceGenerator)*1000.0)/ 1000.0;
	//			output << i + 1 << "," << ranno_demand << "," << ranno_price << "," << i+1 ;
	//			for (int k = 0; k < numSuppliers; k++) {
	//				int availableInScenario=0;
	//				double ranno = round(((double)rand() / RAND_MAX)*1000.0)/1000.0;
	//				if (ranno < (*it)->prob_tmp[k]) {
	//					availableInScenario = 1;
	//				}
	//				output << "," << availableInScenario;
	//			}
	//			output << endl;
	//		}
	//		output.close();
	//	}
	//	
	//}

	// write the basic information of each instance in an excel file 
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
			stringstream convert_t, convert_ins, convert_set;
			convert_t << numSuppliers;
			convert_ins << (*it)->instanceNo;
			string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + ".csv";
			ofstream output;
			output.open(out_file.c_str());
			output << "insNo," << (*it)->instanceNo << endl;
			output << "n," << numSuppliers << endl;
			output << "ec," << (*it)->ec_tmp << endl;
			output << "supplierNum,cost,prob,qi" << endl;
			for (int k = 0; k < numSuppliers; k++) {
				output << k + 1 << "," << (*it)->cost_tmp[k] << "," << (*it)->prob_tmp[k] << endl;
			}
			output << "KR" << endl;
			output.close();

	}

}

double DataGenerator::CalcScenarioCostOfGivenSolution(ProblemInstance* pIn, double scenariodemand, double scenarioSpotPrice, double QtotalIn, double KR) {
	double result = 0;
	double o_m;
	double s_m;

	double Iplus = std::fmax(0.0, QtotalIn - scenariodemand);
	double Iminus = std::fmax(0.0, scenariodemand - QtotalIn);
	if (pIn->ec_tmp < scenarioSpotPrice) {
		if (Iminus <= KR) {
			o_m = Iminus;
			s_m = 0.0;
		}
		else {
			o_m = KR;
			s_m = Iminus - KR;
		}
	}
	else {
		o_m = 0.0;
		s_m = Iminus;
	}
	result += ((double)rc) * KR + ((double)hCost)*Iplus + ((double)pIn->ec_tmp)*o_m + scenarioSpotPrice*s_m;
	return result;
}

void DataGenerator::EstimateStatisticalUpperBoundForAll(string FolderName) {

	for (int k = 0; k < numInstances; k++) {
			stringstream convert_t, convert_ins, convert_set;
			convert_t << numSuppliers;
			convert_ins << k;
			string instanceFileName = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + ".csv";
			DataGenerator::EstimateStatisticalUpperBound(FolderName, instanceFileName);
	}

}

void DataGenerator::EstimateStatisticalUpperBound(string FolderNameIn, string FileName) {
	
	int supplierNo;

	int numSuppliersIn=numSuppliers;
	int executiveCostIn;
	int insNumberIn=0;
	int* costIn;
	double* probabilityIn;
	
	double* orderQuantities;
	double Qtotal;
	double KR_In;

	ifstream input_file(FileName.c_str());

	if (!input_file.is_open())
	{
		cout << "File " << FileName << " not found!" << endl;
		system("pause");
		return;
	}
	string x;
	
	getline(input_file, x, '\n');
	
	string token;
	stringstream s(x);
	getline(s, token, ',');
	getline(s, token, '\n');
	insNumberIn = atoi(token.c_str());
	getline(input_file, x, '\n');
	s.clear();
	s.str(x);
	getline(s, token, ',');
	getline(s, token, '\n');
	numSuppliersIn = atoi(token.c_str());
	getline(input_file, x, '\n');
	s.clear();
	s.str(x);
	getline(s, token, ',');
	getline(s, token, ',');
	executiveCostIn = atoi(token.c_str());
	getline(input_file, x, '\n');
	getline(input_file, x, '\n');
	s.clear();
	s.str(x);
	costIn = new int[numSuppliersIn];
	probabilityIn = new double[numSuppliersIn];
	orderQuantities = new double[numSuppliersIn];
	for (int m = 0; m < numSuppliersIn; m++) {
		getline(s, token, ',');
		supplierNo = atoi(token.c_str());
		getline(s, token, ',');
		costIn[m] = atoi(token.c_str());
		getline(s, token, ',');
		probabilityIn[m] = atof(token.c_str());
		getline(s, token, ',');
		orderQuantities[m] = atof(token.c_str());
		getline(input_file, x, '\n');
		s.clear();
		s.str(x);
	}
	getline(s, token, ',');
	getline(s, token, '\n');
	KR_In= atof(token.c_str());

	ProblemInstance* p = new ProblemInstance(insNumberIn, numSuppliersIn);
	p->AssignCostArray(costIn);
	p->AssignProbArray(probabilityIn);
	p->AssignExecutionCost(executiveCostIn);

	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine demandGenerator(seed1);
	std::normal_distribution<double> distribution1(demandMean, demandSD);
	std::default_random_engine priceGenerator(seed2);
	std::normal_distribution<double> distribution2(spotMarketMean, spotMarketSD);
		
	stringstream convert_t, convert_ins;
	convert_t << numSuppliersIn;
	convert_ins << insNumberIn;
	
	string out_file = FolderNameIn + "/" + "UB_Scenarios_" + convert_ins.str() + ".csv";
	ofstream output;
	output.open(out_file.c_str());
	output << " ,Demand,SP,";
	for (int k = 0; k < numSuppliers; k++) {
		output << "," << k + 1;
	}
	output << ",Cost" << endl;
	for (int i = 0; i < numScenariosForUB; i++)
	{
		double ranno_demand = round(distribution1(demandGenerator)*1000.0) / 1000.0;
		double ranno_price = round(distribution2(priceGenerator)*1000.0) / 1000.0;
		output << i + 1 << "," << ranno_demand << "," << ranno_price << "," << i + 1;
		double sumOrderingCost = 0.0;
		Qtotal = 0;
		for (int k = 0; k < numSuppliers; k++) {
			int availableInScenario = 0;
			double ranno = round(((double)rand() / RAND_MAX)*1000.0) / 1000.0;
			if (ranno < p->prob_tmp[k]) {
				availableInScenario = 1;
				sumOrderingCost += ((double)costIn[k]) * orderQuantities[k];
				Qtotal += orderQuantities[k];
			}
			output << "," << availableInScenario;
			
		}
		double partialCost = CalcScenarioCostOfGivenSolution(p, ranno_demand, ranno_price, Qtotal, KR_In);
		double scenarioCost = sumOrderingCost + partialCost;
		output << "," << scenarioCost;
		output << endl;
	}
	output.close();
	p->~ProblemInstance();


	
}


// First set of tests , includes calculating KR 

void DataGenerator::WriteExactModelsInGamsWithKRIPOPT(string FolderName)
{
	vector<ProblemInstance*>::iterator it;
	DataGenerator::GenerateAllPossibleScenarios();
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		//for (int j = 0; j < numInstances; j++) {
			stringstream convert_t, convert_ins, convert_set;
			convert_t << numSuppliers;
			convert_ins << (*it)->instanceNo;
			//convert_set << j;
			string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_Exact_WithKR.gms";
			ofstream output;
			output.open(out_file.c_str());
			output << "Sets" << endl;
			output << "i suppliers /1*" << numSuppliers << "/" << endl;
			output << "u scenarios /1*" << possbileScenarios << "/" << endl;
			output << endl;
			output << "Parameters" << endl;
			output << "c(i) selling price" << endl;
			output << "/" << endl;
			for (int i = 0; i < numSuppliers; i++)
			{
				output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
			}
			output << "/" << endl;
			output << "pi(i) probability of each supplier's availability" << endl;
			output << "/" << endl;
			for (int i = 0; i < numSuppliers; i++)
			{
				output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
			}
			output << "/" << endl;
			output << "Prob(u) Scenario probability" << endl;
			output << "/" << endl;
			for (int i = 0; i < possbileScenarios; i++) {
				output << i + 1 << '\t' << instanceProbabilities[i] << endl;
			}
			output << "/" << endl;
			output << "Table" << endl;
			output << "R(u,i)  availability" << endl;
			output << '\t';
			for (int i = 0; i < numSuppliers; i++) {
				output << i + 1 << '\t';
			}
			output << endl;
			for (int ns = 0; ns < possbileScenarios; ns++) {
				output << ns + 1 << '\t';
				for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
					output << matrixOfScenarios[ns][supplier_tmp] << '\t';
				}
				output << endl;
			}
			output << endl;
			output << "Scalars" << endl;
			output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
			output << "rc reservation cost at the backup supplier /" << rc << " /" << endl;
			output << "ec execution cost at the backup supplier /" << (*it)->ec_tmp << " /" << endl;
			output << "mu demand mean / " << demandMean << " /" << endl;
			output << "sd demand stdev / " << demandSD << " /" << endl;
			output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
			output << "sdSP spot market price stdev / " << spotMarketSD << " /" << endl;
			output << "telapsed;" << endl;
			output << endl;
			output << "Parameter" << endl;
			output << "zSP, pdfSP, cdfSp, delta;" << endl;
			output << "zSP = (ec - muSP) / sdSP;" << endl;
			output << "pdfSP =  exp(-0.5*power(zSp,2))/sqrt(2*3.14156);" << endl;
			output << "cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));" << endl;
			output << "If (zSP >= 0, " << endl;
			output << "         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);" << endl;
			output << "else  " << endl;
			output << "         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);" << endl;
			output << ");" << endl;
			output << endl;
			output << "Variables" << endl;
			output << "z(u) service level for scenario u" << endl;
			output << "zKR(u)" << endl;
			output << "expected expected profit;" << endl;
			output << endl;
			output << "Positive variables" << endl;
			output << "pdf(u)    pdf for scenario u" << endl;
			output << "cdf(u)    cdf for scenario u" << endl;
			output << "TotalQ(u) total order quantity for scenario u" << endl;
			output << "pdfKR(u)  pdf including KR for scenario u" << endl;
			output << "cdfKR(u)  cdf including KR for scenario u" << endl;
			output << "loss(u)   loss value for Qu scenario u" << endl;
			output << "lossKR(u) loss value for Qu + KR for scenario u" << endl;
			output << "cost(u)   cost scenario u" << endl;
			output << "q(i)      order quantity for supplier i" << endl;
			output << "KR        capacity reserved at the backup supplier;" << endl;
			output << endl;
			output << "Equations" << endl;
			output << "objective        define objective function" << endl;
			output << "scenarioCost(u)  scenario cost definition" << endl;
			output << "calcTotalQ(u)    calculate total order for each scenarion" << endl;
			output << "calcLoss(u)      define loss constraint Qu scenario u" << endl;
			output << "calcLossKR(u)    define loss constraint Qu + KR scenario u" << endl;
			output << "calcpdf(u)       calculate pdf of Qu for scenario u" << endl;
			output << "calccdf(u)       calculate cdf of Qu for scenario u" << endl;
			output << "calcpdfKR(u)     calculate pdf of Qu + KR" << endl;
			output << "calccdfKR(u)     calculate cdf of Qu+KR" << endl;
			output << "calcz(u)         calculate z(u)" << endl;
			output << "calczKR(u)       calculate zKR(u);" << endl;
			output << endl;
			output << "objective..              expected=e= sum((u), Prob(u)*cost(u));" << endl;
			output << "scenarioCost(u)..        cost(u)=e=  sum((i),c(i)*q(i)*R(u,i))+ rc*KR +" << endl;
			output << "                         h*(TotalQ(u)-mu) + (h+muSP-delta)* loss(u)+ " << endl;
			output << "                         delta * lossKR(u); " << endl;
			output << "calcTotalQ(u)..          TotalQ(u) =e=  sum((i), q(i)*R(u,i));" << endl;
			output << "calcLoss(u)..            loss(u)=e=sd*(pdf(u) -z(u)*(1-cdf(u)));" << endl;
			output << "calcLossKR(u)..          lossKR(u) =e= sd*(pdfKR(u) -zKR(u)*(1-cdfKR(u)));" << endl;
			output << "calcpdf(u)..             pdf(u)*sqrt(2*3.14156)=e=exp(-0.5*power(z(u),2));" << endl;
			output << "calccdf(u)..             cdf(u)=e= errorf(z(u));" << endl;
			output << "calcpdfKR(u)..           pdfKR(u)*sqrt(2*3.14156)=e=exp(-0.5*power(zKR(u),2));" << endl;
			output << "calccdfKR(u)..           cdfKR(u)=e= errorf(zKR(u));" << endl;
			output << "calcz(u)..               TotalQ(u)=e=mu+sd*z(u); " << endl;
			output << "calczKR(u)..             TotalQ(u)+KR=e=mu+sd*zKR(u);" << endl;
			output << "cdf.up(u)=1;" << endl;
			output << "cdf.lo(u)=0;" << endl;
			output << "cdfKR.up(u)=1; " << endl;
			output << "cdfKR.lo(u)=0;" << endl;
			output << "z.up(u)=20; " << endl;
			output << "z.lo(u)=-20;" << endl;
			output << "zKR.up(u)=20;" << endl;
			output << "zKR.lo(u)=-20;" << endl;
			output << endl;
			output << endl;
			output << "Model AONnew / all / ;" << endl;
			output << "option optcr = 1e-5;" << endl;
			output << "option nlp=IPOPT;" << endl;
			output << "option iterlim = 100000000;" << endl;
			output << "option reslim = 5000;" << endl;
			output << "option domlim = 2000000;" << endl;
			output << "solve  AONnew using nlp minimizing expected;" << endl;
			output << endl;
			output << endl;
			output << endl;
			output << "Display AONnew.modelstat;" << endl;
			output << "Display AONnew.solvestat;" << endl;
			output << "telapsed = TimeElapsed;" << endl;
			output << endl;
			output << "file  Results storing the results / results_" + convert_t.str() + "_Exact_IPOPT.txt / ; " << endl;
			output << "Results.ap=1;" << endl;
			output << "put Results;" << endl;
			output << "*put /'Objective,Elapsed time','AONnew.solvestat','AONnew.modelstat/;" << endl;
			output << "put /expected.l','telapsed','AONnew.solvestat','AONnew.modelstat/;" << endl;
			output << "putclose Results;" << endl;
			output.close();
		//}

	}

}

void DataGenerator::WriteExactModelsInGamsWithKRConopt(string FolderName)
{
	vector<ProblemInstance*>::iterator it;
	DataGenerator::GenerateAllPossibleScenarios();
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		//for (int j = 0; j < numInstances; j++) {
		stringstream convert_t, convert_ins, convert_set;
		convert_t << numSuppliers;
		convert_ins << (*it)->instanceNo;
		//convert_set << j;
		string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_Exact_WithKR.gms";
		ofstream output;
		output.open(out_file.c_str());
		output << "Sets" << endl;
		output << "i suppliers /1*" << numSuppliers << "/" << endl;
		output << "u scenarios /1*" << possbileScenarios << "/" << endl;
		output << endl;
		output << "Parameters" << endl;
		output << "c(i) selling price" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "pi(i) probability of each supplier's availability" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "Prob(u) Scenario probability" << endl;
		output << "/" << endl;
		for (int i = 0; i < possbileScenarios; i++) {
			output << i + 1 << '\t' << instanceProbabilities[i] << endl;
		}
		output << "/" << endl;
		output << "Table" << endl;
		output << "R(u,i)  availability" << endl;
		output << '\t';
		for (int i = 0; i < numSuppliers; i++) {
			output << i + 1 << '\t';
		}
		output << endl;
		for (int ns = 0; ns < possbileScenarios; ns++) {
			output << ns + 1 << '\t';
			for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
				output << matrixOfScenarios[ns][supplier_tmp] << '\t';
			}
			output << endl;
		}
		output << endl;
		output << "Scalars" << endl;
		output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
		output << "rc reservation cost at the backup supplier /" << rc << " /" << endl;
		output << "ec execution cost at the backup supplier /" << (*it)->ec_tmp << " /" << endl;
		output << "mu demand mean / " << demandMean << " /" << endl;
		output << "sd demand stdev / " << demandSD << " /" << endl;
		output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
		output << "sdSP spot market price stdev / " << spotMarketSD << " /" << endl;
		output << "telapsed;" << endl;
		output << endl;
		output << "Parameter" << endl;
		output << "zSP, pdfSP, cdfSp, delta;" << endl;
		output << "zSP = (ec - muSP) / sdSP;" << endl;
		output << "pdfSP =  exp(-0.5*power(zSp,2))/sqrt(2*3.14156);" << endl;
		output << "cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));" << endl;
		output << "If (zSP >= 0, " << endl;
		output << "         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);" << endl;
		output << "else  " << endl;
		output << "         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);" << endl;
		output << ");" << endl;
		output << endl;
		output << "Variables" << endl;
		output << "z(u) service level for scenario u" << endl;
		output << "zKR(u)" << endl;
		output << "expected expected profit;" << endl;
		output << endl;
		output << "Positive variables" << endl;
		output << "pdf(u)    pdf for scenario u" << endl;
		output << "cdf(u)    cdf for scenario u" << endl;
		output << "TotalQ(u) total order quantity for scenario u" << endl;
		output << "pdfKR(u)  pdf including KR for scenario u" << endl;
		output << "cdfKR(u)  cdf including KR for scenario u" << endl;
		output << "loss(u)   loss value for Qu scenario u" << endl;
		output << "lossKR(u) loss value for Qu + KR for scenario u" << endl;
		output << "cost(u)   cost scenario u" << endl;
		output << "q(i)      order quantity for supplier i" << endl;
		output << "KR        capacity reserved at the backup supplier;" << endl;
		output << endl;
		output << "Equations" << endl;
		output << "objective        define objective function" << endl;
		output << "scenarioCost(u)  scenario cost definition" << endl;
		output << "calcTotalQ(u)    calculate total order for each scenarion" << endl;
		output << "calcLoss(u)      define loss constraint Qu scenario u" << endl;
		output << "calcLossKR(u)    define loss constraint Qu + KR scenario u" << endl;
		output << "calcpdf(u)       calculate pdf of Qu for scenario u" << endl;
		output << "calccdf(u)       calculate cdf of Qu for scenario u" << endl;
		output << "calcpdfKR(u)     calculate pdf of Qu + KR" << endl;
		output << "calccdfKR(u)     calculate cdf of Qu+KR" << endl;
		output << "calcz(u)         calculate z(u)" << endl;
		output << "calczKR(u)       calculate zKR(u);" << endl;
		output << endl;
		output << "objective..              expected=e= sum((u), Prob(u)*cost(u));" << endl;
		output << "scenarioCost(u)..        cost(u)=e=  sum((i),c(i)*q(i)*R(u,i))+ rc*KR +" << endl;
		output << "                         h*(TotalQ(u)-mu) + (h+muSP-delta)* loss(u)+ " << endl;
		output << "                         delta * lossKR(u); " << endl;
		output << "calcTotalQ(u)..          TotalQ(u) =e=  sum((i), q(i)*R(u,i));" << endl;
		output << "calcLoss(u)..            loss(u)=e=sd*(pdf(u) -z(u)*(1-cdf(u)));" << endl;
		output << "calcLossKR(u)..          lossKR(u) =e= sd*(pdfKR(u) -zKR(u)*(1-cdfKR(u)));" << endl;
		output << "calcpdf(u)..             pdf(u)*sqrt(2*3.14156)=e=exp(-0.5*power(z(u),2));" << endl;
		output << "calccdf(u)..             cdf(u)=e= errorf(z(u));" << endl;
		output << "calcpdfKR(u)..           pdfKR(u)*sqrt(2*3.14156)=e=exp(-0.5*power(zKR(u),2));" << endl;
		output << "calccdfKR(u)..           cdfKR(u)=e= errorf(zKR(u));" << endl;
		output << "calcz(u)..               TotalQ(u)=e=mu+sd*z(u); " << endl;
		output << "calczKR(u)..             TotalQ(u)+KR=e=mu+sd*zKR(u);" << endl;
		output << "cdf.up(u)=1;" << endl;
		output << "cdf.lo(u)=0;" << endl;
		output << "cdfKR.up(u)=1; " << endl;
		output << "cdfKR.lo(u)=0;" << endl;
		output << "z.up(u)=20; " << endl;
		output << "z.lo(u)=-20;" << endl;
		output << "zKR.up(u)=20;" << endl;
		output << "zKR.lo(u)=-20;" << endl;
		output << endl;
		output << endl;
		output << "Model AONnew / all / ;" << endl;
		output << "option optcr = 1e-5;" << endl;
		output << "option nlp=Conopt;" << endl;
		output << "option iterlim = 100000000;" << endl;
		output << "option reslim = 5000;" << endl;
		output << "option domlim = 2000000;" << endl;
		output << "solve  AONnew using nlp minimizing expected;" << endl;
		output << endl;
		output << endl;
		output << endl;
		output << "Display AONnew.modelstat;" << endl;
		output << "Display AONnew.solvestat;" << endl;
		output << "telapsed = TimeElapsed;" << endl;
		output << endl;
		output << "file  Results storing the results / results_" + convert_t.str() + "_Exact_Conopt.txt / ; " << endl;
		output << "Results.ap=1;" << endl;
		output << "put Results;" << endl;
		output << "*put /'Objective,Elapsed time','AONnew.solvestat','AONnew.modelstat/;" << endl;
		output << "put /expected.l','telapsed','AONnew.solvestat','AONnew.modelstat/;" << endl;
		output << "putclose Results;" << endl;
		output.close();
		//}

	}

}

void DataGenerator::WriteExactModelsInGamsWithKRLindoglobal(string FolderName)
{
	vector<ProblemInstance*>::iterator it;
	DataGenerator::GenerateAllPossibleScenarios();
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		//for (int j = 0; j < numInstances; j++) {
		stringstream convert_t, convert_ins, convert_set;
		convert_t << numSuppliers;
		convert_ins << (*it)->instanceNo;
		//convert_set << j;
		string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_Exact_WithKR.gms";
		ofstream output;
		output.open(out_file.c_str());
		output << "Sets" << endl;
		output << "i suppliers /1*" << numSuppliers << "/" << endl;
		output << "u scenarios /1*" << possbileScenarios << "/" << endl;
		output << endl;
		output << "Parameters" << endl;
		output << "c(i) selling price" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "pi(i) probability of each supplier's availability" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "Prob(u) Scenario probability" << endl;
		output << "/" << endl;
		for (int i = 0; i < possbileScenarios; i++) {
			output << i + 1 << '\t' << instanceProbabilities[i] << endl;
		}
		output << "/" << endl;
		output << "Table" << endl;
		output << "R(u,i)  availability" << endl;
		output << '\t';
		for (int i = 0; i < numSuppliers; i++) {
			output << i + 1 << '\t';
		}
		output << endl;
		for (int ns = 0; ns < possbileScenarios; ns++) {
			output << ns + 1 << '\t';
			for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
				output << matrixOfScenarios[ns][supplier_tmp] << '\t';
			}
			output << endl;
		}
		output << endl;
		output << "Scalars" << endl;
		output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
		output << "rc reservation cost at the backup supplier /" << rc << " /" << endl;
		output << "ec execution cost at the backup supplier /" << (*it)->ec_tmp << " /" << endl;
		output << "mu demand mean / " << demandMean << " /" << endl;
		output << "sd demand stdev / " << demandSD << " /" << endl;
		output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
		output << "sdSP spot market price stdev / " << spotMarketSD << " /" << endl;
		output << "telapsed;" << endl;
		output << endl;
		output << "Parameter" << endl;
		output << "zSP, pdfSP, cdfSp, delta;" << endl;
		output << "zSP = (ec - muSP) / sdSP;" << endl;
		output << "pdfSP =  exp(-0.5*power(zSp,2))/sqrt(2*3.14156);" << endl;
		output << "cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));" << endl;
		output << "If (zSP >= 0, " << endl;
		output << "         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);" << endl;
		output << "else  " << endl;
		output << "         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);" << endl;
		output << ");" << endl;
		output << endl;
		output << "Variables" << endl;
		output << "z(u) service level for scenario u" << endl;
		output << "zKR(u)" << endl;
		output << "expected expected profit;" << endl;
		output << endl;
		output << "Positive variables" << endl;
		output << "pdf(u)    pdf for scenario u" << endl;
		output << "cdf(u)    cdf for scenario u" << endl;
		output << "TotalQ(u) total order quantity for scenario u" << endl;
		output << "pdfKR(u)  pdf including KR for scenario u" << endl;
		output << "cdfKR(u)  cdf including KR for scenario u" << endl;
		output << "loss(u)   loss value for Qu scenario u" << endl;
		output << "lossKR(u) loss value for Qu + KR for scenario u" << endl;
		output << "cost(u)   cost scenario u" << endl;
		output << "q(i)      order quantity for supplier i" << endl;
		output << "KR        capacity reserved at the backup supplier;" << endl;
		output << endl;
		output << "Equations" << endl;
		output << "objective        define objective function" << endl;
		output << "scenarioCost(u)  scenario cost definition" << endl;
		output << "calcTotalQ(u)    calculate total order for each scenarion" << endl;
		output << "calcLoss(u)      define loss constraint Qu scenario u" << endl;
		output << "calcLossKR(u)    define loss constraint Qu + KR scenario u" << endl;
		output << "calcpdf(u)       calculate pdf of Qu for scenario u" << endl;
		output << "calccdf(u)       calculate cdf of Qu for scenario u" << endl;
		output << "calcpdfKR(u)     calculate pdf of Qu + KR" << endl;
		output << "calccdfKR(u)     calculate cdf of Qu+KR" << endl;
		output << "calcz(u)         calculate z(u)" << endl;
		output << "calczKR(u)       calculate zKR(u);" << endl;
		output << endl;
		output << "objective..              expected=e= sum((u), Prob(u)*cost(u));" << endl;
		output << "scenarioCost(u)..        cost(u)=e=  sum((i),c(i)*q(i)*R(u,i))+ rc*KR +" << endl;
		output << "                         h*(TotalQ(u)-mu) + (h+muSP-delta)* loss(u)+ " << endl;
		output << "                         delta * lossKR(u); " << endl;
		output << "calcTotalQ(u)..          TotalQ(u) =e=  sum((i), q(i)*R(u,i));" << endl;
		output << "calcLoss(u)..            loss(u)=e=sd*(pdf(u) -z(u)*(1-cdf(u)));" << endl;
		output << "calcLossKR(u)..          lossKR(u) =e= sd*(pdfKR(u) -zKR(u)*(1-cdfKR(u)));" << endl;
		output << "calcpdf(u)..             pdf(u)*sqrt(2*3.14156)=e=exp(-0.5*power(z(u),2));" << endl;
		output << "calccdf(u)..             cdf(u)=e= errorf(z(u));" << endl;
		output << "calcpdfKR(u)..           pdfKR(u)*sqrt(2*3.14156)=e=exp(-0.5*power(zKR(u),2));" << endl;
		output << "calccdfKR(u)..           cdfKR(u)=e= errorf(zKR(u));" << endl;
		output << "calcz(u)..               TotalQ(u)=e=mu+sd*z(u); " << endl;
		output << "calczKR(u)..             TotalQ(u)+KR=e=mu+sd*zKR(u);" << endl;
		output << "cdf.up(u)=1;" << endl;
		output << "cdf.lo(u)=0;" << endl;
		output << "cdfKR.up(u)=1; " << endl;
		output << "cdfKR.lo(u)=0;" << endl;
		output << "z.up(u)=20; " << endl;
		output << "z.lo(u)=-20;" << endl;
		output << "zKR.up(u)=20;" << endl;
		output << "zKR.lo(u)=-20;" << endl;
		output << endl;
		output << endl;
		output << "Model AONnew / all / ;" << endl;
		output << "option optcr = 1e-5;" << endl;
		output << "option nlp=lindoglobal;" << endl;
		output << "option iterlim = 100000000;" << endl;
		output << "option reslim = 5000;" << endl;
		output << "option domlim = 2000000;" << endl;
		output << "solve  AONnew using nlp minimizing expected;" << endl;
		output << endl;
		output << endl;
		output << endl;
		output << "Display AONnew.modelstat;" << endl;
		output << "Display AONnew.solvestat;" << endl;
		output << "telapsed = TimeElapsed;" << endl;
		output << endl;
		output << "file  Results storing the results / results_" + convert_t.str() + "_Exact_LG.txt / ; " << endl;
		output << "Results.ap=1;" << endl;
		output << "put Results;" << endl;
		output << "*put /'Objective,Elapsed time','AONnew.solvestat','AONnew.modelstat/;" << endl;
		output << "put /expected.l','telapsed','AONnew.solvestat','AONnew.modelstat/;" << endl;
		output << "putclose Results;" << endl;
		output.close();
		//}

	}

}


void DataGenerator::WriteSAAModelsInGams(string FolderName) {
	vector<ProblemInstance*>::iterator it;
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		for (int j = 0; j < numSets; j++) {
			stringstream convert_t, convert_ins, convert_set;
			convert_t << numSuppliers;
			convert_ins << (*it)->instanceNo;
			convert_set << j;
			string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() +"_SAA_Set_" + convert_set.str() + ".gms";
			ofstream output;
			output.open(out_file.c_str());
			output << "Scalar" << endl;
			output << "tcomp, texec, telapsed;" << endl;
			output << "Sets" << endl;
			output << "i suppliers /1*" << numSuppliers << "/" << endl;
			output << "m SAA scenarios /1*" << numScenarios << "/" << endl;
			output << "j   coefficients  /Demand , SP/" << endl;
			output << "u scenarios /1*" << possbileScenarios << "/" << endl;
			output << endl;
			output << endl;
			output << "Scalars" << endl;
			output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
			output << "rc reservation cost at the backup supplier /" << rc << " /" << endl;
			output << "ec execution cost at the backup supplier /" << (*it)->ec_tmp << " /" << endl;
			output << "mu demand mean / " << demandMean << " /" << endl;
			output << "sd demand stdev / " << demandSD << " /" << endl;
			output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
			output << "sdSP spot market price stdev / " << spotMarketSD << " /;" << endl;
			output << endl;
			output << "Parameters" << endl;
			output << "oCost(i) ordering cost" << endl;
			output << "/" << endl;
			for (int i = 0; i < numSuppliers; i++)
			{
				int oCost_temp = 0; // floor(0.25 * ((*it)->cost_tmp[i]));
				output << i + 1 << '\t' << oCost_temp << endl;
			}
			output << "/" << endl;
			output << "c(i) selling price" << endl;
			output << "/" << endl;
			for (int i = 0; i < numSuppliers; i++)
			{
				output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
			}
			output << "/" << endl;
			output << "pi(i) probability of each supplier's availability" << endl;
			output << "/" << endl;
			for (int i = 0; i < numSuppliers; i++)
			{
				output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
			}
			output << "/" << endl;
			output << "Prob(u) Scenario probability" << endl;
			output << "/" << endl;
			for (int i = 0; i < possbileScenarios; i++) {
				output << i + 1 << '\t' << instanceProbabilities[i] << endl;
			}
			output << "/" << endl;
			output << "Table" << endl;
			output << "RPrime(u,i)  availability" << endl;
			output << '\t';
			for (int i = 0; i < numSuppliers; i++) {
				output << i + 1 << '\t';
			}
			output << endl;
			for (int ns = 0; ns < possbileScenarios; ns++) {
				output << ns + 1 << '\t';
				for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
					output << matrixOfScenarios[ns][supplier_tmp] << '\t';
				}
				output << endl;
			}
			output << "Parameters " << endl;
			output << "data(m,j) demand and SP values in each scenario " << endl;
			output << "R(m,i) availability of the suppliers in each scenario" << endl;
			output << "sumOrders, sumSelectedSuppliers; " << endl;
			output << endl;
			output << "Variables" << endl;
			output << "expected expected profit;" << endl;
			output << endl;
			output << "Positive variables" << endl;
			output << "TotalQ(m) total order quantity for scenario m " << endl;
			output << "Ipos(m)   amount of overage in scenarion m " << endl;
			output << "Ineg(m)   amount of underage in scenario m " << endl;
			output << "O(m)      amount of the option contract in scenario m" << endl;
			output << "S(m)      spot market purchase in scenario m" << endl;
			output << "cost(m)  cost scenario m" << endl;
			output << "q(i) order quantity for supplier i" << endl;
			output << "KR   capacity reserved at the backup supplier;" << endl;
			output << endl;
			output << "Equations" << endl;
			output << "objective        define objective function" << endl;
			output << "scenarioCost(m)  scenario cost definition " << endl;
			output << "calcTotalQ(m)    calculate total order for each scenarion" << endl;
			output << "overage(m)  " << endl;
			output << "underage(m) " << endl;
			output << "underage2(m)  " << endl;
			output << "optionContract(m); " << endl;
			output << endl;
			output << "objective..         expected=e= (1/" << numScenarios << ") * sum((m),cost(m));" << endl;
			output << "scenarioCost(m)..   cost(m)=e=  sum((i),q(i)*(c(i)*R(m,i)+oCost(i)))+ rc*KR + h*Ipos(m)+ " << endl;
			output << "                         ec* O(m) + data(m, 'SP')*S(m); " << endl;
			output << "calcTotalQ(m)..   TotalQ(m) =e=  sum((i), q(i)*R(m,i)); " << endl;
			output << "overage(m)..      Ipos(m) =g= TotalQ(m) - data(m, 'Demand');" << endl;
			output << "underage(m)..    Ineg(m) =g=  data(m, 'Demand')- TotalQ(m); " << endl;
			output << "underage2(m)..   Ineg(m) =e= O(m)+S(m);  " << endl;
			output << "optionContract(m)..  O(m) =l= KR;" << endl;
			output << endl;
			output << endl;
			output << "Model AONnew / all / ;" << endl;
			output << "$include demand_" << numSuppliers << "_" << (*it)->instanceNo << "_set_" << j << ".inc" << endl;
			output << "$include SP_" << numSuppliers << "_" << (*it)->instanceNo << "_set_" << j << ".inc" << endl;
			output << "$include Rdata_" << numSuppliers << "_" << (*it)->instanceNo << "_set_" << j << ".inc" << endl;
			output << endl;
			output << endl;
			output << "option optcr = 1e-5;" << endl;
			output << "option iterlim = 100000000;" << endl;
			output << "option reslim = 5000;" << endl;
			output << "option domlim = 20000;" << endl;
			output << "solve  AONnew using lp minimizing expected;" << endl;
			output << "Parameter" << endl;
			output << "zSP, pdfSP, cdfSp, delta;" << endl;
			output << "zSP = (ec - muSP) / sdSP;" << endl;
			output << "pdfSP =  exp(-0.5*power(zSp,2))/sqrt(2*3.14156);" << endl;
			output << "cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));" << endl;
			output << "If (zSP >= 0, " << endl;
			output << "         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);" << endl;
			output << "else  " << endl;
			output << "         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);" << endl;
			output << ");" << endl;
			output << "Parameter" << endl;
			output << "RealObj" << endl;
			output << "zReal(u) service level for scenario u" << endl;
			output << "zKRReal(u)  " << endl;
			output << "lossReal(u)   loss value for Qu scenario u" << endl;
			output << "lossKRReal(u) loss value for Qu + KR for scenario u" << endl;
			output << "pdfReal(u)    pdf for scenario u" << endl;
			output << "cdfReal(u)    cdf for scenario u " << endl;
			output << "TotalQReal(u) total order quantity for scenario u" << endl;
			output << "pdfKRReal(u)  pdf including KR for scenario u " << endl;
			output << "cdfKRReal(u)  cdf including KR for scenario u" << endl;
			output << "costReal(u)   cost scenario u  ;" << endl;
			output << endl;
			output << "Loop(u," << endl;
			output << "TotalQReal(u)=  sum((i), q.l(i)*RPrime(u,i));" << endl;
			output << "zReal(u)=(TotalQReal(u)-mu)/sd;" << endl;
			output << "zKRReal(u)=(TotalQReal(u)+KR.l-mu)/sd;" << endl;
			output << "pdfReal(u) = exp(-0.5*power(zReal(u),2))/sqrt(2*3.14156);" << endl;
			output << "cdfReal(u) = errorf(zReal(u));" << endl;
			output << "lossReal(u) = sd * (pdfReal(u) - zReal(u)* (1-cdfReal(u)));" << endl;
			output << "pdfKRReal(u) = exp(-0.5*power(zKRReal(u), 2)) / sqrt(2 * 3.14156);" << endl;
			output << "cdfKRReal(u)= errorf(zKRReal(u));" << endl;
			output << endl;
			output << "lossKRReal(u) = sd * (pdfKRReal(u) - zKRReal(u)*(1-cdfKRReal(u)));" << endl;
			output << "costReal(u)=  sum((i),q.l(i)*(c(i)*RPrime(u,i)+oCost(i)))+ rc*KR.l + " << endl;
			output << "                         h*(TotalQReal(u)-mu) + (h+muSP-delta)* lossReal(u)+" << endl;
			output << "                          delta * lossKRReal(u); " << endl;
			output << ");" << endl;
			output << endl;
			output << "RealObj = sum((u), Prob(u)*costReal(u));" << endl;			
			output << endl;
			output << "sumOrders = sum((i), q.l(i));    " << endl;
			output << "sumSelectedSuppliers = 0; " << endl;
			output << "loop(i, " << endl;
			output << "   if(q.l(i)>0,  " << endl;
			output << "      sumSelectedSuppliers= sumSelectedSuppliers+1;" << endl;
			output << "   );  " << endl;
			output << ");" << endl;
			output << endl;
			output << "Display AONnew.modelstat;" << endl;
			output << "Display AONnew.solvestat;" << endl;
			output << "tcomp = TimeComp;" << endl;
			output << "texec = TimeExec;" << endl;
			output << "telapsed = TimeElapsed;" << endl;
			output << "*Display tcomp, texec, telapsed;" << endl;
			output << endl;
			output << "file  Results storing the results / results_" + convert_t.str() + "_SAA.txt / ; " << endl;
			//output << "file  Results storing the results / results_" + convert_t.str() + "_" + convert_ins.str() + "_SAA.txt / ; " << endl;
			output << "Results.ap=1;" << endl;
			output << "put Results;" << endl;
			output << "put /expected.l','telapsed','sumSelectedSuppliers','sumOrders','KR.l','RealObj/;" << endl;
			//output << "put /expected.l','telapsed','RealObj/;" << endl;
			output.close();

		}

	}
}

void DataGenerator::WriteExactGAMSBatchFile(string FolderName)
{
	stringstream convert_temp;
	convert_temp << numSuppliers;
	string out_file = FolderName + "/Batch_" + convert_temp.str() + "_Exact.gms";
	ofstream output;
	output.open(out_file.c_str());
	//vector<ProblemInstance*>::iterator it;
	for (int j = 0; j <numInstances; j++)
	{			
			output << "$call gams " <<numSuppliers << "_" << j << "_Exact_WithKR" << endl;
	}
	output.close();	
}

void DataGenerator::WriteSaaGAMSBatchFile(string FolderName) 
{
	stringstream convert_temp;
	convert_temp << numSuppliers;
	for (int j = 0; j < numInstances; j++) {
		stringstream convert_ins;
		convert_ins << j;
		string out_file = FolderName + "/Batch_" + convert_temp.str() +"_" +convert_ins.str()+ "_SAA.gms";
		ofstream output;
		output.open(out_file.c_str());
		for (int k = 0; k < numSets; k++) {
			output << "$call gams " << numSuppliers << "_" << j << "_SAA_Set_" << k << endl;
		}
		output.close();
	}
}

void DataGenerator::WriteSaaGAMSBatchBatchFile(string FolderName)
{
	stringstream convert_temp;
	convert_temp << numSuppliers;
	string out_file = FolderName + "/Total_Batch_" + convert_temp.str() + "_SAA.gms";
	ofstream output;
	output.open(out_file.c_str());
	//vector<ProblemInstance*>::iterator it;
	for (int j = 0; j <numInstances; j++)
	{
		output << "$call gams " << "Batch_" << convert_temp.str() << "_" << j << "_SAA" << endl;
	}
	output.close();
}

void DataGenerator::WriteHeuristicModelsGamsWithKR(string FolderName) {

	vector<ProblemInstance*>::iterator it;
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		stringstream convert_t, convert_ins, convert_set;
		convert_t << numSuppliers;
		convert_ins << (*it)->instanceNo;
		string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_Heuristic.gms";
		ofstream output;
		output.open(out_file.c_str());
		output << "$FuncLibIn stolib stodclib  " << endl;
		output << "function icdfNormalTEMP     /stolib.icdfNormal/" << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "Scalar" << endl;
		output << "tcomp, texec, telapsed;" << endl;
		output << "*------------------------------------------------------------------------------ -" << endl;
		output << "Sets" << endl;
		output << "i suppliers / 1 *" << numSuppliers << " /" << endl;
		output << "iter iterator / 1 *" << numSuppliers + 1 << " /" << endl;
		int combinations = pow(2, numSuppliers);
		output << "u scenarios / 1 *" << combinations << " /" << endl;
		output << "cycle(iter);" << endl;
		output << "alias (iter, k);" << endl;
		output << "alias(iter,l);" << endl;
		output << "parameter A(iter,k);" << endl;
		output << endl;
		output << "*------------------------------------------------------------------------------ -" << endl;
		output << "Parameters" << endl;
		output << "c(i) selling price" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "pi(i) probability of each supplier's availability" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "cr(i) critical ratio " << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			double critical_frac = (spotMarketMean - ((double)((*it)->cost_tmp[i]))) / (spotMarketMean + ((double)((*it)->hCostInstance)));
			output << i + 1 << '\t' << critical_frac << endl;
		}
		output << "/" << endl;
		output << "cdfINV_temp(i) inverse of F for critical values" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			double critical_frac = (spotMarketMean - ((double)((*it)->cost_tmp[i]))) / (spotMarketMean + ((double)((*it)->hCostInstance)));
			double cdfinv = DataGenerator::Norm_CDF_INV(critical_frac, demandMean, demandSD);
			output << i + 1 << '\t' << cdfinv << endl;
		}
		output << "/" << endl;
		output << "index position of the backup supplier according to the critical value" << endl;
		output << "fractileKR, cdfINVKR, piPrime(k),crPrime(k), cdfINVPrime(k), cdfINV(iter); " << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Table" << endl;
		output << "R(u,i)  availability" << endl;
		output << '\t';
		for (int i = 0; i < numSuppliers; i++) {
			output << i + 1 << '\t';
		}
		output << endl;
		for (int ns = 0; ns < possbileScenarios; ns++) {
			output << ns + 1 << '\t';
			for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
				output << matrixOfScenarios[ns][supplier_tmp] << '\t';
			}
			output << endl;
		}
		output << "parameter" << endl;
		output << "Prob(u) Scenario probability" << endl;
		output << "/" << endl;
		for (int i = 0; i < possbileScenarios; i++) {
			output << i + 1 << '\t' << instanceProbabilities[i] << endl;
		}
		output << "/;" << endl;
		output << endl;
		output << "Scalar" << endl;
		output << "tcomp, texec, telapsed; " << endl;
		output << "Scalars" << endl;
		output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
		output << "rc reservation cost at the backup supplier /" << rc << " /" << endl;
		output << "ec execution cost at the backup supplier /" << (*it)->ec_tmp << " /" << endl;
		output << "mu demand mean / " << demandMean << " /" << endl;
		output << "sd demand stdev / " << demandSD << " /" << endl;
		output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
		output << "sdSP spot market price stdev / " << spotMarketSD << " /;" << endl;
		output << endl;
		output << "Parameter" << endl;
		output << "zSP, pdfSP, cdfSp, delta, counter, check, testValKR;" << endl;
		output << "zSP = (ec - muSP) / sdSP;" << endl;
		output << "pdfSP =  exp(-0.5*power(zSp,2))/sqrt(2*3.14156);" << endl;
		output << "cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));" << endl;
		output << "If (zSP >= 0, " << endl;
		output << "         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);" << endl;
		output << "else  " << endl;
		output << "         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);" << endl;
		output << ");" << endl;
		output << "*----------- Find the index (position) of KR among suppliers" << endl;
		output << "fractileKR= (delta-rc) / delta;" << endl;
		output << "cdfINVKR =icdfNormalTEMP(fractileKR,mu,sd);" << endl;
		output << "index=1;" << endl;
		output << "loop(i," << endl;
		output << "if (fractileKR < cr(i)," << endl;
		output << "        index=index+1;" << endl;
		output << "));" << endl;
		output << "*----------- Update pi, cr and cdfINV vectors by inserting KR in its position" << endl;
		output << "loop(k," << endl;
		output << "loop(i," << endl;
		output << "if ((ord(k) < index) and (ord(k)  eq ord(i))," << endl;
		output << "         piPrime(k) = pi(i);" << endl;
		output << "         crPrime(k)= cr(i);" << endl;
		output << "         cdfINVPrime(k)= cdfINV_temp(i); " << endl;
		output << "elseif (ord(k) = index)," << endl;
		output << "         piPrime(k)=delta/(h+muSP);" << endl;
		output << "         cdfINVPrime(k)= cdfINVKR;" << endl;
		output << "         crPrime(k)= fractileKR;" << endl;
		output << "         cdfINVPrime(k)= cdfINVKR;" << endl;
		output << "elseif ((ord(k) > index) and (ord(k) eq ord(i)+1))," << endl;
		output << "         piPrime(k) = pi(i);" << endl;
		output << "         crPrime(k)= cr(i);" << endl;
		output << "         cdfINVPrime(k)= cdfINV_temp(i); " << endl;
		output << ");" << endl;
		output << ");" << endl;
		output << ");" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << endl;
		output << "Variables" << endl;
		output << "z service level for scenario u" << endl;
		output << "zKR" << endl;
		output << "expected expected profit;" << endl;
		output << endl;
		output << "Positive variables" << endl;
		output << "pdf    pdf for scenario u" << endl;
		output << "cdf    cdf for scenario u" << endl;
		output << "TotalQ total order quantity for scenario u" << endl;
		output << "pdfKR  pdf including KR for scenario u" << endl;
		output << "cdfKR  cdf including KR for scenario u" << endl;
		output << "loss   loss value for Qu scenario u" << endl;
		output << "lossKR loss value for Qu + KR for scenario u" << endl;
		output << "q(i)      order quantity for supplier i" << endl;
		output << "qPrime(k) q for the converted problem" << endl;
		output << "KR        capacity reserved at the backup supplier;" << endl;
		output << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Equations" << endl;
		output << "objective        define objective function" << endl;
		output << "calcTotalQ    calculate total order for each scenarion" << endl;
		output << "calcLoss      define loss constraint Q" << endl;
		output << "calcLossKR    define loss constraint Q + KR" << endl;
		output << "calcz         calculate z" << endl;
		output << "calczKR       calculate zKR" << endl;
		output << "calcpdf" << endl;
		output << "calccdf" << endl;
		output << "calcpdfKR" << endl;
		output << "calccdfKR" << endl;
		output << "PrimeToReal1" << endl;
		output << "PrimeToReal2" << endl;
		output << "PrimeToReal3" << endl;
		output << "LinConst; " << endl;
		output << endl;
		output << "objective..              expected=e= sum((i),(c(i)+h)*q(i)*pi(i))+ rc*KR -" << endl;
		output << "                         h*mu + (h+muSP-delta)* loss+ delta * lossKR;" << endl;
		output << "calcTotalQ..             TotalQ =e=  sum((i), pi(i)*q(i));" << endl;
		output << "calcLoss..               loss =e=sd*(pdf -z*(1-cdf));" << endl;
		output << "calcLossKR..             lossKR =e= sd*(pdfKR -zKR*(1-cdfKR));" << endl;
		output << "calcz..                  TotalQ=e=mu+sd*z;" << endl;
		output << "calczKR..                TotalQ+KR=e=mu+sd*zKR;" << endl;
		output << "calcpdf..                pdf*sqrt(2*3.14156)=e=exp(-0.5*power(z,2));" << endl;
		output << "calccdf..                cdf=e= errorf(z);" << endl;
		output << "calcpdfKR..              pdfKR*sqrt(2*3.14156)=e=exp(-0.5*power(zKR,2));" << endl;
		output << "calccdfKR..              cdfKR=e= errorf(zKR);" << endl;
		output << "PrimeToReal1(i,k)$((ord(k) < index) and (ord(i) eq ord(k)))..    q(i)=e=qPrime(k);" << endl;
		output << "PrimeToReal2(k)$(ord(k) eq index)..    KR=e=qPrime(k); " << endl;
		output << "PrimeToReal3(i,k)$((ord(k) > index) and (ord(i)+1 eq ord(k)))..    q(i)=e=qPrime(k);" << endl;
		output << "LinConst(cycle)..        sum((k)$(ord(k) le counter),(A(cycle,k)*qPrime(k))) =e= cdfINV(cycle);" << endl;
		//output << "LinConst(cycle)..        sum((i),(A(i,cycle)*q(i))) =e= cdfINV(cycle);" << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "Model AONnew /All/;" << endl;
		output << "option limrow=100;" << endl;
		output << "cycle(iter)=no;" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "******** Updating the left - hand - side : matrix A" << endl;
		output << "loop(iter," << endl;
		output << "cycle(iter-1)$(ord(iter)<>1)=no;" << endl;
		output << "cycle(iter)=yes;" << endl;
		output << "A(cycle,k)$(ord(iter) <> ord(k))=piPrime(k);" << endl;
		output << "A(cycle,k)$(ord(iter) eq ord(k))=1;" << endl;
		output << "display A;" << endl;
		output << ");" << endl;
		output << "cycle(iter) = no;" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "******** Updating the right-hand-side : parameter " << endl;
		output << "loop(iter," << endl;
		output << "cycle(iter - 1)$(ord(iter)<>1) = no;" << endl;
		output << "cycle(iter) = yes;" << endl;
		output << "   loop(k," << endl;
		output << "      if (ord(iter)=ord(k)," << endl;
		output << "         cdfINV(cycle)=cdfINVPrime(k);" << endl;
		output << "     )); " << endl;
		output << "display cdfINVPrime;" << endl;
		output << ");" << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "********* Solving the problem " << endl;
		output << "parameters " << endl;
		output << endl;
		output << "RealObj" << endl;
		output << "RealObj_temp" << endl;
		output << "zReal(u) service level for scenario u" << endl;
		output << "zKRReal(u) " << endl;
		output << "pdfReal(u)    pdf for scenario u  " << endl;
		output << "cdfReal(u)    cdf for scenario u" << endl;
		output << "TotalQReal(u) total order quantity for scenario u" << endl;
		output << "pdfKRReal(u)  pdf including KR for scenario u " << endl;
		output << "cdfKRReal(u)  cdf including KR for scenario u" << endl;
		output << "lossReal(u)   loss value for Qu scenario u" << endl;
		output << "lossKRReal(u) loss value for Qu + KR for scenario u" << endl;
		output << "costReal(u)   cost scenario u;" << endl;
		output << endl;
		output << "*---- Solve the approximate problem with KR in the linear equation set " << endl;
		output << "RealObj=100000000; " << endl;
		output << "counter=0;" << endl;
		output << "cycle(iter)=no;" << endl;
		output << "loop(iter," << endl;
		output << "		counter=counter+1;" << endl;
		output << "		cycle(iter)=yes; " << endl;
		output << "		qPrime.lo(k)$(ord(k) le ord(iter))=0;" << endl;
		output << "		qPrime.up(k)$(ord(k) le ord(iter))=+inf; " << endl;
		output << "		qPrime.fx(k)$(ord(k) > ord(iter))=0;" << endl;
		output << "		option nlp=Conopt;" << endl;
		output << "		solve  AONnew using nlp minimizing expected;" << endl;
		output << "		display qPrime.l,cdfINV;" << endl;
		output << "TotalQReal(u)=  sum((i), q.l(i)*R(u,i));" << endl;
		output << "zReal(u) = (TotalQReal(u) - mu) / sd;" << endl;
		output << "zKRReal(u) = (TotalQReal(u) + KR.l - mu) / sd;" << endl;
		output << "pdfReal(u) = exp(-0.5*power(zReal(u), 2)) / sqrt(2 * 3.14156);" << endl;
		output << "cdfReal(u) = errorf(zReal(u)); " << endl;
		output << "Loop(u," << endl;
		output << "		 lossReal(u) = sd * (pdfReal(u) - zReal(u)*(1 - cdfReal(u)));" << endl;
		output << ");" << endl;
		output << "pdfKRReal(u) = exp(-0.5*power(zKRReal(u), 2)) / sqrt(2 * 3.14156);" << endl;
		output << "cdfKRReal(u)= errorf(zKRReal(u));" << endl;
		output << "Loop(u," << endl;
		output << "		lossKRReal(u) = sd * (pdfKRReal(u) - zReal(u)*(1 - cdfKRReal(u)));" << endl;
		output << ");" << endl;
		output << endl;
		output << "costReal(u)=  sum((i),c(i)*q.l(i)*R(u,i))+ rc*KR.l  +" << endl;
		output << "                         h*(TotalQReal(u)-mu) + (h+muSP-delta)* lossReal(u)+" << endl;
		output << "                         delta * lossKRReal(u); " << endl;
		output << "RealObj_temp= sum((u), Prob(u)*costReal(u)); " << endl;
		output << "Display RealObj; " << endl;
		output << "If (RealObj_temp < RealObj," << endl;
		output << "         RealObj=RealObj_temp;" << endl;
		output << ");  " << endl;
		output << ");  " << endl;
		output << "*----- Solve the approximate problem with KR set to a value =0 " << endl;
		output << "testValKR=0;" << endl;
		output << "counter=0;" << endl;
		output << "cycle(iter)=no;" << endl;
		output << "loop(iter," << endl;
		output << "		counter=counter+1;" << endl;
		output << "		cycle(iter)=yes; " << endl;
		output << "		qPrime.lo(k)$(ord(k) le ord(iter))=0;" << endl;
		output << "		qPrime.up(k)$(ord(k) le ord(iter))=+inf; " << endl;
		output << "		qPrime.fx(k)$(ord(k) > ord(iter))=0;" << endl;
		output << "     KR.fx=testValKR; " << endl;
		output << "		option nlp=Conopt;" << endl;
		output << "		solve  AONnew using nlp minimizing expected;" << endl;
		output << "		display qPrime.l,cdfINV;" << endl;
		output << "TotalQReal(u)=  sum((i), q.l(i)*R(u,i));" << endl;
		output << "zReal(u) = (TotalQReal(u) - mu) / sd;" << endl;
		output << "zKRReal(u) = (TotalQReal(u) + KR.l - mu) / sd;" << endl;
		output << "pdfReal(u) = exp(-0.5*power(zReal(u), 2)) / sqrt(2 * 3.14156);" << endl;
		output << "cdfReal(u) = errorf(zReal(u)); " << endl;
		output << "Loop(u," << endl;
		output << "		 lossReal(u) = sd * (pdfReal(u) - zReal(u)*(1 - cdfReal(u)));" << endl;
		output << ");" << endl;
		output << "pdfKRReal(u) = exp(-0.5*power(zKRReal(u), 2)) / sqrt(2 * 3.14156);" << endl;
		output << "cdfKRReal(u)= errorf(zKRReal(u));" << endl;
		output << "Loop(u," << endl;
		output << "		lossKRReal(u) = sd * (pdfKRReal(u) - zReal(u)*(1 - cdfKRReal(u)));" << endl;
		output << ");" << endl;
		output << endl;
		output << "costReal(u)=  sum((i),c(i)*q.l(i)*R(u,i))+ rc*KR.l  +" << endl;
		output << "                         h*(TotalQReal(u)-mu) + (h+muSP-delta)* lossReal(u)+" << endl;
		output << "                         delta * lossKRReal(u); " << endl;
		output << "RealObj_temp= sum((u), Prob(u)*costReal(u)); " << endl;
		output << "Display RealObj; " << endl;
		output << "If (RealObj_temp < RealObj," << endl;
		output << "         RealObj=RealObj_temp;" << endl;
		output << ");  " << endl;
		output << ");  " << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Display RealObj;" << endl;
		output << "tcomp = TimeComp;" << endl;
		output << "texec = TimeExec;" << endl;
		output << "telapsed = TimeElapsed;" << endl;
		//output << "*Display tcomp, texec, telapsed;" << endl;
		output << endl;
		output << "file  Results storing the results / results_" + convert_t.str() + "_" + convert_ins.str() + "_Heuristic.txt / ; " << endl;
		output << "Results.ap=1;" << endl;
		output << "put Results;" << endl;
		output << "put /Realobj','telapsed/;" << endl;
		output << "putclose Results;" << endl;
		output.close();
	}
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		(*it)->~ProblemInstance();
	}
}

void DataGenerator::WriteHeuristicModelsGamsWithKR_ApproximateIterative(string FolderName) {

	vector<ProblemInstance*>::iterator it;
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		stringstream convert_t, convert_ins, convert_set;
		convert_t << numSuppliers;
		convert_ins << (*it)->instanceNo;
		string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_HeuristicApproxIterative.gms";
		ofstream output;
		output.open(out_file.c_str());
		output << "$FuncLibIn stolib stodclib  " << endl;
		output << "function icdfNormalTEMP     /stolib.icdfNormal/" << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "Scalar" << endl;
		output << "tcomp, texec, telapsed;" << endl;
		output << "*------------------------------------------------------------------------------ -" << endl;
		output << "Sets" << endl;
		output << "i suppliers / 1 *" << numSuppliers << " /" << endl;
		int combinations = pow(2, numSuppliers);
		output << "u scenarios / 1 *" << combinations << " /" << endl;
		output << "alias(i, j);" << endl;
		output << "alias(i,iter,k);" << endl;
		output << "set cycle(iter);" << endl;
		output << "*------------------------------------------------------------------------------ -" << endl;
		output << "Parameters" << endl;
		output << "c(i) selling price" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "pi(i) probability of each supplier's availability" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "cdfINV_temp(i) inverse of F for critical values" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			double critical_frac = (spotMarketMean - ((double)((*it)->cost_tmp[i]))) / (spotMarketMean + ((double)((*it)->hCostInstance)));
			double cdfinv = DataGenerator::Norm_CDF_INV(critical_frac, demandMean, demandSD);
			output << i + 1 << '\t' << cdfinv << endl;
		}
		output << "/" << endl;
		output << "cdfINV(iter) " << endl;
		output << "partX(i) X part of the equations" << endl;
		output << "cdfINV_partX(i) inverse of F for fractions including X" << endl;
		output << "fractileKR, cdfINVKR, A(iter,i); " << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Table" << endl;
		output << "R(u,i)  availability" << endl;
		output << '\t';
		for (int i = 0; i < numSuppliers; i++) {
			output << i + 1 << '\t';
		}
		output << endl;
		for (int ns = 0; ns < possbileScenarios; ns++) {
			output << ns + 1 << '\t';
			for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
				output << matrixOfScenarios[ns][supplier_tmp] << '\t';
			}
			output << endl;
		}
		output << "parameter" << endl;
		output << "Prob(u) Scenario probability" << endl;
		output << "/" << endl;
		for (int i = 0; i < possbileScenarios; i++) {
			output << i + 1 << '\t' << instanceProbabilities[i] << endl;
		}
		output << "/;" << endl;
		output << endl;
		output << "Scalars" << endl;
		output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
		output << "rc reservation cost at the backup supplier /" << rc << " /" << endl;
		output << "ec execution cost at the backup supplier /" << (*it)->ec_tmp << " /" << endl;
		output << "mu demand mean / " << demandMean << " /" << endl;
		output << "sd demand stdev / " << demandSD << " /" << endl;
		output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
		output << "sdSP spot market price stdev / " << spotMarketSD << " /;" << endl;
		output << endl;
		output << "Parameter" << endl;
		output << "zSP, pdfSP, cdfSp, delta, counter, check, gapAllowed, iterationCountLimit" << endl;
		output << "sumOrders, sumSelectedSuppliers; " << endl;
		output << endl;
		output << "gapAllowed= 0.0001;" << endl;
		output << "iterationCountLimit=10;" << endl;
		output << endl;
		output << "zSP = (ec - muSP) / sdSP;" << endl;
		output << "pdfSP =  exp(-0.5*power(zSp,2))/sqrt(2*3.14156);" << endl;
		output << "cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));" << endl;
		output << "If (zSP >= 0, " << endl;
		output << "         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);" << endl;
		output << "else  " << endl;
		output << "         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);" << endl;
		output << ");" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << endl;
		output << "Variables" << endl;
		output << "z service level " << endl;
		output << "expected expected cost;" << endl;
		output << endl;
		output << "Positive variables" << endl;
		output << "pdf" << endl;
		output << "cdf" << endl;
		output << "loss" << endl;
		output << "TotalQ" << endl;
		output << "q(i)      order quantity for supplier i;" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Equations" << endl;
		output << "objective        define objective function" << endl;
		output << "calcTotalQ    calculate total order for each scenarion" << endl;
		output << "calcLoss      define loss constraint Q" << endl;
		output << "calcz         calculate z" << endl;
		output << "calcpdf" << endl;
		output << "calccdf" << endl;
		output << "LinConst; " << endl;
		output << endl;
		output << "objective..              expected=e= sum((i),c(i)*q(i)*pi(i))+h*(TotalQ-mu) + (h+muSP)* loss;" << endl;
		output << "calcTotalQ..             TotalQ =e=  sum((i), pi(i)*q(i));" << endl;
		output << "calcLoss..               loss =e=sd*(pdf -z*(1-cdf));" << endl;
		output << "calcz..                  TotalQ=e=mu+sd*z; " << endl;
		output << "calcpdf..                pdf*sqrt(2*3.14156)=e=exp(-0.5*power(z,2));" << endl;
		output << "calccdf..                cdf=e= errorf(z);" << endl;
		output << "LinConst(cycle)..        sum((i)$(ord(i) le counter),(A(cycle,i)*q(i))) =e= cdfINV(cycle);" << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "Model AONnew /All/;" << endl;
		output << "option limrow=100;" << endl;
		output << "cycle(iter)=no;" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "******** Updating the left-hand-side : matrix A" << endl;
		output << "loop(iter," << endl;
		output << "cycle(iter-1)$(ord(iter)<>1)=no; " << endl;
		output << "cycle(iter)=yes;" << endl;
		output << "A(cycle, i)$(ord(iter) <> ord(i)) = pi(i);" << endl;
		output << "A(cycle, i)$(ord(iter) eq ord(i)) = 1; " << endl;
		output << "display A;" << endl;
		output << ");" << endl;
		output << "cycle(iter) = no;" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "******** Updating the right-hand-side : parameter" << endl;
		output << "loop(iter," << endl;
		output << "cycle(iter - 1)$(ord(iter)<>1) = no;" << endl;
		output << "cycle(iter) = yes;" << endl;
		output << "         loop(i," << endl;
		output << "                 if (ord(iter) = ord(i), " << endl;
		output << "                    cdfINV(cycle) = cdfINV_temp(i)" << endl;
		output << "                 ));" << endl;
		output << "display cdfINV;" << endl;
		output << ");" << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "********* Solving the problem " << endl;
		output << "parameters " << endl;
		output << endl;
		output << "RealObj" << endl;
		output << "RealObj_iter" << endl;
		output << "RealObj_new" << endl;
		output << "KR" << endl;
		output << "gapPercent" << endl;
		output << "iterationCount" << endl;
		output << "zReal(u) service level for scenario u" << endl;
		output << "zKRReal(u)" << endl;
		output << "pdfReal(u)    pdf for scenario u" << endl;
		output << "cdfReal(u)    cdf for scenario u" << endl;
		output << "TotalQReal(u) total order quantity for scenario u" << endl;
		output << "pdfKRReal(u)  pdf including KR for scenario u " << endl;
		output << "cdfKRReal(u)  cdf including KR for scenario u" << endl;
		output << "lossReal(u)   loss value for Qu scenario u " << endl;
		output << "lossKRReal(u) loss value for Qu + KR for scenario u " << endl;
		output << "costReal(u)   cost scenario u " << endl;
		output << "temp" << endl;
		output << "temp1" << endl;
		output << "temp2" << endl;
		output << "temp3" << endl;
		output << "lowerB        for generating intitial X_i values randomly" << endl;
		output << "upperB        for generating intitial X_i values randomly; " << endl;
		output << endl;
		output << "*---------- Start the iteration loop " << endl;
		output << "RealObj=100000000;  " << endl;
		output << "counter=0; " << endl;
		output << "cycle(iter)=no;" << endl;
		output << endl;
		output << "loop(iter," << endl;
		output << endl;
		output << "RealObj_iter=100000000;" << endl;
		output << "iterationCount=0;" << endl;
		output << "counter=counter+1;" << endl;
		output << "cycle(iter)=yes; " << endl;
		output << "gapPercent=1; " << endl;
		output << endl;
		output << "*----- Find initial values of cdf inverse for the fractional values including X " << endl;
		output << "fractileKR= (delta-rc) / delta;" << endl;
		output << "cdfINVKR =icdfNormalTEMP(fractileKR,mu,sd);" << endl;
		output << "upperB = fractileKR;" << endl;
		output << endl;
		output << "loop(i," << endl;
		output << "         lowerB= fractileKR * pi(i);" << endl;
		output << "         temp1= fractileKR;" << endl;
		output << "         temp2= pi(i)*(muSP-c(i))/delta;" << endl;
		output << "         if (temp1<temp2," << endl;
		output << "                 upperB=temp1;" << endl;
		output << "         else" << endl;
		output << "                 upperB=temp2;" << endl;
		output << "         );  " << endl;
		output << "         partX(i)= Uniform(lowerB,upperB);" << endl;
		output << "         temp = (pi(i) * (muSP - c(i)) - delta * partX(i))/(pi(i)*(h + muSP - delta)); " << endl;
		output << "         cdfINV_partX(i) =  icdfNormalTEMP(temp,mu,sd); " << endl;
		output << "         display partX, cdfINV_partX; " << endl;
		output << ");  " << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "while ((gapPercent>  gapAllowed) or (iterationCount<iterationCountLimit), " << endl;
		output << "loop(k,  " << endl;
		output << "         if (ord(k) le counter, " << endl;
		output << "            cdfINV(cycle(k)) = cdfINV_partX(k)" << endl;
		output << "         );" << endl;
		output << ");" << endl;
		output << "*** Set the values of q_i by solving the Linear equations" << endl;
		output << endl;
		output << "q.lo(i)$(ord(i) le ord(iter))=0;" << endl;
		output << "q.up(i)$(ord(i) le ord(iter))=+inf;" << endl;
		output << "q.fx(i)$(ord(i) > ord(iter))=0;" << endl;
		output << "option nlp=Conopt; " << endl;
		output << "solve  AONnew using nlp minimizing expected; " << endl;
		output << "display q.l; " << endl;
		output << endl;
		output << "KR = cdfINVKR - sum((i), q.l(i)*pi(i));" << endl;
		output << endl;
		output << "display KR,iterationCount;  " << endl;
		output << "*** Calc real objetive value for this solution " << endl;
		output << endl;
		output << "TotalQReal(u)=  sum((i), q.l(i)*R(u,i));" << endl;
		output << endl;
		output << "zReal(u)=(TotalQReal(u)-mu)/sd; " << endl;
		output << "zKRReal(u)=(TotalQReal(u)+KR-mu)/sd; " << endl;
		output << "pdfReal(u) = exp(-0.5*power(zReal(u),2))/sqrt(2*3.14156); " << endl;
		output << "cdfReal(u) = errorf(zReal(u));" << endl;
		output << "Loop(u,  " << endl;
		output << "         lossReal(u) = sd * (pdfReal(u) - zReal(u)* (1-cdfReal(u)));  " << endl;
		output << "); " << endl;
		output << "pdfKRReal(u) = exp(-0.5*power(zKRReal(u), 2)) / sqrt(2 * 3.14156);" << endl;
		output << "cdfKRReal(u)= errorf(zKRReal(u));" << endl;
		output << "Loop(u," << endl;
		output << "         lossKRReal(u) = sd * (pdfKRReal(u) - zKRReal(u)*(1-cdfKRReal(u)));" << endl;
		output << ");" << endl;
		output << endl;
		output << "costReal(u)=  sum((i),c(i)*q.l(i)*R(u,i))+ rc*KR +" << endl;
		output << "                         h*(TotalQReal(u)-mu) + (h+muSP-delta)* lossReal(u)+ " << endl;
		output << "                          delta * lossKRReal(u);  " << endl;
		output << endl;
		output << "RealObj_new= sum((u), Prob(u)*costReal(u));" << endl;
		output << endl;
		output << "If (RealObj_new < RealObj_iter, " << endl;
		output << "         gapPercent = (RealObj_iter-RealObj_new)/RealObj_new;" << endl;
		output << "         RealObj_iter=RealObj_new;" << endl;
		output << "else " << endl;
		output << "gapPercent=0;   " << endl;
		output << ");" << endl;
		output << "iterationCount=iterationCount+1;" << endl;
		output << endl;
		output << "loop(i," << endl;
		output << endl;
		output << "temp1= sum((u)$(R(u,i) eq 1), Prob(u)* cdfKRReal(u)); " << endl;
		output << "temp2= pi(i)*(muSP-c(i))/delta; " << endl;
		output << "temp3= fractileKR;  " << endl;
		output << "if (temp2<temp3,   " << endl;
		output << "        upperB=temp2; " << endl;
		output << "else" << endl;
		output << "        upperB=temp3; " << endl;
		output << "); " << endl;
		output << "lowerB=pi(i)*(delta-rc)/delta;" << endl;
		output << endl;
		output << "if ((temp1<upperB)and(temp1>lowerB)," << endl;
		output << "         partX(i)=temp1;   " << endl;
		output << "elseif (temp1 ge upperB), " << endl;
		output << "         partX(i)=upperB;" << endl;
		output << "else" << endl;
		output << "         partX(i)=lowerB;" << endl;
		output << ");  " << endl;
		output << "temp = (pi(i) * (muSP - c(i)) - delta * partX(i))/(pi(i)*(h + muSP - delta));  " << endl;
		output << "cdfINV_partX(i) =  icdfNormalTEMP(temp,mu,sd); " << endl;
		output << "); " << endl;
		output << "); " << endl;
		output << "If (RealObj_iter < RealObj, " << endl;
		output << "        RealObj=RealObj_iter; " << endl;
		output << "sumOrders = sum((i), q.l(i));    " << endl;
		output << "sumSelectedSuppliers = 0; " << endl;
		output << "loop(i, " << endl;
		output << "   if(q.l(i)>0,  " << endl;
		output << "      sumSelectedSuppliers= sumSelectedSuppliers+1;" << endl;
		output << "   );  " << endl;
		output << ");" << endl;
		output << "); " << endl;
		output << "); " << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Display RealObj;" << endl;
		output << "tcomp = TimeComp;" << endl;
		output << "texec = TimeExec;" << endl;
		output << "telapsed = TimeElapsed;" << endl;
		//output << "*Display tcomp, texec, telapsed;" << endl;
		output << endl;
		output << "file  Results storing the results / results_" + convert_t.str() + "_HeuristicApproxIterative.txt / ; " << endl;
		output << "Results.ap=1;" << endl;
		output << "put Results;" << endl;
		output << "put /Realobj','telapsed','sumSelectedSuppliers','sumOrders/;" << endl;
		output << "putclose Results;" << endl;
		output.close();
	}
	//for (it = Instances->begin(); it < Instances->end(); it++)
	//{
	//	(*it)->~ProblemInstance();
	//}
}

void DataGenerator::WriteHeuristicModelsGamsWithKR_ExactIterative(string FolderName) {

	vector<ProblemInstance*>::iterator it;
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		stringstream convert_t, convert_ins, convert_set;
		convert_t << numSuppliers;
		convert_ins << (*it)->instanceNo;
		string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_HeuristicExactIterative.gms";
		ofstream output;
		output.open(out_file.c_str());
		output << "$FuncLibIn stolib stodclib  " << endl;
		output << "function icdfNormalTEMP     /stolib.icdfNormal/" << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "Scalar" << endl;
		output << "tcomp, texec, telapsed;" << endl;
		output << "*------------------------------------------------------------------------------ -" << endl;
		output << "Sets" << endl;
		output << "i suppliers / 1 *" << numSuppliers << " /" << endl;
		int combinations = pow(2, numSuppliers);
		output << "u scenarios / 1 *" << combinations << " /" << endl;
		output << "alias(i,iter,k);" << endl;
		output << "set cycle(iter);" << endl;
		output << "*------------------------------------------------------------------------------ -" << endl;
		output << "Parameters" << endl;
		output << "c(i) selling price" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "pi(i) probability of each supplier's availability" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "cdfINV_temp(i) inverse of F for critical values" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			double critical_frac = (spotMarketMean - ((double)((*it)->cost_tmp[i]))) / (spotMarketMean + ((double)((*it)->hCostInstance)));
			double cdfinv = DataGenerator::Norm_CDF_INV(critical_frac, demandMean, demandSD);
			output << i + 1 << '\t' << cdfinv << endl;
		}
		output << "/" << endl;
		output << "cdfVal(i)" << endl;
		output << "partX(i) X part of the equations" << endl;
		output << "cdf_partX(i) inverse of F for fractions including X" << endl;
		output << "fractileKR, cdfINVKR, A(iter,i);" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Table" << endl;
		output << "R(u,i)  availability" << endl;
		output << '\t';
		for (int i = 0; i < numSuppliers; i++) {
			output << i + 1 << '\t';
		}
		output << endl;
		for (int ns = 0; ns < possbileScenarios; ns++) {
			output << ns + 1 << '\t';
			for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
				output << matrixOfScenarios[ns][supplier_tmp] << '\t';
			}
			output << endl;
		}
		output << "parameter" << endl;
		output << "Prob(u) Scenario probability" << endl;
		output << "/" << endl;
		for (int i = 0; i < possbileScenarios; i++) {
			output << i + 1 << '\t' << instanceProbabilities[i] << endl;
		}
		output << "/;" << endl;
		output << endl;
		output << "Scalars" << endl;
		output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
		output << "rc reservation cost at the backup supplier /" << rc << " /" << endl;
		output << "ec execution cost at the backup supplier /" << (*it)->ec_tmp << " /" << endl;
		output << "mu demand mean / " << demandMean << " /" << endl;
		output << "sd demand stdev / " << demandSD << " /" << endl;
		output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
		output << "sdSP spot market price stdev / " << spotMarketSD << " /;" << endl;
		output << endl;
		output << "Parameter" << endl;
		output << "zSP, pdfSP, cdfSp, delta, counter, check, gapAllowed, iterationCountLimit;" << endl;
		output << endl;
		output << "gapAllowed= 0.0001;" << endl;
		output << "iterationCountLimit=10;" << endl;
		output << endl;
		output << "zSP = (ec - muSP) / sdSP;" << endl;
		output << "pdfSP =  exp(-0.5*power(zSp,2))/sqrt(2*3.14156);" << endl;
		output << "cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));" << endl;
		output << "If (zSP >= 0, " << endl;
		output << "         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);" << endl;
		output << "else  " << endl;
		output << "         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);" << endl;
		output << ");" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << endl;
		output << "Variables" << endl;
		output << "z(u) service level for scenario u " << endl;
		output << "expected expected cost;" << endl;
		output << endl;
		output << "Positive variables" << endl;
		output << "pdf(u)    pdf for scenario u" << endl;
		output << "cdf(u)    cdf for scenario u" << endl;
		output << "TotalQ(u) total order quantity for scenario u" << endl;
		output << "loss(u)   loss value for Qu scenario u " << endl;
		output << "cost(u)   cost scenario u" << endl;
		output << "q(i)      order quantity for supplier i;" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Equations" << endl;
		output << "objective        define objective function" << endl;
		output << "scenarioCost(u)  scenario cost definitio " << endl;
		output << "calcTotalQ(u)    calculate total order for each scenarion " << endl;
		output << "calcLoss(u)      define loss constraint Qu scenario u " << endl;
		output << "calcpdf(u)       calculate pdf of Qu for scenario u" << endl;
		output << "calccdf(u)       calculate cdf of Qu for scenario u" << endl;
		output << "calcz(u)         calculate z(u)" << endl;
		output << "KKTconst(i); " << endl;
		output << endl;
		output << "objective..              expected=e= sum((u), Prob(u)*cost(u));" << endl;
		output << "scenarioCost(u)..        cost(u)=e=  sum((i),c(i)*q(i)*R(u,i))+" << endl;
		output << "                         h*(TotalQ(u)-mu) + (h+muSP)* loss(u);" << endl;
		output << "calcTotalQ(u)..          TotalQ(u) =e=  sum((i), q(i)*R(u,i)); " << endl;
		output << "calcLoss(u)..            loss(u)=e=sd*(pdf(u) -z(u)*(1-cdf(u)));" << endl;
		output << "calcpdf(u)..             pdf(u)*sqrt(2*3.14156)=e=exp(-0.5*power(z(u),2));" << endl;
		output << "calccdf(u)..             cdf(u)=e= errorf(z(u));" << endl;
		output << "calcz(u)..               TotalQ(u)=e=mu+sd*z(u); " << endl;
		output << "KKTconst(i)$(ord(i) le counter)..            sum((u)$(R(u,i) eq 1),(cdf(u)*Prob(u))) =e= cdfVal(i);" << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "Model AONnew /All/;" << endl;
		output << "option limrow=100;" << endl;
		output << "cycle(iter)=no;" << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "********* Solving the problem " << endl;
		output << "parameters " << endl;
		output << endl;
		output << "RealObj" << endl;
		output << "RealObj_iter" << endl;
		output << "RealObj_new" << endl;
		output << "KR" << endl;
		output << "gapPercent" << endl;
		output << "iterationCount" << endl;
		output << "zReal(u) service level for scenario u" << endl;
		output << "zKRReal(u)" << endl;
		output << "pdfReal(u)    pdf for scenario u" << endl;
		output << "cdfReal(u)    cdf for scenario u" << endl;
		output << "TotalQReal(u) total order quantity for scenario u" << endl;
		output << "pdfKRReal(u)  pdf including KR for scenario u " << endl;
		output << "cdfKRReal(u)  cdf including KR for scenario u" << endl;
		output << "lossReal(u)   loss value for Qu scenario u " << endl;
		output << "lossKRReal(u) loss value for Qu + KR for scenario u " << endl;
		output << "costReal(u)   cost scenario u " << endl;
		output << "temp" << endl;
		output << "temp1" << endl;
		output << "temp2" << endl;
		output << "temp3" << endl;
		output << "lowerB        for generating intitial X_i values randomly" << endl;
		output << "upperB        for generating intitial X_i values randomly; " << endl;
		output << endl;
		output << "*---------- Start the iteration loop " << endl;
		output << "RealObj=100000000;  " << endl;
		output << "counter=0; " << endl;
		output << "cycle(iter)=no;" << endl;
		output << endl;
		output << "loop(iter," << endl;
		output << endl;
		output << "RealObj_iter=100000000;" << endl;
		output << "iterationCount=0;" << endl;
		output << "counter=counter+1;" << endl;
		output << "cycle(iter)=yes; " << endl;
		output << "gapPercent=1; " << endl;
		output << endl;
		output << "*----- Find initial values of cdf inverse for the fractional values including X " << endl;
		output << "fractileKR= (delta-rc) / delta;" << endl;
		output << "cdfINVKR =icdfNormalTEMP(fractileKR,mu,sd);" << endl;
		output << "upperB = fractileKR;" << endl;
		output << endl;
		output << "loop(i," << endl;
		output << "         lowerB= fractileKR * pi(i);" << endl;
		output << "         temp1= fractileKR;" << endl;
		output << "         temp2= pi(i)*(muSP-c(i))/delta;" << endl;
		output << "         if (temp1<temp2," << endl;
		output << "                 upperB=temp1;" << endl;
		output << "         else" << endl;
		output << "                 upperB=temp2;" << endl;
		output << "         );  " << endl;
		output << "         partX(i)= Uniform(lowerB,upperB);" << endl;
		output << "         cdf_partX(i) = (pi(i) * (muSP - c(i)) - delta * partX(i))/(h + muSP - delta);" << endl;
		output << ");  " << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "while ((gapPercent>  gapAllowed) or (iterationCount<iterationCountLimit), " << endl;
		output << "loop(k,  " << endl;
		output << "            cdfVal(k) = cdf_partX(k)  " << endl;
		output << ");" << endl;
		output << "*** Set the values of q_i by solving the Linear equations" << endl;
		output << endl;
		output << "q.lo(i)$(ord(i) le ord(iter))=0;" << endl;
		output << "q.up(i)$(ord(i) le ord(iter))=+inf;" << endl;
		output << "q.fx(i)$(ord(i) > ord(iter))=0;" << endl;
		output << "option nlp=Conopt; " << endl;
		output << "solve  AONnew using nlp minimizing expected; " << endl;
		output << "display q.l; " << endl;
		output << endl;
		output << "KR = cdfINVKR - sum((i), q.l(i)*pi(i));" << endl;
		output << endl;
		output << "display KR,iterationCount;  " << endl;
		output << "*** Calc real objetive value for this solution " << endl;
		output << endl;
		output << "TotalQReal(u)=  sum((i), q.l(i)*R(u,i));" << endl;
		output << endl;
		output << "zReal(u)=(TotalQReal(u)-mu)/sd; " << endl;
		output << "zKRReal(u)=(TotalQReal(u)+KR-mu)/sd; " << endl;
		output << "pdfReal(u) = exp(-0.5*power(zReal(u),2))/sqrt(2*3.14156); " << endl;
		output << "cdfReal(u) = errorf(zReal(u));" << endl;
		output << "Loop(u,  " << endl;
		output << "         lossReal(u) = sd * (pdfReal(u) - zReal(u)* (1-cdfReal(u)));  " << endl;
		output << "); " << endl;
		output << "pdfKRReal(u) = exp(-0.5*power(zKRReal(u), 2)) / sqrt(2 * 3.14156);" << endl;
		output << "cdfKRReal(u)= errorf(zKRReal(u));" << endl;
		output << "Loop(u," << endl;
		output << "         lossKRReal(u) = sd * (pdfKRReal(u) - zKRReal(u)*(1-cdfKRReal(u)));" << endl;
		output << ");" << endl;
		output << endl;
		output << "costReal(u)=  sum((i),c(i)*q.l(i)*R(u,i))+ rc*KR +" << endl;
		output << "                         h*(TotalQReal(u)-mu) + (h+muSP-delta)* lossReal(u)+ " << endl;
		output << "                          delta * lossKRReal(u);  " << endl;
		output << endl;
		output << "RealObj_new= sum((u), Prob(u)*costReal(u));" << endl;
		output << endl;
		output << "If (RealObj_new < RealObj_iter, " << endl;
		output << "         gapPercent = (RealObj_iter-RealObj_new)/RealObj_new;" << endl;
		output << "         RealObj_iter=RealObj_new;" << endl;
		output << "else " << endl;
		output << "gapPercent=0;   " << endl;
		output << ");" << endl;
		output << "iterationCount=iterationCount+1;" << endl;
		output << endl;
		output << "loop(i," << endl;
		output << endl;
		output << "temp1= sum((u)$(R(u,i) eq 1), Prob(u)* cdfKRReal(u)); " << endl;
		output << "temp2= pi(i)*(muSP-c(i))/delta; " << endl;
		output << "temp3= fractileKR;  " << endl;
		output << "if (temp2<temp3,   " << endl;
		output << "        upperB=temp2; " << endl;
		output << "else" << endl;
		output << "        upperB=temp3; " << endl;
		output << "); " << endl;
		output << "lowerB=pi(i)*(delta-rc)/delta;" << endl;
		output << endl;
		output << "if ((temp1<upperB)and(temp1>lowerB)," << endl;
		output << "         partX(i)=temp1;   " << endl;
		output << "elseif (temp1 ge upperB), " << endl;
		output << "         partX(i)=upperB;" << endl;
		output << "else" << endl;
		output << "         partX(i)=lowerB;" << endl;
		output << ");  " << endl;
		output << "cdf_partX(i) =(pi(i) * (muSP - c(i)) - delta * partX(i))/(h + muSP - delta);" << endl;
		output << "); " << endl;
		output << "); " << endl;
		output << "If (RealObj_iter < RealObj, " << endl;
		output << "        RealObj=RealObj_iter; " << endl;
		output << "); " << endl;
		output << "); " << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Display RealObj;" << endl;
		output << "tcomp = TimeComp;" << endl;
		output << "texec = TimeExec;" << endl;
		output << "telapsed = TimeElapsed;" << endl;
		//output << "*Display tcomp, texec, telapsed;" << endl;
		output << endl;
		output << "file  Results storing the results / results_" + convert_t.str() + "_" + convert_ins.str() + "_HeuristicExactIterative.txt / ; " << endl;
		output << "Results.ap=1;" << endl;
		output << "put Results;" << endl;
		output << "put /Realobj','telapsed/;" << endl;
		output << "putclose Results;" << endl;
		output.close();
	}
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		(*it)->~ProblemInstance();
	}
}

void DataGenerator::WriteApproxHeuristicGAMSBatchFile(string FolderName)
{
	stringstream convert_temp;
	convert_temp << numSuppliers;
	string out_file = FolderName + "/Batch_" + convert_temp.str() + "_ApproxHeuristic.gms";
	ofstream output;
	output.open(out_file.c_str());
	//vector<ProblemInstance*>::iterator it;
	for (int j = 0; j <numInstances; j++)
	{
		output << "$call gams " << numSuppliers << "_" << j << "_HeuristicApproxIterative" << endl;
	}
	output.close();
}


// First set of tests , setting KR equal to zero

void DataGenerator::WriteExactModelsInGamsNoKRIPOPT(string FolderName)
{
	vector<ProblemInstance*>::iterator it;
	DataGenerator::GenerateAllPossibleScenarios();
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		stringstream convert_t, convert_ins, convert_set;
		convert_t << numSuppliers;
		convert_ins << (*it)->instanceNo;
		string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_Exact_NoKR.gms";
		ofstream output;
		output.open(out_file.c_str());
		output << "Sets" << endl;
		output << "i suppliers /1*" << numSuppliers << "/" << endl;
		output << "u scenarios /1*" << possbileScenarios << "/" << endl;
		output << endl;
		output << "Parameters" << endl;
		output << "c(i) selling price" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "pi(i) probability of each supplier's availability" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "Prob(u) Scenario probability" << endl;
		output << "/" << endl;
		for (int i = 0; i < possbileScenarios; i++) {
			output << i + 1 << '\t' << instanceProbabilities[i] << endl;
		}
		output << "/" << endl;
		output << "Table" << endl;
		output << "R(u,i)  availability" << endl;
		output << '\t';
		for (int i = 0; i < numSuppliers; i++) {
			output << i + 1 << '\t';
		}
		output << endl;
		for (int ns = 0; ns < possbileScenarios; ns++) {
			output << ns + 1 << '\t';
			for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
				output << matrixOfScenarios[ns][supplier_tmp] << '\t';
			}
			output << endl;
		}
		output << endl;
		output << "Scalars" << endl;
		output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
		output << "mu demand mean / " << demandMean << " /" << endl;
		output << "sd demand stdev / " << demandSD << " /" << endl;
		output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
		output << "sdSP spot market price stdev / " << spotMarketSD << " /" << endl;
		output << "telapsed;" << endl;
		output << endl;
		output << endl;
		output << "Variables" << endl;
		output << "z(u) service level for scenario u" << endl;
		output << "expected expected profit;" << endl;
		output << endl;
		output << "Positive variables" << endl;
		output << "pdf(u)    pdf for scenario u" << endl;
		output << "cdf(u)    cdf for scenario u" << endl;
		output << "TotalQ(u) total order quantity for scenario u" << endl;
		output << "loss(u)   loss value for Qu scenario u" << endl;
		output << "cost(u)   cost scenario u" << endl;
		output << "q(i)      order quantity for supplier i;" << endl;
		output << endl;
		output << "Equations" << endl;
		output << "objective        define objective function" << endl;
		output << "scenarioCost(u)  scenario cost definition" << endl;
		output << "calcTotalQ(u)    calculate total order for each scenarion" << endl;
		output << "calcLoss(u)      define loss constraint Qu scenario u" << endl;
		output << "calcpdf(u)       calculate pdf of Qu for scenario u" << endl;
		output << "calccdf(u)       calculate cdf of Qu for scenario u" << endl;
		output << "calcz(u)         calculate z(u);" << endl;
		output << endl;
		output << "objective..              expected=e= sum((u), Prob(u)*cost(u));" << endl;
		output << "scenarioCost(u)..        cost(u)=e=  sum((i),c(i)*q(i)*R(u,i))+ " << endl;
		output << "                         h*(TotalQ(u)-mu) + (h+muSP)* loss(u); " << endl;
		output << "calcTotalQ(u)..          TotalQ(u) =e=  sum((i), q(i)*R(u,i));" << endl;
		output << "calcLoss(u)..            loss(u)=e=sd*(pdf(u) -z(u)*(1-cdf(u)));" << endl;
		output << "calcpdf(u)..             pdf(u)*sqrt(2*3.14156)=e=exp(-0.5*power(z(u),2));" << endl;
		output << "calccdf(u)..             cdf(u)=e= errorf(z(u));" << endl;
		output << "calcz(u)..               TotalQ(u)=e=mu+sd*z(u); " << endl;
		output << "z.up(u)=20; " << endl;
		output << "z.lo(u)=-20;" << endl;
		output << endl;
		output << "Model AONnew / all / ;" << endl;
		output << "option optcr = 1e-5;" << endl;
		output << "option nlp=IPOPT; " << endl;
		output << "option iterlim = 100000000;" << endl;
		output << "option reslim = 5000;" << endl;
		output << "option domlim = 1000;" << endl;
		output << "solve  AONnew using nlp minimizing expected;" << endl;
		output << endl;	
		output << endl;
		output << "Display AONnew.modelstat;" << endl;
		output << "Display AONnew.solvestat;" << endl;
		output << "telapsed = TimeElapsed;" << endl;
		output << endl;
		output << "file  Results storing the results / results_" + convert_t.str() + "_Exact_IPOPT.txt / ; " << endl;
		output << "Results.ap=1;" << endl;
		output << "put Results;" << endl;
		output << "*put /'Objective,Elapsed time'/;" << endl;
		output << "put /expected.l','telapsed','AONnew.solvestat','AONnew.modelstat/;" << endl;
		output << "putclose Results;" << endl;
		output.close();

	}

}

void DataGenerator::WriteExactModelsInGamsNoKRConopt(string FolderName)
{
	vector<ProblemInstance*>::iterator it;
	DataGenerator::GenerateAllPossibleScenarios();
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		stringstream convert_t, convert_ins, convert_set;
		convert_t << numSuppliers;
		convert_ins << (*it)->instanceNo;
		string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_Exact_NoKR.gms";
		ofstream output;
		output.open(out_file.c_str());
		output << "Sets" << endl;
		output << "i suppliers /1*" << numSuppliers << "/" << endl;
		output << "u scenarios /1*" << possbileScenarios << "/" << endl;
		output << endl;
		output << "Parameters" << endl;
		output << "c(i) selling price" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "pi(i) probability of each supplier's availability" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "Prob(u) Scenario probability" << endl;
		output << "/" << endl;
		for (int i = 0; i < possbileScenarios; i++) {
			output << i + 1 << '\t' << instanceProbabilities[i] << endl;
		}
		output << "/" << endl;
		output << "Table" << endl;
		output << "R(u,i)  availability" << endl;
		output << '\t';
		for (int i = 0; i < numSuppliers; i++) {
			output << i + 1 << '\t';
		}
		output << endl;
		for (int ns = 0; ns < possbileScenarios; ns++) {
			output << ns + 1 << '\t';
			for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
				output << matrixOfScenarios[ns][supplier_tmp] << '\t';
			}
			output << endl;
		}
		output << endl;
		output << "Scalars" << endl;
		output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
		output << "mu demand mean / " << demandMean << " /" << endl;
		output << "sd demand stdev / " << demandSD << " /" << endl;
		output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
		output << "sdSP spot market price stdev / " << spotMarketSD << " /" << endl;
		output << "telapsed;" << endl;
		output << endl;
		output << endl;
		output << "Variables" << endl;
		output << "z(u) service level for scenario u" << endl;
		output << "expected expected profit;" << endl;
		output << endl;
		output << "Positive variables" << endl;
		output << "pdf(u)    pdf for scenario u" << endl;
		output << "cdf(u)    cdf for scenario u" << endl;
		output << "TotalQ(u) total order quantity for scenario u" << endl;
		output << "loss(u)   loss value for Qu scenario u" << endl;
		output << "cost(u)   cost scenario u" << endl;
		output << "q(i)      order quantity for supplier i;" << endl;
		output << endl;
		output << "Equations" << endl;
		output << "objective        define objective function" << endl;
		output << "scenarioCost(u)  scenario cost definition" << endl;
		output << "calcTotalQ(u)    calculate total order for each scenarion" << endl;
		output << "calcLoss(u)      define loss constraint Qu scenario u" << endl;
		output << "calcpdf(u)       calculate pdf of Qu for scenario u" << endl;
		output << "calccdf(u)       calculate cdf of Qu for scenario u" << endl;
		output << "calcz(u)         calculate z(u);" << endl;
		output << endl;
		output << "objective..              expected=e= sum((u), Prob(u)*cost(u));" << endl;
		output << "scenarioCost(u)..        cost(u)=e=  sum((i),c(i)*q(i)*R(u,i))+ " << endl;
		output << "                         h*(TotalQ(u)-mu) + (h+muSP)* loss(u); " << endl;
		output << "calcTotalQ(u)..          TotalQ(u) =e=  sum((i), q(i)*R(u,i));" << endl;
		output << "calcLoss(u)..            loss(u)=e=sd*(pdf(u) -z(u)*(1-cdf(u)));" << endl;
		output << "calcpdf(u)..             pdf(u)*sqrt(2*3.14156)=e=exp(-0.5*power(z(u),2));" << endl;
		output << "calccdf(u)..             cdf(u)=e= errorf(z(u));" << endl;
		output << "calcz(u)..               TotalQ(u)=e=mu+sd*z(u); " << endl;
		output << "z.up(u)=20; " << endl;
		output << "z.lo(u)=-20;" << endl;
		output << endl;
		output << "Model AONnew / all / ;" << endl;
		output << "option optcr = 1e-5;" << endl;
		output << "option nlp=Conopt; " << endl;
		output << "option iterlim = 100000000;" << endl;
		output << "option reslim = 5000;" << endl;
		output << "option domlim = 1000;" << endl;
		output << "solve  AONnew using nlp minimizing expected;" << endl;
		output << endl;
		output << endl;
		output << "Display AONnew.modelstat;" << endl;
		output << "Display AONnew.solvestat;" << endl;
		output << "telapsed = TimeElapsed;" << endl;
		output << endl;
		output << "file  Results storing the results / results_" + convert_t.str() + "_Exact_Conopt.txt / ; " << endl;
		output << "Results.ap=1;" << endl;
		output << "put Results;" << endl;
		output << "*put /'Objective,Elapsed time'/;" << endl;
		output << "put /expected.l','telapsed','AONnew.solvestat','AONnew.modelstat/;" << endl;
		output << "putclose Results;" << endl;
		output.close();

	}

}

void DataGenerator::WriteExactModelsInGamsNoKRLindoglobal(string FolderName)
{
	vector<ProblemInstance*>::iterator it;
	DataGenerator::GenerateAllPossibleScenarios();
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		stringstream convert_t, convert_ins, convert_set;
		convert_t << numSuppliers;
		convert_ins << (*it)->instanceNo;
		string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_Exact_NoKR.gms";
		ofstream output;
		output.open(out_file.c_str());
		output << "Sets" << endl;
		output << "i suppliers /1*" << numSuppliers << "/" << endl;
		output << "u scenarios /1*" << possbileScenarios << "/" << endl;
		output << endl;
		output << "Parameters" << endl;
		output << "c(i) selling price" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "pi(i) probability of each supplier's availability" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "Prob(u) Scenario probability" << endl;
		output << "/" << endl;
		for (int i = 0; i < possbileScenarios; i++) {
			output << i + 1 << '\t' << instanceProbabilities[i] << endl;
		}
		output << "/" << endl;
		output << "Table" << endl;
		output << "R(u,i)  availability" << endl;
		output << '\t';
		for (int i = 0; i < numSuppliers; i++) {
			output << i + 1 << '\t';
		}
		output << endl;
		for (int ns = 0; ns < possbileScenarios; ns++) {
			output << ns + 1 << '\t';
			for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
				output << matrixOfScenarios[ns][supplier_tmp] << '\t';
			}
			output << endl;
		}
		output << endl;
		output << "Scalars" << endl;
		output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
		output << "mu demand mean / " << demandMean << " /" << endl;
		output << "sd demand stdev / " << demandSD << " /" << endl;
		output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
		output << "sdSP spot market price stdev / " << spotMarketSD << " /" << endl;
		output << "telapsed;" << endl;
		output << endl;
		output << endl;
		output << "Variables" << endl;
		output << "z(u) service level for scenario u" << endl;
		output << "expected expected profit;" << endl;
		output << endl;
		output << "Positive variables" << endl;
		output << "pdf(u)    pdf for scenario u" << endl;
		output << "cdf(u)    cdf for scenario u" << endl;
		output << "TotalQ(u) total order quantity for scenario u" << endl;
		output << "loss(u)   loss value for Qu scenario u" << endl;
		output << "cost(u)   cost scenario u" << endl;
		output << "q(i)      order quantity for supplier i;" << endl;
		output << endl;
		output << "Equations" << endl;
		output << "objective        define objective function" << endl;
		output << "scenarioCost(u)  scenario cost definition" << endl;
		output << "calcTotalQ(u)    calculate total order for each scenarion" << endl;
		output << "calcLoss(u)      define loss constraint Qu scenario u" << endl;
		output << "calcpdf(u)       calculate pdf of Qu for scenario u" << endl;
		output << "calccdf(u)       calculate cdf of Qu for scenario u" << endl;
		output << "calcz(u)         calculate z(u);" << endl;
		output << endl;
		output << "objective..              expected=e= sum((u), Prob(u)*cost(u));" << endl;
		output << "scenarioCost(u)..        cost(u)=e=  sum((i),c(i)*q(i)*R(u,i))+ " << endl;
		output << "                         h*(TotalQ(u)-mu) + (h+muSP)* loss(u); " << endl;
		output << "calcTotalQ(u)..          TotalQ(u) =e=  sum((i), q(i)*R(u,i));" << endl;
		output << "calcLoss(u)..            loss(u)=e=sd*(pdf(u) -z(u)*(1-cdf(u)));" << endl;
		output << "calcpdf(u)..             pdf(u)*sqrt(2*3.14156)=e=exp(-0.5*power(z(u),2));" << endl;
		output << "calccdf(u)..             cdf(u)=e= errorf(z(u));" << endl;
		output << "calcz(u)..               TotalQ(u)=e=mu+sd*z(u); " << endl;
		output << "z.up(u)=20; " << endl;
		output << "z.lo(u)=-20;" << endl;
		output << endl;
		output << "Model AONnew / all / ;" << endl;
		output << "option optcr = 1e-5;" << endl;
		output << "option nlp=lindoglobal; " << endl;
		output << "option iterlim = 100000000;" << endl;
		output << "option reslim = 5000;" << endl;
		output << "option domlim = 1000;" << endl;
		output << "solve  AONnew using nlp minimizing expected;" << endl;
		output << endl;
		output << endl;
		output << "Display AONnew.modelstat;" << endl;
		output << "Display AONnew.solvestat;" << endl;
		output << "telapsed = TimeElapsed;" << endl;
		output << endl;
		output << "file  Results storing the results / results_" + convert_t.str() + "_Exact_LG.txt / ; " << endl;
		output << "Results.ap=1;" << endl;
		output << "put Results;" << endl;
		output << "*put /'Objective,Elapsed time'/;" << endl;
		output << "put /expected.l','telapsed','AONnew.solvestat','AONnew.modelstat/;" << endl;
		output << "putclose Results;" << endl;
		output.close();

	}

}


void DataGenerator::WriteExactGAMSBatchFileNoKR(string FolderName) {
	stringstream convert_temp;
	convert_temp << numSuppliers;
	string out_file = FolderName + "/Batch_" + convert_temp.str() + "_Exact.gms";
	ofstream output;
	output.open(out_file.c_str());
	//vector<ProblemInstance*>::iterator it;
	for (int j = 0; j <numInstances; j++)
	{
		output << "$call gams " << numSuppliers << "_" << j << "_Exact_NoKR" << endl;
	}
	output.close();
}

void DataGenerator::WriteHeuristicModelsGamsNoKR(string FolderName) {

	// this is incomplete
	vector<ProblemInstance*>::iterator it;
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		stringstream convert_t, convert_ins, convert_set;
		convert_t << numSuppliers;
		convert_ins << (*it)->instanceNo;
		string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_Heuristic.gms";
		ofstream output;
		output.open(out_file.c_str());
		output << "Scalar" << endl;
		output << "tcomp, texec, telapsed;" << endl;
		output << "Sets" << endl;
		output << "i suppliers / 1 *" << numSuppliers << " /" << endl;
		output << "iter iterator / 1 *" << numSuppliers << " /" << endl;
		int combinations = pow(2, numSuppliers);
		output << "u scenarios / 1 *" << combinations << " /" << endl;
		output << "cycle(iter);" << endl;
		output << "alias(i, j);" << endl;
		output << "parameter A(iter,i);" << endl;
		output << "*------------------------------------------------------------------------------ -" << endl;
		output << "Parameters" << endl;		
		output << "c(i) selling price" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "pi(i) probability of each supplier's availability" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
		}
		output << "/" << endl;
		output << "cr(i) critical ratio " << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			double critical_frac = (spotMarketMean - ((double)((*it)->cost_tmp[i]))) / (spotMarketMean + ((double)((*it)->hCostInstance)));
			output << i + 1 << '\t' << critical_frac << endl;
		}
		output << "/" << endl;
		output << "cdfINV_temp(i) inverse of F for critical values" << endl;
		output << "/" << endl;
		for (int i = 0; i < numSuppliers; i++)
		{
			double critical_frac = (spotMarketMean - ((double)((*it)->cost_tmp[i]))) / (spotMarketMean + ((double)((*it)->hCostInstance)));
			double cdfinv = DataGenerator::Norm_CDF_INV(critical_frac, demandMean, demandSD);
			output << i + 1 << '\t' << cdfinv << endl;
		}
		output << "/" << endl;
		output << "cdfINV(iter)" << endl;
		output << "*------------------------------------------------------------------------------ -" << endl;
		output << "Table" << endl;
		output << "R(u,i)  availability" << endl;
		output << '\t';
		for (int i = 0; i < numSuppliers; i++) {
			output << i + 1 << '\t';
		}
		output << endl;
		for (int ns = 0; ns < possbileScenarios; ns++) {
			output << ns + 1 << '\t';
			for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
				output << matrixOfScenarios[ns][supplier_tmp] << '\t';
			}
			output << endl;
		}
		output << endl;
		output << "parameter" << endl;
		output << "Prob(u) Scenario probability" << endl;
		output << "/" << endl;
		for (int i = 0; i < possbileScenarios; i++) {
			output << i + 1 << '\t' << instanceProbabilities[i] << endl;
		}
		output << "/;" << endl;		
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Scalars" << endl;
		output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
		output << "mu demand mean / " << demandMean << " /" << endl;
		output << "sd demand stdev / " << demandSD << " /" << endl;
		output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
		output << "sdSP spot market price stdev / " << spotMarketSD << " /;" << endl;
		output << endl;
		output << "Parameter" << endl;
		output << "counter, check;" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << endl;
		output << "Variables" << endl;
		output << "z service level" << endl;
		output << "expected expected cost;" << endl;
		output << endl;
		output << "Positive variables" << endl;
		output << "pdf" << endl;
		output << "cdf" << endl;
		output << "TotalQ total order quantity for scenario u" << endl;
		output << "loss" << endl;
		output << "q(i)      order quantity for supplier i;" << endl;
		output << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Equations" << endl;
		output << "objective     define objective function" << endl;
		output << "calcTotalQ    calculate total order for each scenarion" << endl;
		output << "calcLoss      define loss constraint Q" << endl;
		output << "calcz         calculate z" << endl;
		output << "calcpdf" << endl;
		output << "calccdf" << endl;
		output << "LinConst;" << endl;
		output << endl;
		output << "objective..              expected=e= sum((i),c(i)*q(i)*pi(i))+h*(TotalQ-mu) + (h+muSP)* loss;" << endl;
		output << "calcTotalQ..             TotalQ =e=  sum((i), pi(i)*q(i));" << endl;
		output << "calcLoss..               loss =e= sd * (pdf - z*(1-cdf)); " << endl;
		output << "calcz..                  TotalQ=e=mu+sd*z;" << endl;
		output << "calcpdf..                pdf*sqrt(2*3.14156)=e=exp(-0.5*power(z,2)); " << endl;
		output << "calccdf..                cdf=e= errorf(z);" << endl;
		output << "LinConst(cycle)..        sum((i)$(ord(i) le counter),(A(cycle,i)*q(i))) =e= cdfINV(cycle);" << endl;
		output << "*------------------------------------------------------------------------------- " << endl;
		output << "Model AONnew /All/;" << endl;
		output << "option limrow=100;" << endl;
		output << "cycle(iter)=no;" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "******** Updating the left - hand - side : matrix A" << endl;
		output << "loop(iter," << endl;
		output << "cycle(iter - 1)$(ord(iter)<>1) = no;" << endl;
		output << "cycle(iter) = yes;" << endl;
		output << endl;
		output << "A(cycle, i)$(ord(iter) <> ord(i)) = pi(i);" << endl;
		output << "A(cycle, i)$(ord(iter) eq ord(i)) = 1;" << endl;
		output << "Display A;" << endl;
		output << ");" << endl;
		output << "cycle(iter) = no;" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "******** Updating the right-hand-side : parameter " << endl;
		output << "loop(iter," << endl;
		output << "cycle(iter - 1)$(ord(iter)<>1) = no;" << endl;
		output << "cycle(iter) = yes;" << endl;
		output << "		loop(i," << endl;
		output << "			if (ord(iter) = ord(i)," << endl;
		output << "			cdfINV(cycle) = cdfINV_temp(i)" << endl;
		output << "		)); " << endl;
		output << "display cdfINV;" << endl;
		output << ");" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "********* Solving the problem " << endl;
		output << "Parameters" << endl;
		output << endl;
		output << "RealObj  " << endl;
		output << "RealObj_temp  " << endl;
		output << "zReal(u)			service level for scenario u" << endl;
		output << "pdfReal(u)		pdf for scenario u" << endl;
		output << "cdfReal(u)		cdf for scenario u" << endl;
		output << "TotalQReal(u)	total order quantity for scenario u" << endl;
		output << "lossReal(u)		loss value for Qu scenario u" << endl;
		output << "costReal(u)   cost scenario u;" << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "RealObj=100000000; " << endl;
		output << "counter=0;" << endl;
		output << "cycle(iter) = no;" << endl;
		output << "loop(iter," << endl;
		output << "		counter=counter+1;" << endl;
		output << "		cycle(iter) = yes;" << endl;
		output << "		q.lo(i)$(ord(i) le ord(iter)) = 0;" << endl;
		output << "		q.up(i)$(ord(i) le ord(iter)) = +inf;" << endl;
		output << "		q.fx(i)$(ord(i)> ord(iter)) = 0;" << endl;
		output << "		option nlp=Conopt;  " << endl;
		output << "		option optcr=0; " << endl;
		output << "		solve  AONnew using nlp minimizing expected;  " << endl;
		output << "		display q.l;" << endl;
		output << "TotalQReal(u)=  sum((i), q.l(i)*R(u,i));" << endl;
		output << "zReal(u)=(TotalQReal(u)-mu)/sd;" << endl;
		output << "pdfReal(u) = exp(-0.5*power(zReal(u), 2)) / sqrt(2 * 3.14156);" << endl;
		output << "cdfReal(u) =  errorf(zReal(u));" << endl;
		output << "Loop(u," << endl;
		output << "lossReal(u) = sd * (pdfReal(u) - zReal(u)* (1-cdfReal(u)));" << endl;
		output << ");" << endl;
		output << endl;
		output << "costReal(u)=  sum((i),c(i)*q.l(i)*R(u,i))+ h*(TotalQReal(u)-mu) + (h+muSP)* lossReal(u); " << endl;
		output << "RealObj_temp= sum((u), Prob(u)*costReal(u)); " << endl;
		output << "Display RealObj_temp; " << endl;
		output << "" << endl;
		output << "If (RealObj_temp < RealObj," << endl;
		output << "         RealObj=RealObj_temp;" << endl;
		output << ");   " << endl;
		output << ");   " << endl;
		output << "*-------------------------------------------------------------------------------" << endl;
		output << "Display RealObj; " << endl;
		output << "tcomp = TimeComp;" << endl;
		output << "texec = TimeExec;" << endl;
		output << "telapsed = TimeElapsed;" << endl;
		output << "*Display tcomp, texec, telapsed;" << endl;
		output << endl;
		output << "file  Results storing the results / results_" + convert_t.str() + "_Heuristic.txt / ; " << endl;
		output << "Results.ap=1;" << endl;
		output << "put Results;" << endl;
		output << "put /Realobj','telapsed/;" << endl;
		output << "putclose Results;" << endl;
		output.close();
	}
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		(*it)->~ProblemInstance();
	}
}

void DataGenerator::WriteSAAModelsInGamsNoKR(string FolderName) {
	vector<ProblemInstance*>::iterator it;
	int possbileScenarios = pow(2, numSuppliers);
	for (it = Instances->begin(); it < Instances->end(); it++)
	{
		vector<double> instanceProbabilities = DataGenerator::CalcScenarioProbabilities((*it)->prob_tmp);
		for (int j = 0; j < numSets; j++) {
			stringstream convert_t, convert_ins, convert_set;
			convert_t << numSuppliers;
			convert_ins << (*it)->instanceNo;
			convert_set << j;
			string out_file = FolderName + "/" + convert_t.str() + "_" + convert_ins.str() + "_SAA_Set_" + convert_set.str() + ".gms";
			ofstream output;
			output.open(out_file.c_str());
			output << "Scalar" << endl;
			output << "tcomp, texec, telapsed;" << endl;
			output << "Sets" << endl;
			output << "i suppliers /1*" << numSuppliers << "/" << endl;
			output << "m SAA scenarios /1*" << numScenarios << "/" << endl;
			output << "j   coefficients  /Demand , SP/" << endl;
			output << "u scenarios /1*" << possbileScenarios << "/" << endl;
			output << endl;
			output << "Parameters" << endl;
			output << "c(i) selling price" << endl;
			output << "/" << endl;
			for (int i = 0; i < numSuppliers; i++)
			{
				output << i + 1 << '\t' << (*it)->cost_tmp[i] << endl;
			}
			output << "/" << endl;
			output << "pi(i) probability of each supplier's availability" << endl;
			output << "/" << endl;
			for (int i = 0; i < numSuppliers; i++)
			{
				output << i + 1 << '\t' << (*it)->prob_tmp[i] << endl;
			}
			output << "/" << endl;
			output << "Prob(u) Scenario probability" << endl;
			output << "/" << endl;
			for (int i = 0; i < possbileScenarios; i++) {
				output << i + 1 << '\t' << instanceProbabilities[i] << endl;
			}
			output << "/" << endl;
			output << "Table" << endl;
			output << "RPrime(u,i)  availability" << endl;
			output << '\t';
			for (int i = 0; i < numSuppliers; i++) {
				output << i + 1 << '\t';
			}
			output << endl;
			for (int ns = 0; ns < possbileScenarios; ns++) {
				output << ns + 1 << '\t';
				for (int supplier_tmp = 0; supplier_tmp < numSuppliers; supplier_tmp++) {
					output << matrixOfScenarios[ns][supplier_tmp] << '\t';
				}
				output << endl;
			}
			output << endl;
			output << "Scalars" << endl;
			output << "h cost of unsold item / " << (*it)->hCostInstance << " /" << endl;
			output << "rc reservation cost at the backup supplier /" << rc << " /" << endl;
			output << "ec execution cost at the backup supplier /" << (*it)->ec_tmp << " /" << endl;
			output << "mu demand mean / " << demandMean << " /" << endl;
			output << "sd demand stdev / " << demandSD << " /" << endl;
			output << "muSP spot market price mean / " << spotMarketMean << " /" << endl;
			output << "sdSP spot market price stdev / " << spotMarketSD << " /;" << endl;
			output << endl;
			output << "Parameters " << endl;
			output << "data(m,j) demand and SP values in each scenario " << endl;
			output << "R(m,i) availability of the suppliers in each scenario" << endl;
			output << "sumOrders, sumSelectedSuppliers; " << endl;
			output << endl;
			output << "Variables" << endl;
			output << "expected expected profit;" << endl;
			output << endl;
			output << "Positive variables" << endl;
			output << "TotalQ(m) total order quantity for scenario m " << endl;
			output << "Ipos(m)   amount of overage in scenarion m " << endl;
			output << "Ineg(m)   amount of underage in scenario m " << endl;
			output << "O(m)      amount of the option contract in scenario m" << endl;
			output << "S(m)      spot market purchase in scenario m" << endl;
			output << "cost(m)  cost scenario m" << endl;
			output << "q(i) order quantity for supplier i" << endl;
			output << "KR   capacity reserved at the backup supplier;" << endl;
			output << endl;
			output << "Equations" << endl;
			output << "objective        define objective function" << endl;
			output << "scenarioCost(m)  scenario cost definition " << endl;
			output << "calcTotalQ(m)    calculate total order for each scenarion" << endl;
			output << "overage(m)  " << endl;
			output << "underage(m) " << endl;
			output << "underage2(m)  " << endl;
			output << "optionContract(m)" << endl;
			output << "ConstraintKR;" << endl;
			output << endl;
			output << "objective..         expected=e= (1/" << numScenarios << ") * sum((m),cost(m));" << endl;
			output << "scenarioCost(m)..   cost(m)=e=  sum((i),c(i)*q(i)*R(m,i))+ rc*KR + h*Ipos(m)+ " << endl;
			output << "                         ec* O(m) + data(m, 'SP')*S(m); " << endl;
			output << "calcTotalQ(m)..   TotalQ(m) =e=  sum((i), q(i)*R(m,i)); " << endl;
			output << "overage(m)..      Ipos(m) =g= TotalQ(m) - data(m, 'Demand');" << endl;
			output << "underage(m)..    Ineg(m) =g=  data(m, 'Demand')- TotalQ(m); " << endl;
			output << "underage2(m)..   Ineg(m) =e= O(m)+S(m);  " << endl;
			output << "optionContract(m)..  O(m) =l= KR;" << endl;
			output << "ConstraintKR..       KR =e= 0;" << endl;
			output << endl;
			output << endl;
			output << "Model AONnew / all / ;" << endl;
			output << "$include " << numSuppliers << "_" << (*it)->instanceNo << "_set_" << j << ".inc" << endl;
			output << endl;
			output << endl;
			output << "option optcr = 1e-5;" << endl;
			output << "option iterlim = 100000000;" << endl;
			output << "option reslim = 5000;" << endl;
			output << "option domlim = 20000;" << endl;
			output << "solve  AONnew using lp minimizing expected;" << endl;
			output << endl;
			output << "Parameter" << endl;
			output << "zSP, pdfSP, cdfSp, delta;" << endl;
			output << "zSP = (ec - muSP) / sdSP;" << endl;
			output << "pdfSP =  exp(-0.5*power(zSp,2))/sqrt(2*3.14156);" << endl;
			output << "cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));" << endl;
			output << "If (zSP >= 0, " << endl;
			output << "         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);" << endl;
			output << "else  " << endl;
			output << "         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);" << endl;
			output << ");" << endl;
			output << "Parameter" << endl;
			output << "RealObj" << endl;
			output << "zReal(u) service level for scenario u" << endl;
			output << "zKRReal(u)  " << endl;
			output << "lossReal(u)   loss value for Qu scenario u" << endl;
			output << "lossKRReal(u) loss value for Qu + KR for scenario u" << endl;
			output << "pdfReal(u)    pdf for scenario u" << endl;
			output << "cdfReal(u)    cdf for scenario u " << endl;
			output << "TotalQReal(u) total order quantity for scenario u" << endl;
			output << "pdfKRReal(u)  pdf including KR for scenario u " << endl;
			output << "cdfKRReal(u)  cdf including KR for scenario u" << endl;
			output << "costReal(u)   cost scenario u  ;" << endl;
			output << endl;
			output << "Loop(u," << endl;
			output << "TotalQReal(u)=  sum((i), q.l(i)*RPrime(u,i));" << endl;
			output << "zReal(u)=(TotalQReal(u)-mu)/sd;" << endl;
			output << "zKRReal(u)=(TotalQReal(u)+KR.l-mu)/sd;" << endl;
			output << "pdfReal(u) = exp(-0.5*power(zReal(u),2))/sqrt(2*3.14156);" << endl;
			output << "cdfReal(u) = errorf(zReal(u));" << endl;
			output << "lossReal(u) = sd * (pdfReal(u) - zReal(u)* (1-cdfReal(u)));" << endl;
			output << "pdfKRReal(u) = exp(-0.5*power(zKRReal(u), 2)) / sqrt(2 * 3.14156);" << endl;
			output << "cdfKRReal(u)= errorf(zKRReal(u));" << endl;
			output << endl;
			output << "lossKRReal(u) = sd * (pdfKRReal(u) - zKRReal(u)*(1-cdfKRReal(u)));" << endl;
			output << "costReal(u)=  sum((i),c(i)*q.l(i)*RPrime(u,i))+ rc*KR.l +" << endl;
			output << "                         h*(TotalQReal(u)-mu) + (h+muSP-delta)* lossReal(u)+" << endl;
			output << "                          delta * lossKRReal(u); " << endl;
			output << ");" << endl;
			output << endl;
			output << "RealObj = sum((u), Prob(u)*costReal(u));" << endl;
			output << endl;
			output << "sumOrders = sum((i), q.l(i));    " << endl;
			output << "sumSelectedSuppliers = 0; " << endl;
			output << "loop(i, " << endl;
			output << "   if(q.l(i)>0,  " << endl;
			output << "      sumSelectedSuppliers= sumSelectedSuppliers+1;" << endl;
			output << "   );  " << endl;
			output << ");" << endl;
			output << endl;
			output << "Display AONnew.modelstat;" << endl;
			output << "Display AONnew.solvestat;" << endl;
			output << "tcomp = TimeComp;" << endl;
			output << "texec = TimeExec;" << endl;
			output << "telapsed = TimeElapsed;" << endl;
			output << "*Display tcomp, texec, telapsed;" << endl;
			output << endl;
			output << "file  Results storing the results / results_" + convert_t.str()  + "_SAA.txt / ; " << endl;
			output << "Results.ap=1;" << endl;
			output << "put Results;" << endl;
			//output << "put /expected.l','telapsed/;" << endl;
			output << "put /expected.l','telapsed','RealObj/;" << endl;
			output << "putclose Results;" << endl;
			output.close();

		}

	}
}

void DataGenerator::WriteNoKRHeuristicGAMSBatchFile(string FolderName)
{
	stringstream convert_temp;
	convert_temp << numSuppliers;
	string out_file = FolderName + "/Batch_" + convert_temp.str() + "_NorKRHeuristic.gms";
	ofstream output;
	output.open(out_file.c_str());
	//vector<ProblemInstance*>::iterator it;
	for (int j = 0; j <numInstances; j++)
	{
		output << "$call gams " << numSuppliers << "_" << j << "_Heuristic" << endl;
	}
	output.close();
}