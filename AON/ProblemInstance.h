#ifndef PROBLEMINSTANCE_H
#define PROBLEMINSTANCE_H

class ProblemInstance {
public:
	int instanceNo;
	int numSuppliers_tmp;
	int ec_tmp; //randomly generated execution cost
	int* cost_tmp; // unit cost
	double* prob_tmp; // probability of being available
	int hCostInstance;
	


	ProblemInstance(int instanceNoIn, int numSuppliersIn) {
		instanceNo = instanceNoIn;
		numSuppliers_tmp = numSuppliersIn;
		cost_tmp = new int[numSuppliersIn];
		prob_tmp = new double[numSuppliersIn];
	}

	~ProblemInstance() {
		delete[] cost_tmp;
		delete[] prob_tmp;
	}
	void AssignExecutionCost(int ecIn) { ec_tmp = ecIn; }
	void AssignCost(int t, int c_tmp) { cost_tmp[t] = c_tmp; }
	void AssignCostArray(int* costArrayIn) { cost_tmp = costArrayIn; }
	void AssignProb(int t, double p_tmp) { prob_tmp[t] = p_tmp; }
	void AssignProbArray(double* probArrayIn) { prob_tmp = probArrayIn; }
	void Assignh(int hIn) { hCostInstance = hIn; }

};
#endif