#pragma once
#include <vector>
#include <string>
#include <iostream>

struct Rate
{
	Rate(std::string ID_, std::vector<std::string> parameters_, double A_, double n_, double E_)
	{
		//ID = ID_;
		//parameters.resize(parameters_.size());
		//for (int i = 0; i < parameters.size(); i++)
		//	parameters[i] = parameters_[i];
		//A = A_;
		//n = n_;
		//E = E_;
		Rate(ID_, parameters_, A_, n_, E_, "");
	}

	Rate(std::string ID_, std::vector<std::string> parameters_, double A_, double n_, double E_, std::string ref)
	{
		ID = ID_;
		parameters.resize(parameters_.size());
		for (int i = 0; i < parameters.size(); i++)
			parameters[i] = parameters_[i];
		A = A_;
		n = n_;
		E = E_;
		reference = ref;
	}
	std::string ID;
	std::vector<std::string> parameters;
	double A;
	double n;
	double E;
	std::string reference = "";
};


class RateRule
{
public:
	RateRule(int numDim, bool isSymm);

	std::vector<double> returnRates(std::vector<std::string> param);
	std::vector<double> returnRates(std::vector<std::string> param, std::string *ident);
	void addRate(std::string ID, std::vector<std::string> param, double A, double n, double E);
	void addRate(std::string ID, std::vector<std::string> param, double A, double n, double E, std::string ref);
	int dimensionality() { return numDimensions; };
private:
	std::vector<Rate> rates;
	bool isSymmetric = false;
	int numDimensions = 0;
	bool zeroDim = false;
};

