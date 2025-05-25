#pragma once

#include "molec.h"
#include <vector>

struct reactionComment
{
	std::string rateRule;
	int HMultiplier = 1;
	int mutliPathMultiplier = 1;

	std::string comment()
	{
		std::string comm;
		comm.append(rateRule);
		if (HMultiplier > 1)
		{
			comm.append(", X");
			comm.append(std::to_string(HMultiplier));
			comm.append(" H multiplier");
		}
		if (mutliPathMultiplier > 1)
		{
			comm.append(", X");
			comm.append(std::to_string(mutliPathMultiplier));
			comm.append(" multiple path multiplier");
		}
	}
};

class Reaction
{
public:
	//Reaction(Molecola              reac, Molecola              prod, double kinPar[3]);
	//Reaction(Molecola              reac, std::vector<Molecola> prod, double kinPar[3]);
	//Reaction(std::vector<Molecola> reac, std::vector<Molecola> prod, double kinPar[3]);
	//Reaction(std::vector<Molecola> reac, Molecola              prod, double kinPar[3]);
	
	Reaction(Molecola              reac, Molecola              prod, double kinPar[3], std::string name_);
	Reaction(Molecola              reac, std::vector<Molecola> prod, double kinPar[3], std::string name_);
	Reaction(std::vector<Molecola> reac, std::vector<Molecola> prod, double kinPar[3], std::string name_);
	Reaction(std::vector<Molecola> reac, Molecola              prod, double kinPar[3], std::string name_);


	double rateConstant(double Temperature);				// return the rate constant at the temperature
	void print(ostream& stream);
	std::vector<Molecola*> reactantList();
	std::vector<Molecola*> productList();
	void updateReactant(int index, Molecola* mol);
	void updateProduct(int index, Molecola* mol);
	double A() { return kineticParameters[0]; };
	double n() { return kineticParameters[1]; };
	double E() { return kineticParameters[2]; };
	std::string reactionLabel() { return label; };
	void setA(double A_) { kineticParameters[0] = A_; };
	friend ostream& operator<<(ostream& stream, Reaction& reaction);
	int operator==(Reaction b);

private:
	std::vector<Molecola> reactants;
	std::vector<Molecola> products;
	double kineticParameters[3];			// [A, n, E]
	std::string label;
};