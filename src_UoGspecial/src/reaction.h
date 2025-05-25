#pragma once

#include "molec.h"
#include <vector>
#include "ReactionComment.h"

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

	Reaction(Molecola              reac, Molecola              prod, double kinPar[3], std::string name_, reactionComment comment);
	Reaction(Molecola              reac, std::vector<Molecola> prod, double kinPar[3], std::string name_, reactionComment comment);
	Reaction(std::vector<Molecola> reac, std::vector<Molecola> prod, double kinPar[3], std::string name_, reactionComment comment);
	Reaction(std::vector<Molecola> reac, Molecola              prod, double kinPar[3], std::string name_, reactionComment comment);


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
	bool weakEquality(Reaction b);  // return true if the two reactions are equal apart from the pre-exponential factor
	//std::string comment;
	void setMultiplePathMultiplier(int mult)
	{
		comment.multiPathMultiplier = mult;
	}
	std::string printComment()
	{
		return comment.comment();
	}
	//bool isPLOG() { return isPLOG; }
	//void makePLOG(); 
private:
	std::vector<Molecola> reactants;
	std::vector<Molecola> products;
	double kineticParameters[3];			// [A, n, E]
	std::string label;
	reactionComment comment;
	//bool hasPLOG = False;
	//std::vector<float> plogValues;
};