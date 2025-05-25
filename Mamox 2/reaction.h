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
	void print(std::ostream& stream);
	std::vector<Molecola*> reactantList();
	std::vector<Molecola*> productList();

	std::vector<Molecola>& Reactants() { return reactants; }; // here, add & to make sure the assignment is working. 
	std::vector<Molecola>& Products() { return products; };


	void updateReactant(int index, Molecola* mol);
	void updateProduct(int index, Molecola* mol);
	double A() { return kineticParameters[0]; };
	double n() { return kineticParameters[1]; };
	double E() { return kineticParameters[2]; };
	std::string reactionLabel() { return label; };
	void setA(double A_) { kineticParameters[0] = A_; };
	friend std::ostream& operator<<(std::ostream& stream, Reaction& reaction);
	int operator==(Reaction b);
	void setMultiplePathMultiplier(int mult)
	{
		comment.multiPathMultiplier = mult;
	}
	std::string printComment()
	{
		return comment.comment();
	}
	bool isDuplicate() { return duplicate; };
	void setDuplicate() { duplicate = true; };
	bool weakEquality(Reaction other);		// return true if other has same reactants 
											// and products (don't care about rate)
	Molecola parentFuel();	
	Molecola parentFuelOH();
	Molecola parentFuel_complete();
	Molecola reac_size();
	Molecola prod_size();
	//Molecola reac_clean();

	//std::vector<Molecola*> Reactants() { return reactants; };
	//std::vector<Molecola*> Products() { return products; };


	void setReversible(bool value) { reversible = value; };
	bool isReversible() { return reversible; };

private:
	std::vector<Molecola> reactants;
	std::vector<Molecola> products;
	double kineticParameters[3];			// [A, n, E]
	std::string label;
	reactionComment comment;
	bool duplicate = false;
	bool reversible = false;
};