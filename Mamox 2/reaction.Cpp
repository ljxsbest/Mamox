#include "reaction.h"

Reaction::Reaction(Molecola reac, Molecola prod, double kinPar[3], std::string label_)
{
	reactants.push_back(reac);
	products.push_back(prod);
	kineticParameters[0] = kinPar[0];
	kineticParameters[1] = kinPar[1];
	kineticParameters[2] = kinPar[2];
	label = label_;
}

Reaction::Reaction(Molecola reac, std::vector<Molecola> prod, double kinPar[3], std::string label_)
{
	reactants.push_back(reac);
	for(int i = 0; i<prod.size(); i++)
		products.push_back(prod[i]);

	kineticParameters[0] = kinPar[0];
	kineticParameters[1] = kinPar[1];
	kineticParameters[2] = kinPar[2];
	label = label_;
}

Reaction::Reaction(std::vector<Molecola> reac, std::vector<Molecola> prod, double kinPar[3], std::string label_)
{
	for (int i = 0; i < reac.size(); i++)
		reactants.push_back(reac[i]);
	for (int i = 0; i < prod.size(); i++)
		products.push_back(prod[i]);

	kineticParameters[0] = kinPar[0];
	kineticParameters[1] = kinPar[1];
	kineticParameters[2] = kinPar[2];
	label = label_;
}

Reaction::Reaction(std::vector<Molecola> reac, Molecola prod, double kinPar[3], std::string label_)
{
	for (int i = 0; i < reac.size(); i++)
		reactants.push_back(reac[i]);
	
	products.push_back(prod);

	kineticParameters[0] = kinPar[0];
	kineticParameters[1] = kinPar[1];
	kineticParameters[2] = kinPar[2];
	label = label_;
}

Reaction::Reaction(Molecola reac, Molecola prod, double kinPar[3], std::string label_, reactionComment comm)
{
	reactants.push_back(reac);
	products.push_back(prod);
	kineticParameters[0] = kinPar[0];
	kineticParameters[1] = kinPar[1];
	kineticParameters[2] = kinPar[2];
	label = label_;
	comment = comm;
}
Reaction::Reaction(Molecola reac, std::vector<Molecola> prod, double kinPar[3], std::string label_, reactionComment comm)
{
	reactants.push_back(reac);
	for (int i = 0; i < prod.size(); i++)
		products.push_back(prod[i]);

	kineticParameters[0] = kinPar[0];
	kineticParameters[1] = kinPar[1];
	kineticParameters[2] = kinPar[2];
	label = label_;
	comment = comm;
}
Reaction::Reaction(std::vector<Molecola> reac, std::vector<Molecola> prod, double kinPar[3], std::string label_, reactionComment comm)
{
	for (int i = 0; i < reac.size(); i++)
		reactants.push_back(reac[i]);
	for (int i = 0; i < prod.size(); i++)
		products.push_back(prod[i]);

	kineticParameters[0] = kinPar[0];
	kineticParameters[1] = kinPar[1];
	kineticParameters[2] = kinPar[2];
	label = label_;
	comment = comm;
}
Reaction::Reaction(std::vector<Molecola> reac, Molecola prod, double kinPar[3], std::string label_, reactionComment comm)
{
	for (int i = 0; i < reac.size(); i++)
		reactants.push_back(reac[i]);

	products.push_back(prod);

	kineticParameters[0] = kinPar[0];
	kineticParameters[1] = kinPar[1];
	kineticParameters[2] = kinPar[2];
	label = label_;
	comment = comm;
}

void Reaction::print(std::ostream& stream)
{
	for (int i = 0; i < reactants.size(); i++)
	{
		stream << reactants[i];
		if (i != reactants.size() - 1) stream << " + ";
	}
	if(reversible == false)
		stream << " => ";
	else
		stream << " = ";
	for (int i = 0; i < products.size(); i++)
	{
		stream << products[i];
		if (i != products.size() - 1) stream << " + ";
	}
	stream << std::endl;
}

double Reaction::rateConstant(double Temperature)
{
	return kineticParameters[0] * pow(Temperature, kineticParameters[1]) * exp(-kineticParameters[2] / 1.987 / Temperature);
}

std::vector<Molecola*> Reaction::reactantList()
{
	std::vector<Molecola*> vec;
	for (int i = 0; i < reactants.size(); i++)
		vec.push_back(&(reactants[i]));

	return vec;
}

std::vector<Molecola*> Reaction::productList()
{
	std::vector<Molecola*> vec;
	for (int i = 0; i < products.size(); i++)
		vec.push_back(&(products[i]));

	return vec;
}

void Reaction::updateReactant(int index, Molecola* mol)
{
	if (*mol == reactants[index]) reactants[index] = *mol;
	else
	{
		std::cerr << "ERROR: trying to replace " << reactants[index] << " with " << *mol << " they are not equivalent!";
		exit(0);
	}
}

void Reaction::updateProduct(int index, Molecola* mol)
{
	if (*mol == products[index]) products[index] = *mol;
	else
	{
		std::cerr << "ERROR: trying to replace " << products[index] << " with " << *mol << " they are not equivalent!";
		exit(0);
	}
}

std::ostream& operator<<(std::ostream& stream, Reaction& reaction)
{
	std::vector<Molecola*> reactants = reaction.reactantList();
	for (int i = 0; i < reactants.size(); i++)
	{
		stream << *(reactants[i]);
		if (i != reactants.size() - 1)
			stream << " + ";
	}
	stream << " => ";
	std::vector<Molecola*> products = reaction.productList();
	for (int i = 0; i < products.size(); i++)
	{
		stream << *(products[i]);
		if (i != products.size() - 1)
			stream << " + ";
	}
	stream << "         " << reaction.A() << "   " << reaction.n() << "   " << reaction.E();
	return stream;

}

bool areMoleculeListEquivalent(std::vector<Molecola*> vec1, std::vector<Molecola*> vec2)
{
	if (vec1.size() != vec2.size())
		return false;

	std::vector<bool> alreadyMatched(vec1.size(), false);
	for (int i = 0; i < vec1.size(); i++)
	{
		for (int j = 0; j < vec2.size(); j++)
		{
			if (*(vec1[i]) == *(vec2[j]) && alreadyMatched[j] == false)
			{
				alreadyMatched[j] = true;
			}
		}
	}

	for (int i = 0; i < alreadyMatched.size(); i++)
		if (alreadyMatched[i] == false)
			return false;

	return true;
}

int Reaction::operator==(Reaction b)
{
	if (A() != b.A() || n() != b.n() || E() != b.E())	// if the kinetic parameters are different they are not equal so continue
		return 0;
	
	if (reactantList().size() != b.reactantList().size())		// if they have not the same number of reactants they are not equal so continue
		return 0;

	if (productList().size() != b.productList().size())			// if they have not the same number of products they are not equal so continue
		return 0;

	// if they have the same kinetics parameters, and the same number of products and reactants then check if the products and reactants are equal
	// TO DO N.B. this should be improved because it works only if the products and reactants are in the same order in the two reactions 
	if (areMoleculeListEquivalent(reactantList(), b.reactantList()) == false)
		return 0;
	if (areMoleculeListEquivalent(productList(), b.productList()) == false)
		return 0;

	// if all the criteria are met then the two reactions are equal
	return 1;
}

bool Reaction::weakEquality(Reaction other)
{
	if (isReversible() and other.isReversible())
	{
		if(areMoleculeListEquivalent(reactantList(), other.productList()) && areMoleculeListEquivalent(productList(), other.reactantList()))
			return true;
	}
	if (areMoleculeListEquivalent(reactantList(), other.reactantList()) && areMoleculeListEquivalent(productList(), other.productList()))
		return true;
	return false;
	
	//if (areMoleculeListEquivalent(reactantList(), other.reactantList()) == false)
	//	return false;
	//if (areMoleculeListEquivalent(productList(), other.productList()) == false)
	//	return false;
	//
	//
	//return true;
}

Molecola Reaction::parentFuel()
{
	if (reactants.size() == 1)
	{
		return reactants[0].parentFuel();
	}
	else if (reactants.size() == 2)
	{
		if (reactants[0].isSpecialMolecule() != 0)
			return reactants[1].parentFuel();
		else if(reactants[1].isSpecialMolecule() != 0)
			return reactants[0].parentFuel();
		else if (reactants[0].size() > reactants[1].size())
			return reactants[0].parentFuel();
		else if (reactants[0].size() < reactants[1].size())
			return reactants[1].parentFuel();
		else
		{
			UTL::error("In Reaction::parentFuel an unexpected condition occured (1)");
		}
	}
	else 
		UTL::error("In Reaction::parentFuel an unexpected condition occured (2)");
	return Molecola(1);
}

Molecola Reaction::parentFuelOH()
{
	if (reactants.size() == 1)
	{
		return reactants[0].parentFuel();
	}
	else if (reactants.size() == 2)
	{
		if (reactants[0].isSpecialMolecule() != 0)
			return reactants[1].parentFuelOH();
		else if (reactants[1].isSpecialMolecule() != 0)
			return reactants[0].parentFuelOH();
		else if (reactants[0].size() > reactants[1].size())
			return reactants[0].parentFuelOH();
		else if (reactants[0].size() < reactants[1].size())
			return reactants[1].parentFuelOH();
		else
		{
			UTL::error("In Reaction::parentFuel an unexpected condition occured (1)");
		}
	}
	else
		UTL::error("In Reaction::parentFuel an unexpected condition occured (2)");
	return Molecola(1);
}

Molecola Reaction::parentFuel_complete()
{
	if (reactants.size() == 1)
	{
		return reactants[0].parentFuel();
	}
	else if (reactants.size() == 2)
	{
		if (reactants[0].isSpecialMolecule() != 0)
			return reactants[1].parentFuel_complete();
		else if (reactants[1].isSpecialMolecule() != 0)
			return reactants[0].parentFuel_complete();
		else if (reactants[0].size() > reactants[1].size())
			return reactants[0].parentFuel_complete();
		else if (reactants[0].size() < reactants[1].size())
			return reactants[1].parentFuel_complete();
		else
		{
			UTL::error("In Reaction::parentFuel an unexpected condition occured (1)");
		}
	}
	else
		UTL::error("In Reaction::parentFuel an unexpected condition occured (2)");
	return Molecola(1);
}

Molecola Reaction::reac_size()
{
	if (reactants.size() == 1)
	{
		return reactants[0].size();
	}
	else if (reactants.size() == 2)
	{
		int a = reactants[0].size();
		int b = reactants[1].size();
		int c = std::max(a, b);
		return c;
	}
	else
	{
			UTL::error("In Reaction::reac_size an unexpected condition occured (1)");
	}
	return Molecola(1);
}

Molecola Reaction::prod_size()
{
	if (products.size() == 1)
	{
		return products[0].size();
	}
	else if (products.size() == 2)
	{
		int a = products[0].size();
		int b = products[1].size();
		int c = std::max(a, b);
		return c;
	}
	else
	{
		UTL::error("In Reaction::reac_size an unexpected condition occured (1)");
	}
	return Molecola(1);
}


