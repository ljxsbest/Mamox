#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "molec.h"
#include "kinox.h"
#include "reaction.h"
#include "CKMechReader.h"
#include "utilities.h"


#define CHEMKIN_MAX_CHAR_PER_ROW  250



class ChemkinOut
{
public:
	ChemkinOut(CKMechReader* baseMechanism);	// constructor

	std::string molToName(Molecola mol); // return the name of the molecule, quit the program if the molecule cannot be named
	std::string molToNameLump(Molecola mol);
	std::string checkIfNamed(Molecola* mol);
	void printLongNameSpeciesMessage();

private:
	std::string alkaneStructurePrefix(Molecola* mol);
	std::string nameCOBranch(Molecola branch);
	CKMechReader* baseMech;
	void removeH(std::string* name);
	std::string specialMolName(Molecola mol);
	std::vector<Molecola> longNamedSpecies;
	std::vector<std::string> longNamedSpeciesNames;
	std::vector<std::string> longNamedSpeciesOriginalNames;
};