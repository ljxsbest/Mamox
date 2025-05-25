// HEADER ONLY CLASS
#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "utilities.h"

struct baseReaction
{
	baseReaction()
	{

	}
	baseReaction(std::vector<std::string> reactants_,
		std::vector<std::string> products_, bool isDup, bool isRev)
	{
		reactants = reactants_;
		products = products_;
		isDuplicate = isDup;
		isReversible = isRev;
	}

	
	void print()
	{
		for (int i = 0; i < reactants.size(); i++)
		{
			std::cout << reactants[i];
			if (i != reactants.size() - 1)
				std::cout << " + ";
		}
		std::cout << " => ";
		for (int i = 0; i < products.size(); i++)
		{
			std::cout << products[i];
			if (i != products.size() - 1)
				std::cout << " + ";
		}
		if (isDuplicate == true)
			std::cout << "            DUP" << std::endl;
		else
			std::cout << std::endl;
	}

	std::vector<std::string> reactants = {};
	std::vector<std::string> products = {};
	bool isDuplicate = false;
	bool isReversible = false;
};

class CKMechReader
{
public:
	// constructor:
	CKMechReader(std::string mechPath, std::string thermoPath, std::string glossaryPath);


	std::vector<std::string> mechSpecies;		// list of all the species in the kinetic file
	std::vector<std::string> thermoSpecies;		// list of all the species in the thermo file
	std::vector<std::string> glossarySpecies;	// list of all the species in the glosasry file
	std::vector<std::string> glossaryInChIs;	// list of all the InChIs with same order of 
	                                            // glossarySpecies
	std::stringstream speciesText;	    // stringstream containing just the species sector of the kinetic file
	std::stringstream reactionsText;	// stringstream containing just the reactions sector of the kinetic file
	std::stringstream thermoText;		// stringstream containin just the tehrmo sector of the thermo file
	std::vector<baseReaction> reacs;	// list of all the reactions

private:
	std::ifstream mech;					// file containing the kinetic mechanism
	std::ifstream thermo;				// file containing the thermodynamic data
	std::ifstream glossary;				// file containing the glossary
	bool isReaction(std::string line, baseReaction* reac);
};

inline CKMechReader::CKMechReader(std::string mechPath, std::string thermoPath, std::string glossaryPath)
	:mech(mechPath, std::ifstream::in)
	, thermo(thermoPath, std::ifstream::in)
	, glossary(glossaryPath, std::ifstream::in)
{
	// ################## READ GLOSSARY SPECIES #####################################
	{
		std::string line;
		std::string word;
		while (!glossary.eof())
		{
			std::getline(glossary, line);
			if (line == "")
				continue;
			std::stringstream lineSS(line);
			std::getline(lineSS, word, ' ');
			if (word == "")
				continue;
			glossarySpecies.push_back(word);
			while (!lineSS.eof())
			{
				std::getline(lineSS, word, ' ');
				if (word.size() > 0)
					break;
			}
			glossaryInChIs.push_back(word);
		}
		glossary.clear();
		glossary.seekg(0);
	}

	// ################## READ THERMO SPECIES #######################################
	{
		std::vector<std::string> glossarySpecCopy = glossarySpecies;
		std::string line;
		std::stringstream noCommThermo;
		while (!thermo.eof())
		{
			std::getline(thermo, line);
			std::stringstream lineSS(line);
			std::getline(lineSS, line, '!');
			noCommThermo << line << std::endl;
		}
		//std::cout << noCommThermo.str();
		while (!noCommThermo.eof())
		{
			std::getline(noCommThermo, line);
			if (line[79] == '1')
			{
				std::string word;
				std::stringstream lineSS(line);
				std::getline(lineSS, word, ' ');
				UTL::addUniqueString(&thermoSpecies, word);
				if (UTL::addUniqueString(&glossarySpecCopy, word) == -1)
					UTL::error("Species " + word + " is not present in glossary!");
			}
		}
		thermo.clear();
		thermo.seekg(0);
		bool startCopying = false;
		while (!thermo.eof())
		{
			std::getline(thermo, line);
			if (UTL::removeSpaces(line) == "END" || 
				UTL::removeSpaces(line) == "end")
				break;
			if (startCopying == true)
				thermoText << line << std::endl;
			if (UTL::removeSpaces(line).substr(0, 6) == "THERMO" ||
				UTL::removeSpaces(line).substr(0, 6) == "thermo")
			{
				std::getline(thermo, line);
				startCopying = true;
			}
		}
		thermo.clear();
		thermo.seekg(0);
	}

	// ################## READ MECH SPECIES #######################################
	{
		std::vector<std::string> thermoSpecCopy = thermoSpecies;
		std::string line;
		std::stringstream noCommMech;
		while (!mech.eof())
		{
			std::getline(mech, line);
			std::stringstream lineSS(line);
			std::getline(lineSS, line, '!');
			noCommMech << line << std::endl;
		}
		//std::cout << noCommMech.str();
		noCommMech.clear();
		noCommMech.seekg(0);
		while (!noCommMech.eof()) // go to start of species definition
		{
			std::getline(noCommMech, line);
			if (line.substr(0, 7) == "SPECIES" || line.substr(0, 7) == "species")
				break;
		}
		while (!noCommMech.eof())
		{
			std::getline(noCommMech, line);
			std::stringstream lineSS(line);
			std::string word;
			bool speciesEnd = false;
			while (!lineSS.eof())
			{
				std::getline(lineSS, word, ' ');
				if (word.substr(0, 3) == "END" || word.substr(0, 3) == "end")
				{
					speciesEnd = true;
					break;
				}
				if (word == "" || word == " ")
					continue;
				int indx = UTL::addUniqueString(&mechSpecies, word);
				if (indx != -1)
					UTL::warning("In reading species from base mech, duplicate species has been found");
				indx = UTL::addUniqueString(&thermoSpecCopy, word);
				if (indx == -1)
					UTL::error("Species " + word + " is not present in base thermo file!");

			}
			if (speciesEnd)
				break;
		}
		mech.clear();
		mech.seekg(0);
		bool startCopying = false;
		bool speciesSectionStarted = false;
		while (!mech.eof())
		{
			std::getline(mech, line);
			if ((UTL::removeSpaces(line) == "END" ||
				UTL::removeSpaces(line) == "end") && speciesSectionStarted)
				break;
			if (startCopying == true)
				speciesText << line << std::endl;
			if (UTL::removeSpaces(line).substr(0, 7) == "SPECIES"
				|| UTL::removeSpaces(line).substr(0, 7) == "species")
			{
				speciesSectionStarted = true;
				startCopying = true;
			}
		}
		mech.clear();
		mech.seekg(0);
	}
	
	//std::cout <<"Glossary InChIs  size: " << glossaryInChIs.size()<< std::endl;

	// ##################### READ REACTIONS #######################################
	{
		{
			baseReaction tempReac;
			std::string line;
			while (std::getline(mech, line))		// find the beginning of the species list
			{
				if (line.substr(0, 9) == "REACTIONS" 
					|| line.substr(0, 9) == "reactions")
					break;
			}
			while (std::getline(mech, line))
			{
				if (line.substr(0, 3) == "END"
					|| line.substr(0, 3) == "end")	// break if it finds end keyword
					break;
				if (isReaction(line, &tempReac))
					reacs.push_back(tempReac);

				if (line.substr(0, 3) == "DUP")
				{
					//std::cout << reactionsList.size() << std::endl;
					reacs[reacs.size() - 1].isDuplicate = true;
				}
			}
		}
		mech.clear();
		mech.seekg(0);
		bool startCopying = false;
		bool reactionSectionStarted = false;
		std::string line;
		while (!mech.eof())
		{
			std::getline(mech, line);
			if ((UTL::removeSpaces(line) == "END" || 
				UTL::removeSpaces(line) == "end") && reactionSectionStarted)
				break;
			if (startCopying == true)
				reactionsText << line << std::endl;
			if (UTL::removeSpaces(line).substr(0, 9) == "REACTIONS"
				|| UTL::removeSpaces(line).substr(0, 9) == "reactions")
			{
				reactionSectionStarted = true;
				startCopying = true;
			}
		}
		mech.clear();
		mech.seekg(0);
	}
	std::cout << " - " << glossarySpecies.size() << " glossary species imported" << std::endl;
	std::cout << " - " << mechSpecies.size() << " mechanism species imported" << std::endl;
	std::cout << " - " << thermoSpecies.size() << " species' thermodynamic data imported" << std::endl;
	std::cout << " - " << reacs.size() << " reactions imported" << std::endl;
}

inline bool CKMechReader::isReaction(std::string line, baseReaction* reac)
{
	if (line.length() == 0)
		return false;
	if (line.substr(0, 3) == "DUP")
		return false;
	if (line.substr(0, 4) == "PLOG")
		return false;
	if (line[0] == '!')
		return false;

	std::stringstream strStream(line);
	//remove comments if there are
	std::getline(strStream, line, '!');
	// check if there is the equal sign
	bool foundEqual = false;
	for (int i = 0; i < line.length(); i++)
	{
		if (line[i] == '=')
		{
			foundEqual = true;
			break;
		}
	}
	if (foundEqual == false)
		return false;
	// remove the kinetic parameters
	strStream = std::stringstream(line);
	strStream.seekp(0);
	std::vector<std::string> words;
	while (!strStream.eof())
	{
		std::string word;
		std::getline(strStream, word, ' ');
		if (word != "" && word != " ")
			words.push_back(word);
	}
	line = "";
	if (words.size() < 4)
		return false;
	for (int i = 0; i < words.size() - 3; i++)
		line = line + words[i];

	std::string reactantsString;
	std::string productsString;
	strStream = std::stringstream(line);

	std::getline(strStream, reactantsString, '=');
	std::getline(strStream, productsString, '=');

	bool rev = false;
	if (reactantsString[reactantsString.length() - 1] == '<')
	{
		reactantsString.pop_back();
		rev = true;
	}
	if (productsString[0] == '>')
		productsString.erase(0, 1);
	else
		rev = true;

	std::vector<std::string> reactants;
	std::vector<std::string> products;

	std::stringstream reactantsStrStream(reactantsString);
	std::stringstream productsStrStream(productsString);

	while (!reactantsStrStream.eof())
	{
		std::string word;
		std::getline(reactantsStrStream, word, '+');
		reactants.push_back(word);
	}
	while (!productsStrStream.eof())
	{
		std::string word;
		std::getline(productsStrStream, word, '+');
		products.push_back(word);
	}

	baseReaction r(reactants, products, false, rev);
	*reac = r;

	return true;
}