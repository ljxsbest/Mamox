#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

inline std::string reverseString(std::string str)
{
	std::string revStr = str;
	for (int i = 0; i < str.size(); i++)
		revStr[revStr.size() - i - 1] = str[i];
	return revStr;
}

inline std::string removeSpaces(std::string str)
{
	str.erase(remove(str.begin(), str.end(), ' '), str.end());
	return str;
}

inline int parseReaction(std::string line, std::vector<std::string>* reactants, std::vector<std::string>* products)
{
	std::vector<std::string> unwantedSubstrings = { "/", "PLOG", "DUP", "LOW", "TROE", "(S)" };

	int exitFlag = 1;			// 0 no reaction found
								// 1 irreversible reaction found 
								// 2 reversible reaction found
								// 3 third body reaction found

	std::stringstream lineSS(line);
	std::getline(lineSS, line, '!');

	if (line.size() == 0)			// if it was all a comment skip
		return 0;


	bool unwantedPresent = false;
	for (int i = 0; i < unwantedSubstrings.size(); i++)
	{
		if (line.find(unwantedSubstrings[i]) != std::string::npos)
		{
			unwantedPresent = true;
			break;
		}
	}
	if (unwantedPresent)
		return 0;

	bool isReversible = false;
	std::string arrowSign;
	if (line.find("<=>") != std::string::npos)
	{
		isReversible = true;
		arrowSign = "<=>";
		exitFlag = 2;
	}
	else if (line.find("=>") != std::string::npos)
	{
		arrowSign = "=>";
		exitFlag = 1;
	}
	else
	{
		//std::cout << "WARNING: reaction doesn't have <=> nor =>!" << std::endl;
		return 0;
	}

	if (line.find("(+M)") != std::string::npos)
		exitFlag = 3;

	std::stringstream reactantsSide(line.substr(0, line.find(arrowSign)));

	std::string word;
	
	while (!reactantsSide.eof())
	{
		std::getline(reactantsSide, word, '+');
		word = removeSpaces(word);	// remove spaces
		reactants->push_back(word);
	}


	std::stringstream productsSide(line.substr(line.find(arrowSign) + arrowSign.size(), line.size() - (line.find(arrowSign) + arrowSign.size())));
	//reverse the productsSide
	std::string revProductsSide = reverseString(productsSide.str());

	std::stringstream revProdSideSS(revProductsSide);
	// remove the first 3 numbers
	for (int i = 0; i < 3; i++)
	{
		bool wordFound = false;
		while (wordFound == false)
		{
			std::getline(revProdSideSS, word, ' ');
			if (word.size() > 0)
				wordFound = true;
		}
		word = reverseString(word);
		word = removeSpaces(word);

		try
		{
			std::stod(word);
		}
		catch (const std::invalid_argument& ia)
		{
			std::cout << "WARNING: the reaction doesn't end with 3 numbers!" << std::endl;
			return 0;
		}
	}

	{
		std::stringstream temp;
		temp << revProdSideSS.rdbuf();
		productsSide.str(reverseString(removeSpaces(temp.str())));
	}

	while (!productsSide.eof())
	{
		std::getline(productsSide, word, '+');
		word = removeSpaces(word);	// remove spaces
		products->push_back(word);
	}
	return exitFlag;
}

inline int parseReaction(std::string line, std::vector<std::vector< std::string >> *reaction)
{
	std::vector<std::string> reactants;
	std::vector<std::string> products;
	int flag = parseReaction(line, &reactants, &products);
	reaction->push_back(reactants);
	reaction->push_back(products);
	return flag;
}

inline int getNextReaction(std::ifstream* kinFile, std::vector<std::vector<std::string>>* reac, std::string* reacLine)
{
	std::string line;
	while (!kinFile->eof())
	{
		std::getline(*kinFile, line);
		std::vector<std::string> reactants;
		std::vector<std::string> products;
		int flag = parseReaction(line, &reactants, &products);
		if (flag == 0)
			continue;
		*reacLine = line;
		reac->push_back(reactants);
		reac->push_back(products);
		return flag;
	}
	return 0;
}

inline void parseKinFile(std::ifstream* inFile,
	std::vector<std::vector<std::vector<std::string>>>* irrReacList,
	std::vector<std::vector<std::vector<std::string>>>* revReacList)
{
	irrReacList->resize(0);
	revReacList->resize(0);
	while (!(*inFile).eof())
	{
		std::vector<std::vector<string>> reaction;
		std::string reacLine;
		int flag = getNextReaction(inFile, &reaction, &reacLine);
		if (flag == 1)
			irrReacList->push_back(reaction);
		if (flag == 2)
			revReacList->push_back(reaction);
	}
	inFile->clear();
	inFile->seekg(0);
}

inline void parseKinFile(std::string inFilePath,
	std::vector<std::vector<std::vector<std::string>>>* irrReacList,
	std::vector<std::vector<std::vector<std::string>>>* revReacList)
{
	std::ifstream inFile;
	inFile.open(inFilePath);
	parseKinFile(&inFile, irrReacList, revReacList);
}

inline bool isDuplicate(std::vector<std::vector<std::string>>* reaction,
	std::vector<std::vector<std::vector<std::string>>>* irrReacList,
	std::vector<std::vector<std::vector<std::string>>>* revReacList)
{
	bool isDup = false;

	// check with irreversible reactions
	for (int i = 0; i < irrReacList->size(); i++)
	{
		//std::cout << (*reaction)[0].size() << std::endl;
		if ((*reaction)[0].size() != (*irrReacList)[i][0].size())		// if they have different number of reactants they are not the same
			continue;
		if ((*reaction)[1].size() != (*irrReacList)[i][1].size())		// if they have different number of products they are not the same
			continue;

		bool notDup = false;
		// check the reactants
		for (int j = 0; j < (*reaction)[0].size(); j++)
		{
			bool thereIsAMatch = false;
			for (int k = 0; k < (*irrReacList)[i][0].size(); k++)
				if ((*reaction)[0][j] == (*irrReacList)[i][0][k])
					thereIsAMatch = true;
			if (thereIsAMatch == false)
			{
				notDup = true;
				break;
			}
		}
		if (notDup)
			continue;
		// check the products
		for (int j = 0; j < (*reaction)[1].size(); j++)
		{
			bool thereIsAMatch = false;
			for (int k = 0; k < (*irrReacList)[i][1].size(); k++)
				if ((*reaction)[1][j] == (*irrReacList)[i][1][k])
					thereIsAMatch = true;
			if (thereIsAMatch == false)
			{
				notDup = true;
				break;
			}
		}
		if (notDup)
			continue;
		// if all the criteria are met then they are duplciates
		return true;
	}

	// check with direct reversible reactions
	for (int i = 0; i < revReacList->size(); i++)
	{
		if ((*reaction)[0].size() != (*revReacList)[i][0].size())		// if they have different number of reactants they are not the same
			continue;
		if ((*reaction)[1].size() != (*revReacList)[i][1].size())		// if they have different number of products they are not the same
			continue;

		bool notDup = false;
		// check the reactants
		for (int j = 0; j < (*reaction)[0].size(); j++)
		{
			bool thereIsAMatch = false;
			for (int k = 0; k < (*revReacList)[i][0].size(); k++)
				if ((*reaction)[0][j] == (*revReacList)[i][0][k])
					thereIsAMatch = true;
			if (thereIsAMatch == false)
			{
				notDup = true;
				break;
			}
		}
		if (notDup)
			continue;
		// check the products
		for (int j = 0; j < (*reaction)[1].size(); j++)
		{
			bool thereIsAMatch = false;
			for (int k = 0; k < (*revReacList)[i][1].size(); k++)
				if ((*reaction)[1][j] == (*revReacList)[i][1][k])
					thereIsAMatch = true;
			if (thereIsAMatch == false)
			{
				notDup = true;
				break;
			}
		}
		if (notDup)
			continue;
		// if all the criteria are met then they are duplciates
		return true;
	}

	// check with reverse reversible reactions
	for (int i = 0; i < revReacList->size(); i++)
	{
		if ((*reaction)[0].size() != (*revReacList)[i][1].size())		// if they have different number of reactants they are not the same
			continue;
		if ((*reaction)[1].size() != (*revReacList)[i][0].size())		// if they have different number of products they are not the same
			continue;

		bool notDup = false;
		// check the reactants
		for (int j = 0; j < (*reaction)[0].size(); j++)
		{
			bool thereIsAMatch = false;
			for (int k = 0; k < (*revReacList)[i][1].size(); k++)
				if ((*reaction)[0][j] == (*revReacList)[i][1][k])
					thereIsAMatch = true;
			if (thereIsAMatch == false)
			{
				notDup = true;
				break;
			}
		}
		if (notDup)
			continue;
		// check the products
		for (int j = 0; j < (*reaction)[1].size(); j++)
		{
			bool thereIsAMatch = false;
			for (int k = 0; k < (*revReacList)[i][0].size(); k++)
				if ((*reaction)[1][j] == (*revReacList)[i][0][k])
					thereIsAMatch = true;
			if (thereIsAMatch == false)
			{
				notDup = true;
				break;
			}
		}
		if (notDup)
			continue;
		// if all the criteria are met then they are duplciates
		return true;
	}

	// if none of the criteria is met then there are no duplicates
	return false;
}

inline void parseSpecies(std::ifstream* inFile, std::vector<std::string>* species)
{
	inFile->clear();
	inFile->seekg(0);
	
	species->resize(0);

	std::string line;

	while (!(*inFile).eof())
	{
		std::getline(*inFile, line);
		if (line.substr(0, 8) == "SPECIES" || line.substr(0, 8) == "species")
			break;
	}

	while (!(*inFile).eof())
	{
		std::getline(*inFile, line);
	
		if (line.substr(0, 3) == "END" || line.substr(0, 3) == "end")
			break;
		
		if (line.size() == 0)
			continue;
		std::stringstream lineSS(line);
		std::getline(lineSS, line, '!'); // remove comments
		if (line.size() == 0)
			continue;
		lineSS.str(line);
		std::string word;
		while (!lineSS.eof())
		{
			std::getline(lineSS, word);
			if (word.size() != 0)
				species->push_back(word);
		}
	}

	inFile->clear();
	inFile->seekg(0);
}

inline void parseSpecies(std::string inFilePath, std::vector<std::string>* species)
{
	std::ifstream inFile;
	inFile.open(inFilePath);
	parseSpecies(&inFile, species);
	inFile.close();
}