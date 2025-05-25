#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "molec.h"

class XML
{
public:

	// constructor
	//		if mode is 'n' a new file is made;
	//		if mode is 'r' the file is opened in read mode
	XML(std::string fileName_, char mode);

	void addMolecule(Molecola* mol, std::string name);
	void getMolecules(std::vector<Molecola>* vec, std::vector<std::string>* names);
	void close();

private:
	std::string fileName;
	std::fstream file;
	char mode;
	std::string getNextKeyword();
	void goToKeyword(std::string keyword);

};

