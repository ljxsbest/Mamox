#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "molec.h"
#include "kinox.h"
#include "reaction.h"
#include "definitions.h"

#define CHEMKIN_MAX_CHAR_PER_ROW  250



class ChemkinOut
{
public:
	
	ChemkinOut(std::string fileNameDetailed);	// constructor
	void wrireaDetailed(int nreaz, double A, double n, double E, Molecola m1, ...);
	void writeDetReacName(std::string comm);
	void writeHeadingDetailed(std::stringstream *sstr);
	void writeHeadingLumped(std::stringstream* sstr);
	void endlineDetailed();
	void printAllListsDetailed();
	void completeDetailed();
	void coutLumpedreactions() { std::cout << lumpedReactions.str() << std::endl; };
	void completeLumped();
	void wrireaDetailed(Reaction rea);
	void writeDetailedReactions(std::vector<Reaction>* vec);
	void wrireaLump(std::string label, std::vector<Molecola*> reactants, std::vector<double> stoicCoeff, std::vector<Molecola*> products, double A, double n, double E);
	void wrireaLump(std::string label, std::vector<Molecola*> reactants, std::vector<Molecola*> products, double A, double n, double E);
	void wrireaLump(std::string label, std::vector<Molecola*> reactants, std::vector<double> stoicCoeff, std::vector<Molecola*> products, vector<double> Press, vector<double> A, vector<double> n, vector<double> E, int indexP);
	void wrireaLump(std::string label, std::vector<Molecola*> reactants, std::vector<Molecola*> products, vector<double> Press, vector<double> A, vector<double> n, vector<double> E, int indexP);
	void wrireaLump(Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	void wrireaLump1(Molecola reag, HAbsRad rad, int numC, double A, double n, double E);
	void wrireaLump3(Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	void wrireaLump4(int numC, double A, double n, double E);     // R'numC' + O2 -> R'numC'OO
	void wrireaLump5(int numC, double A, double n, double E);     // R'numC'OO -> R'numC' + O2
	void wrireaLump6(int numC, double A, double n, double E);     // R'numC' + O2 -> OLE'numC' + HO2
	void wrireaLump7(int numC, double A, double n, double E);     // R'numC'OO -> Q'numC'OOH
	void wrireaLump8(int numC, double A, double n, double E);     // Q'numC'OOH -> R'numC'OO
	void wrireaLump9(int numC, double A, double n, double E);     // Q'numC'OOH -> ETER'numC' + OH
	void wrireaLump10(int numC, double A, double n, double E);     // Q'numC'OOH -> OLE'numC' + HO2
	void wrireaLump11(Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	void wrireaLump11b(Molecola reag, Molecola prod[], std::vector<double> stoicCoeff1, double A1, double n1, double E1,
						std::vector<double> stoicCoeff2, double A2, double n2, double E2);
	void wrireaLump12(int numC, double A, double n, double E);     // Q'numC'OOH + O2 -> OOQ'numC'OOH
	void wrireaLump13(int numC, double A, double n, double E);     // OOQ'numC'OOH -> Q'numC'OOH + O2
	void wrireaLump14(int numC, double A, double n, double E);     // OOQ'numC'OOH -> OQ'numC'OOH + OH
	void wrireaLump15(Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	std::string nameHAbsRad(HAbsRad rad);				// return a string containing the name of the radical rad
	std::string nameHAbsRadPlusH(HAbsRad rad);			// return a string containing the name of the radical rad when an H is added to it
	void printLumpName(std::string name);				// print the name of the reaction to lumpedReactions
	void setFuel(Molecola* fuel_) { fuel = molToName(*fuel_); };
	std::string molToName(Molecola mol);
	std::string molToNameLump(Molecola mol);
	void makeThermoOut(std::string path, std::vector<Molecola*>* specList);
	void makeThermoOutLumped(std::string path);
	std::string fuelName() { return fuel; };
	void addSpeciesToKinFile(Molecola mol);

	int sizeLumpedSpeciesList() { return lumpedMolecules.size(); }; // return the number of lumped species
	std::string ithLumpedSpeciesName(int i) { return lumpedNames[i]; }; // return the name of the i-th lumped species
	Molecola ithLumpedSpecies(int i) { return *(lumpedMolecules[i]); }; // return the i-th lumped species

private:
	int max_species = 12;
	ofstream outFileDetailed;
	ofstream outSpeciesMap;
	std::stringstream detailedReactions;
	std::stringstream detailedHeading;
	std::stringstream lumpedReactions;
	std::stringstream lumpedHeading;
	std::string fuel;
	std::vector<std::string> genericSpecies;
	std::vector<std::string> genericRadicals;
	std::vector<std::string> RList;
	std::vector<std::string> ROOList;
	std::vector<std::string> QOOHList;
	std::vector<std::string> OLEList;
	std::vector<std::string> COList;
	std::vector<std::string> EtherList;
	std::vector<std::string> OOQOOHList;
	std::vector<std::string> ROList;
	std::vector<std::string> KHPList;
	std::vector<std::string> ROOHList;
	std::vector<std::string> POOH2List;
	std::vector<std::string> oleOOHList;
	std::vector<std::string> etherOOHList;
	std::vector<std::string> oleRList;
	std::vector<std::string> oleCOList;
	std::vector<std::string> etherRList;
	std::vector<std::string> etherCOList;
	std::vector<std::string> linEtheROList;
	std::vector<std::string> alkROList;
	std::vector<std::string> lumpedSpeciesList;
	std::vector<std::vector<std::string>> detReacSpecies;
	std::vector<std::string> detReacConst;
	std::vector<std::string> detReacComments;
	std::vector<bool> isDetailedReactionDuplicate;
	std::vector<std::string> detReacName;
	std::vector<int> detReacNamePosition;
	// lists of named isomers used to find the number in the name that identifies the isomer
	std::vector<Molecola> fuelIsomList;					
	std::vector<Molecola> linEtheROIsomersList;
	std::vector<Molecola> RIsomList;			
	std::vector<Molecola> ROOIsomList;			
	std::vector<Molecola> QOOHIsomList;			
	std::vector<Molecola> OOQOOHIsomList;		
	std::vector<Molecola> oleIsomList;			
	std::vector<Molecola> COIsomList;			
	std::vector<Molecola> etherIsomList;		
	std::vector<Molecola> ROIsomList;			
	std::vector<Molecola> KHPIsomList;			
	std::vector<Molecola> ROOHIsomList;			
	std::vector<Molecola> POOH2IsomList;		
	std::vector<Molecola> oleOOHIsomList;		
	std::vector<Molecola> cycOOHIsomList;		
	std::vector<Molecola> oleRIsomList;			
	std::vector<Molecola> oleCOIsomList;		
	std::vector<Molecola> etherRIsomList;		
	std::vector<Molecola> etherCOIsomList;		
	std::vector<Molecola> linEtherCOIsomList;	
	std::vector<Molecola> alkROIsomList;	
	
	std::vector<std::string> lumpedNames;		// all unique names of lumped molecules
	std::vector<Molecola*> lumpedMolecules;		// molecules with the same order of lumpedNames

	std::vector<Molecola> namedSpecies;
	std::vector<string> namesList;

	std::vector<double> roundStoicCoeff(std::vector<double>* vect, int numOfDigit);		// return a vector with the element of vect rounded leaving numOfDifit of decimal digits
	bool addUniqueName(std::vector<std::string>* vect, std::string name);	// append name to vect only if name is not already present in vect
																			// return true if name was not present, otherwise return false
	std::string alkaneStructurePrefix(Molecola* mol);
	void addSpeciesToList(species type, std::string name);			// add the species name to the list of species at the beggining of th e file
	int Add(Molecola newMol, std::vector<Molecola>* vec);
	std::vector<int> searchDuplicates(std::vector<std::vector<std::string>>* vect, std::vector<std::string> reac );	// return the list of index of the elements in vect that are equal to reac
	int findPosition(std::vector<int>* vect, int num);		// if num is present in vect return its position, otherwise return -1
	int addUniqueLumpedSpecies(Molecola* mol);
	std::string checkIfNamed(Molecola* mol);				// check if the molecule is present in the list of named species, if yes return the name, otherwise return "not found"
};