#pragma once
#include <stdlib.h>
#include<iostream>
#include <string>
#include <vector>
#include "molec.h"
#include "chemkinOut.h"
#include <direct.h>
#include "Eigen/dense"
#include "nlopt.hpp"

struct Group
{
	std::string name;	// name of the group
	std::string units;	// unit of measure
	double Hf;			// contribution to the heat of formation
	double S;			// contribution to the entropy
	double Cp300;		// contribution to the specific heat at 300 K
	double Cp400;		// contribution to the specific heat at 400K
	double Cp500;		// contribution to the specific heat at 500K
	double Cp600;		// contribution to the specific heat at 600K
	double Cp800;		// contribution to the specific heat at 80K
	double Cp1000;		// contribution to the specific heat at 1000K
	double Cp1500;		// contribution to the specific heat at 1500K
};

class ThermoOut
{
public:
	ThermoOut(std::string groupsPath, std::string radicalCorrectionsPath, std::string knownValuesPath,
	ChemkinOut* chemout_);	
			// groupPath: the path of the file containing the group contributions
			// radicalCorrectionsPath: the path of the file containing the radicals corrections
			// knownValuiesPath: the path to the file of known values
	ThermoOut(std::string groupsPath, std::string radicalCorrectionsPath, ChemkinOut* chemout_);
		// groupPath: the path of the file containing the group contributions
		// radicalCorrectionsPath: the path of the file containing the radicals corrections

	std::string NASAOutputLumped(std::string name, std::vector<Molecola> mols, 
		std::vector<double> temperatures, std::vector<std::vector<double>> distr);
		// return the NASA formatted thermodynamic properties of the lumped molecule representing the
		// molecules in the vector mols distributed with the distributions in the vector distr

	std::string NASAOutput(Molecola* mol);		// return the NASA formatted thermodynamic parameters of the moelcules mol
	std::string NASAOutput(Molecola* mol, std::string name);
	double Hf(Molecola mol, double T);	// return the heat of formation at the temperature T in [cal/mol]
	double S(Molecola mol, double T);	// return the entropy at the temperature T in [cal/mol/K]

private:

	std::ifstream groupsFile;						// file containing the contributions of the groups
	std::ifstream correctionsFile;					// file containing the corrections of the radicals
	std::ifstream knownFile;						// file containing the thermodynamic values of the known species
	bool knownFileProvided;							// true if the knownFile is provided, otherwise false
	std::vector<std::string> knownMoleculesNames;	// list of the names of all the molecules we know the thermo data of
	std::vector<std::string> knownMoleculesData;	// the complete text of the thermo data for each species of knownMoleculesNames (same order)
	ChemkinOut* chemOut;							// the pointer to the ChemkinOut object used by the main script
	std::vector<Group> groupList;					// list of the groups
	std::vector<Group> correctionsList;				// list of the corrections

	int isInKnownList(Molecola* mol);				// return -1 if the molecule is not in the list of known values, otherwise return it's intdex in the vector of knownValues
	int isInKnownList(std::string molName);			// return -1 if the molecule is not in the list of known values, otherwise return it's intdex in the vector of knownValues

	void importGroups();
	void importCorrections();
	void addGroupToList(std::string groupName, std::vector<std::string>* groupVec, 
		std::vector<int>* freqVec);
	Group findGroup(std::string groupName, std::vector<Group>* groupVector);
	void fitNASAParam(Eigen::VectorXd temps, Eigen::VectorXd Cps, double Hf, double T_Hf, double S, double T_S, Eigen::VectorXd* params);
	void fitNASAParam(std::vector<double> temps, std::vector<double> Cps, double Hf, double T_Hf, double S, double T_S, std::vector<double>* params);
	std::string NASAParametersToString(std::string name, int numC, int numH, int numO, std::vector<double> Tlims, Eigen::VectorXd lowTParams, Eigen::VectorXd highTParams);
	std::string NASAParametersToString(std::string name, int numC, int numH, int numO, std::vector<double> Tlims, std::vector<double> lowTParams, std::vector<double> highTParams);
	double Hf_T(double T, Eigen::VectorXd params);	// return the heat of formation at the temperature T in [cal/mol]
	double S_T(double T, Eigen::VectorXd params);	// return the entropy at the temperature T in [cal/mol/K]
	double Cp_T(double T, Eigen::VectorXd params);	// return the specific heat at the temperature T in [cal/mol/K]
	double Hf_T(double T, std::vector<double> params);	// return the heat of formation at the temperature T in [cal/mol]
	double S_T(double T,  std::vector<double> params);	// return the entropy at the temperature T in [cal/mol/K]
	double Cp_T(double T, std::vector<double> params);	// return the specific heat at the temperature T in [cal/mol/K]	
	void knownMoleculeData(std::string molName, std::vector<double>* parLowT, 
		std::vector<double>* parHighT, std::vector<double>* Temps);
	void knownMoleculeData(Molecola* mol, std::vector<double>* parLowT,
		std::vector<double>* parHighT, std::vector<double>* Temps);
	void knownMoleculeData(std::string molName, Eigen::VectorXd* parLowT,
		Eigen::VectorXd* parHighT, std::vector<double>* Temps);
	void knownMoleculeData(Molecola* mol, Eigen::VectorXd* parLowT,
		Eigen::VectorXd* parHighT, std::vector<double>* Temps);
	void NASAParameters(Molecola* mol, std::vector<double>* TLimsVec, 
		std::vector<double>* paramLowT, std::vector<double>*  paramHighT);
};