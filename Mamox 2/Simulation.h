#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <thread>
#include "molec.h"
#include "chemkinOut.h"
#include "utilities.h"
#include "thermoOut.h"
#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h"
//#include "cantera/base/ctml.h"
#include "cantera/zeroD/flowControllers.h"
#include "cantera/zeroD/Reservoir.h"

class Simulation
{
public:
	Simulation(std::string CKmechPath, std::string CKthermoPath,
		std::vector<double> temperatures, std::vector<double> pressures, 
		std::vector<Molecola>* speciesList,	ChemkinOut* chemkinOut, int numberOfThreads);
	Simulation(std::vector<double> temperatures, std::vector<Molecola>* speciesList,
		ChemkinOut* chemkinOut, ThermoOut* thermoOut, int numberOfThreads);
	void setBatchParameters(Molecola fuel_, double eqRatio, double simTime, double maxSimTime);
	void setCSTRParameters(Molecola fuel_, double eqRatio_, double tau_, double simTime_,
		double maxTime_);
	void solve();
	void printResults(std::string folderPath); // print the results in the provided folder
	std::vector<double> Temperatures() { return Ts; };
	std::vector<double> Pressures() { return Ps; };
	//double specConcentration(int Tind, int Pind, std::string name);
	//double specConcentration(int Tind, int Pind, Molecola mol);
	double specDistribution (int Tind, int Pind, std::string name);
	double specDistribution (int Tind, int Pind, Molecola mol);

private:
	std::vector<double> Ts;			// simulations' temperatures [K]
	std::vector<double> Ps;			// simulations' pressures [atm]
	//std::vector<std::vector<std::vector<double>>> concentrations;	// [mol/m3]
	std::vector<std::vector<std::vector<double>>> distributions;	
	std::vector<Molecola>* species = NULL;
	std::vector<std::string> speciesNames;
	bool useBatch;
	bool isThermoBased = false;
	std::string modelPath;
	ChemkinOut* chemOut = NULL;
	ThermoOut* thermOut = NULL;
	bool hasParameterBeenSet = false;
	bool isSimFinished = false;
	Molecola fuel;
	double eqRatio;
	double tau;
	double simTime;
	double maxTime;
	int numThreads = 1;
	void solveCSTR(double T, double P, std::vector<double>* distr, int* status);
	void solveBatch(double T, double P, std::vector<double>* distr, int* status);
	void solveThermo(double T, std::vector<double>* distr, int* status);
	std::vector<std::string>  lumpedSpecies; // list of all the different lumped species
	std::vector<std::vector<int>> lumpSpecMap;	// for each lumped species of lumpedSpecies
												// list the index of species in species
												// vector belonging to that class
	void initializeDistributions();		// initialize the matrix of distributions 
										// with the right size filled with 0 values
	void initializeLumpedSpecies();		// generate the list of lumped species and 
										// the map between lumped species and normal 
										// species
};

