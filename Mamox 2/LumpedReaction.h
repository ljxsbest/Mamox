#pragma once
#include <vector>
#include <iostream>
#include <string>
#include "Reaction.h"
#include "chemkinOut.h"
#include "Simulation.h"
#include "thermoOut.h"
#include "utilities.h"
#include "nlopt.h"

class LumpedReaction
{
public:
	LumpedReaction(std::vector<Reaction> reactions, Simulation* sim,
		ThermoOut* thermOut,	ChemkinOut* chemkinout);
	LumpedReaction(std::vector<Reaction> reactions, ThermoOut* thermOut,
		std::vector<double> temperatures, ChemkinOut* chemkinout);
	std::string print();
	std::string relevantReactantName;
	std::string relevantReactantThermo;
	std::string printTargetRates();

private:
	std::vector<Reaction> reacs;
	std::vector<double> Pressures;
	std::vector<double> targetTs;
	std::vector<std::vector<double>> targetRates;
	std::vector<double> As;
	std::vector<double> ns;
	std::vector<double> Es;
	bool errorOccured = false;
	ChemkinOut* chemOut;
	bool needsPLOG = true;
	std::vector<std::string> reacNames = {};
	std::vector<std::string> prodNames = {};
	std::vector<Molecola> sampleProducts = {};
	std::vector<Molecola> sampleReactants = {};
	std::vector<double> stoicCoeff = {};
	void computeFirstGuess(double x[], std::vector<double> temps, std::vector<double> ks);
	//double optimizationFunction(unsigned n,	const double* x, double* grad, 
	//	void* my_func_data);
};

