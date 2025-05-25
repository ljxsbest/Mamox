#include "Simulation.h"



Simulation::Simulation(std::string CKmechPath, std::string CKthermoPath,
	std::vector<double> temperatures, std::vector<double> pressures, 
	std::vector<Molecola>* speciesList, ChemkinOut* chemkinOut, int numberOfThreads)
{
	Ts = temperatures;
	Ps = pressures;
	species = speciesList;
	initializeDistributions();
	
	chemOut = chemkinOut;
	numThreads = numberOfThreads;

	std::cout << "Converting mechanism to yaml format." << std::endl;
	try   // try to convert the kinetic mechanism from Chemkin format to CTI
	{
		//Cantera::ck2cti(CKmechPath, CKthermoPath);
		std::string command = "ck2yaml --input=" + CKmechPath +
			" --thermo=" + CKthermoPath + " --permissive --quiet --no-validate";
		system(command.c_str());
		std::stringstream CKmechpatSS(CKmechPath);
		std::getline(CKmechpatSS, modelPath, '.');
		modelPath.append(".yaml");
	}
	catch (Cantera::CanteraError& err)
	{
		std::cout << err.what() << std::endl;
		UTL::fatalError("Error in generating cantera cti file");
	}

	initializeLumpedSpecies();
}

Simulation::Simulation(std::vector<double> temperatures, std::vector<Molecola>* speciesList,
	ChemkinOut* chemkinOut, ThermoOut* thermoOut, int numberOfThreads)
{
	Ts = temperatures;
	Ps = std::vector<double>{ 1.0 };
	species = speciesList;
	initializeDistributions();
	chemOut = chemkinOut;
	numThreads = numberOfThreads;
	thermOut = thermoOut;
	hasParameterBeenSet = true;
	isThermoBased = true;
	initializeLumpedSpecies();
}

void Simulation::initializeLumpedSpecies()
{
	if (species == NULL)
	{
		UTL::error("Simulation::initializeLumpedSpecies called before defining species vector");
		return;
	}
	if (chemOut == NULL)
	{
		UTL::error("Simulation::initializeLumpedSpecies called before defining the ChemkinOut pointer");
		return;
	}
	// generate list of lumped species and the map vector
	for (auto& spec : *species)
		UTL::addUnique(&lumpedSpecies, chemOut->molToNameLump(spec));
	lumpSpecMap.resize(lumpedSpecies.size());
	for (int i = 0; i < lumpedSpecies.size(); i++)
	{
		std::vector<int> temp;
		for (int j = 0; j < species->size(); j++)
			if (lumpedSpecies[i] == chemOut->molToNameLump((*species)[j]))
				temp.push_back(j);
		lumpSpecMap[i] = temp;
	}

	// generate list of species names
	speciesNames.resize(species->size());
	for (int i = 0; i < species->size(); i++)
		speciesNames[i] = chemOut->molToName((*species)[i]);
}


void Simulation::initializeDistributions()
{
	if (species == NULL)
	{
		UTL::error("Simulation::initializeDistributions called before defining species vector");
		return;
	}
	if (Ps.size() == 0)
	{
		UTL::error("Simulation::initializeDistributions called before defining pressure vector");
		return;
	}
	if (Ts.size() == 0)
	{
		UTL::error("Simulation::initializeDistributions called before defining temperature vector");
		return;
	}
	std::vector<double> tempVec1(species->size(), 0);
	std::vector<std::vector<double>> tempVec2(Ps.size(), tempVec1);
	for (int i = 0; i < Ts.size(); i++)
	{
		//concentrations.push_back(tempVec2);
		distributions.push_back(tempVec2);
	}
}

void Simulation::setCSTRParameters(Molecola fuel_, double eqRatio_, double tau_, 
	double simTime_, double maxTime_)
{
	if (isThermoBased)
	{
		UTL::warning("Simulation::setCSTRParameters cannot be called when simulation is based on thermodynamic.");
		return;
	}
	useBatch = false;
	fuel = fuel_;
	eqRatio = eqRatio_;
	tau = tau_;
	simTime = simTime_;
	hasParameterBeenSet = true;
	maxTime = maxTime_;
}

void Simulation::setBatchParameters(Molecola fuel_, double eqRatio_, double simTime_,
	double maxSimTime)
{
	if (isThermoBased)
	{
		UTL::warning("Simulation::setBatchParameters cannot be called when simulation is based on thermodynamic.");
		return;
	}
	useBatch = true;
	fuel = fuel_;
	eqRatio = eqRatio_;
	simTime = simTime_;
	maxTime = maxSimTime;
	hasParameterBeenSet = true;
}


void Simulation::solve()
{
	if (hasParameterBeenSet == false)
	{
		UTL::warning("Simulation::solve called before setting parameters");
		return;
	}
	std::vector<std::vector<int>> simToTPMap;
	for (int i = 0; i < Ts.size(); i++)
		for (int j = 0; j < Ps.size(); j++)
			simToTPMap.push_back(std::vector<int> {i, j});

	
	std::vector<std::thread> threads(numThreads);
	std::vector<int> threadStatus(numThreads, 0); // 0 empty, 1 running, 2 finished
	std::vector<int> threadJob(numThreads,0);
	int indexSim = 0;
	int completedSims = 0;
	//std::cout << std::endl;
	std::cout << "Running simulations:" << std::endl;
	while (true)
	{
		//bool endLoop = false;
		for (int i = 0; i < threads.size(); i++)
		{
			if (threadStatus[i] == 0 && indexSim < simToTPMap.size())
			{
				//std::cout << "thread " << i << " started with job " << indexSim;
				int indT = simToTPMap[indexSim][0];
				int indP = simToTPMap[indexSim][1];
				//std::cout << "(T=" << Ts[indT] << ", P=" << Ps[indP] << ")" << std::endl;
				threadStatus[i] = 1;
				if (isThermoBased)
				{
					threads[i] = std::thread(&Simulation::solveThermo, this, Ts[indT],
						&(distributions[indT][indP]), &(threadStatus[i]));
				}
				else
				{
					if (useBatch)
					{
						threads[i] = std::thread(&Simulation::solveBatch, this, Ts[indT],
							Ps[indP], &(distributions[indT][indP]),
							&(threadStatus[i]));
					}
					else
					{
						threads[i] = std::thread(&Simulation::solveCSTR, this, Ts[indT],
							Ps[indP], &(distributions[indT][indP]),
							&(threadStatus[i]));
					}
				}
				threadJob[i] = indexSim;
				indexSim++;
			}
			if (threadStatus[i] == 2)
			{
				//std::cout << "thread " << i << " finished" << std::endl;
				threads[i].join();
				threadStatus[i] = 0;
				completedSims++;
				UTL::setCursorToLineStart();
				std::cout << "Simulations completed = " << completedSims << "/" <<
					simToTPMap.size();
			}	
		}
		if (indexSim == simToTPMap.size())
		{
			bool isFinished = true;
			for (int j = 0; j < threadStatus.size(); j++)
				if (threadStatus[j] != 0)
					isFinished = false;
			if (isFinished)
				break;
		}
		
	}
	std::wcout << std::endl;
	isSimFinished = true;
}

void Simulation::solveCSTR(double T, double P, std::vector<double>* dist, int* status)
{
	*status = 1;
	try
	{
		auto sol = Cantera::newSolution(modelPath, "gas", "None");
		//sol->kinetics()->
		auto gas = sol->thermo();

		double molHC = 1.0;
		double numC = double(fuel.numberOfC());
		double molO2Stoic = molHC * (numC + (numC + 1) / 2.0);
		double molO2 = molO2Stoic / eqRatio;
		double molN2 = molO2 / 0.21 * 0.79;
		double totMol = molHC + molO2 + molN2;
		double XHC = molHC / totMol;
		double XO2 = molO2 / totMol;
		double XN2 = molN2 / totMol;

		std::string inletComposition = "";
		inletComposition.append(chemOut->molToName(fuel));
		inletComposition.append(":");
		inletComposition.append(std::to_string(XHC));
		inletComposition.append(", O2:");
		inletComposition.append(std::to_string(XO2));
		inletComposition.append(", N2:");
		inletComposition.append(std::to_string(XN2));

		gas->setState_TPX(T, P * 101325.0, inletComposition);

		Cantera::Reservoir flowIn;
		flowIn.insert(sol);
		Cantera::Reservoir flowOut;
		flowOut.insert(sol);
		Cantera::IdealGasConstPressureReactor reactor;
		reactor.setEnergy(0);
		reactor.insert(sol);
		reactor.setInitialVolume(1.0);

		double mass = reactor.mass();	// [kg]
		double massFlow = mass / tau;

		Cantera::MassFlowController massFlowIn;
		massFlowIn.setMassFlowRate(massFlow);

		Cantera::MassFlowController massFlowOut;
		massFlowOut.setMassFlowRate(massFlow);

		massFlowIn.install(flowIn, reactor);
		massFlowOut.install(reactor, flowOut);


		Cantera::ReactorNet sim;
		//sim.setTolerances(1e-12, 1e-20);
		sim.addReactor(reactor);
		int maxSteps = sim.maxSteps();
		sim.setMaxSteps(maxSteps * 100);

		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		std::chrono::steady_clock::time_point now;

		while (sim.time() < simTime)
		{
			now = std::chrono::steady_clock::now();
			if (std::chrono::duration_cast<std::chrono::seconds>(now - begin).count()
				> maxTime)
			{
				UTL::warning("Maximum simulation time exceeded! Simulated up to " 
					+ std::to_string(sim.time()) + " [s]. (" + std::to_string(T) 
					+ "K," + std::to_string(P) + "atm)");
				break;
			}
			sim.step();
		}
		//std::cout << "Simulation at T = " << T << "[K] and P = " << P << "[atm], completed" << std::endl;
		int numSp = gas->nSpecies();
		Cantera::ThermoPhase& c = reactor.contents();
		double totConc = P * 101325 / 8.3145 / T;
		std::vector<double> concentrations = *dist;
		for (int i = 0; i < species->size(); i++)
		{
			std::string specName = chemOut->molToName((*species)[i]);
			for (int j = 0; j < numSp; j++)
			{
				if (specName == c.speciesName(j))
				{
					concentrations[i] = c.moleFraction(j) * totConc;
					break;
				}
			}
		}


		for (int i = 0; i < lumpedSpecies.size(); i++)
		{
			double sum = 0;
			for (int j = 0; j < lumpSpecMap[i].size(); j++)
				sum += concentrations[lumpSpecMap[i][j]];
			for (int j = 0; j < lumpSpecMap[i].size(); j++)
			{
				if (sum == 0)
					(*dist)[lumpSpecMap[i][j]] = 1.0 /
					lumpSpecMap[i].size();
				else
					(*dist)[lumpSpecMap[i][j]] =
					concentrations[lumpSpecMap[i][j]] / sum;
			}
		}
	}
	catch (Cantera::CanteraError& err) {
		std::cout << err.what() << std::endl;
		UTL::error("Error in solving CSTR simulation");
	}

	*status = 2;
}

void Simulation::solveBatch(double T, double P, std::vector<double>* dist, int* status)
{
	*status = 1;
	try
	{
		auto sol = Cantera::newSolution(modelPath, "gas", "None");
		auto gas = sol->thermo();

		double molHC = 1.0;
		double numC = double(fuel.numberOfC());
		double molO2Stoic = molHC * (numC + (numC + 1) / 2.0);
		double molO2 = molO2Stoic / eqRatio;
		double molN2 = molO2 / 0.21 * 0.79;
		double totMol = molHC + molO2 + molN2;
		double XHC = molHC / totMol;
		double XO2 = molO2 / totMol;
		double XN2 = molN2 / totMol;

		std::string inletComposition = "";
		inletComposition.append(chemOut->molToName(fuel));
		inletComposition.append(":");
		inletComposition.append(std::to_string(XHC));
		inletComposition.append(", O2:");
		inletComposition.append(std::to_string(XO2));
		inletComposition.append(", N2:");
		inletComposition.append(std::to_string(XN2));

		gas->setState_TPX(T, P * 101325.0, inletComposition);

		Cantera::IdealGasReactor reactor;
		reactor.setEnergy(1);
		reactor.insert(sol);
		reactor.setInitialVolume(1.0);

		Cantera::ReactorNet sim;
		sim.addReactor(reactor);
		double simRelTol = sim.rtol();
		//std::cout << sim.rtol() << std::endl;
		//std::cout << sim.atol() << std::endl;
		//sim.setTolerances(simRelTol, 1.0e-18);
		Cantera::ThermoPhase* content;
		content = &(reactor.contents());
		int numSp = gas->nSpecies();
		std::vector<int> speciesMap(species->size());
		for (int i = 0; i < species->size(); i++)
		{
			std::string specName = chemOut->molToName((*species)[i]);
			for (int j = 0; j < numSp; j++)
			{
				if (specName == content->speciesName(j))
				{
					speciesMap[i] = j;
					break;
				}
			}
		}
		int OHind = 0;
		for (int j = 0; j < numSp; j++)
		{
			if ("OH" == content->speciesName(j))
			{
				OHind = j;
				break;
			}
		}
		double prevTime = 0.0;
		double maxOHconc = 0.0;

		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		std::chrono::steady_clock::time_point now;
		
		std::vector<double> times; // debug
		std::vector<double> OHs;   // debug
		std::vector<double> Ts;   // debug
		
		while (true)
		{
			sim.step();
			double timeStep = sim.time() - prevTime;
			prevTime = sim.time();
			double OHconc = content->moleFraction(OHind);
			times.push_back(sim.time()); // debug
			OHs.push_back(OHconc);		 // debug
			Ts.push_back(content->temperature());		 // debug
			double deltaT = content->temperature() - T;
			if (OHconc > maxOHconc)
				maxOHconc = OHconc;
			if (OHconc < maxOHconc/(1+simRelTol) && OHconc > 0.0 && deltaT > 10.0)
				break;
			if (deltaT > 500.0)
				break;
			for (int i = 0; i < lumpedSpecies.size(); i++)
			{
				double sum = 0;
				for (int j = 0; j < lumpSpecMap[i].size(); j++)
				{
					double molF = content->moleFraction(speciesMap[lumpSpecMap[i][j]]);
					if (molF < 0.0)
						molF = 0.0;
					sum += molF;
				}
				for (int j = 0; j < lumpSpecMap[i].size(); j++)
				{
					if (sum == 0)
						(*dist)[lumpSpecMap[i][j]] += 1.0 /
						lumpSpecMap[i].size() * timeStep;
					else
					{
						double molF = content->moleFraction(speciesMap[lumpSpecMap[i][j]]);
						if (molF < 0.0)
							molF = 0.0;

						(*dist)[lumpSpecMap[i][j]] += molF / sum * timeStep;
					}
				}
			}
			if (sim.time() > simTime)
			{
				UTL::warning("After " + std::to_string(sim.time()) + " simulated"
					+ "seconds ignition has not been detected. Simulation stopped.T:" 
					+ std::to_string(T) + ", P : " + std::to_string(P));
			}
			now = std::chrono::steady_clock::now();
			if (std::chrono::duration_cast<std::chrono::seconds>(now - begin).count()
				> maxTime)
			{
				UTL::warning("Maximum simulation time exceeded! Simulated up to "
					+ std::to_string(sim.time()) + " [s]. (" + std::to_string(T)
					+ "K," + std::to_string(P) + "atm)");
				break;
			}
		}
		for (int i = 0; i < dist->size(); i++)
			(*dist)[i] = (*dist)[i] / sim.time();
	}
	catch (Cantera::CanteraError& err) {
		std::cout << err.what() << std::endl;
		UTL::error("Error in solving batch simulation");
	}

	*status = 2;
}


void Simulation::solveThermo(double T, std::vector<double>* distr, int* status)
{
	*status = 1;
	for (int lumInd = 0; lumInd < lumpedSpecies.size(); lumInd++)
	{
		std::vector<Molecola> groupedSpecies(lumpSpecMap[lumInd].size());
		for (int i = 0; i < lumpSpecMap[lumInd].size(); i++ )
			groupedSpecies[i] = (*species)[lumpSpecMap[lumInd][i]];

		std::vector<double> DGri(groupedSpecies.size() - 1);
		std::vector<double> Keqi(groupedSpecies.size() - 1);

		for (int i = 0; i < groupedSpecies.size() - 1; i++)
		{
			double Hfip1 = thermOut->Hf(groupedSpecies[i + 1], T);
			double Hfi = thermOut->Hf(groupedSpecies[i], T);
			double Sip1 = thermOut->S(groupedSpecies[i + 1], T);
			double Si = thermOut->S(groupedSpecies[i], T);
			DGri[i] = (Hfip1 - Hfi) - T * (Sip1 - Si);
			Keqi[i] = std::exp(-DGri[i] / 1.987 / T);
		}
		std::vector<double> distrLumped(groupedSpecies.size());
		distrLumped[0] = 0.0;
		for (int i = 0; i < distrLumped.size(); i++)
		{
			double productory = 1;
			for (int j = 0; j < i; j++)
			{
				productory *= Keqi[j];
			}
			distrLumped[0] += productory;
		}
		distrLumped[0] = 1.0 / distrLumped[0];
		for (int i = 0; i < distrLumped.size() - 1; i++)
		{
			distrLumped[i + 1] = distrLumped[i] * Keqi[i];
		}
		for (int i = 0; i < lumpSpecMap[lumInd].size(); i++)
		{
			(*distr)[lumpSpecMap[lumInd][i]] = distrLumped[i];
		}
	}
	*status = 2;
}

void Simulation::printResults(std::string folderPath)
{
	if (UTL::dirExists(folderPath) == false)
	{
		UTL::warning("In Simulation::printResults the provided folder does not exist, it will be gneerated");
		UTL::createDirectory(folderPath);
	}
	for (int Tind = 0; Tind < Ts.size(); Tind++)
	{
		for (int Pind = 0; Pind < Ps.size(); Pind++)
		{
			std::ofstream outFile(folderPath + "\\T" + std::to_string(Ts[Tind]) +
				"_P" + std::to_string(Ps[Pind]) + ".csv");
			outFile << "Species,Lumped species,Distribution[-]"
				<< std::endl;
			for (int i = 0; i < species->size(); i++)
			{
				char buffer[100];
				sprintf_s(buffer, 100, "\"%16s\",\"%16s\",%15e", 
					chemOut->molToName((*species)[i]), 
					chemOut->molToNameLump((*species)[i]),
					distributions[Tind][Pind][i]);
				outFile << buffer << std::endl;
			}
			outFile.close();
		}
	}
}


//double Simulation::specConcentration(int Tind, int Pind, std::string name)
//{
//	for (int i = 0; i < speciesNames.size(); i++)
//		if (speciesNames[i] == name)
//			return concentrations[Tind][Pind][i];
//	UTL::warning("Simulation::specConcentration, species not found. 0 returned.");
//	std::cout << name << std::endl;
//	return 0;
//}
//double Simulation::specConcentration(int Tind, int Pind, Molecola mol)
//{
//	std::string name = chemOut->molToName(mol);
//	return specConcentration(Tind, Pind, name);
//}
double Simulation::specDistribution(int Tind, int Pind, std::string name)
{
	for (int i = 0; i < speciesNames.size(); i++)
		if (speciesNames[i] == name)
			return distributions[Tind][Pind][i];
	UTL::warning("Simulation::specConcentration, species not found. 0 returned.");
	std::cout << name << std::endl;
	return 0;
}
double Simulation::specDistribution(int Tind, int Pind, Molecola mol)
{
	std::string name = chemOut->molToName(mol);
	return specDistribution(Tind, Pind, name);
}