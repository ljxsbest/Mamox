#include "LumpedReaction.h"
double optimizationFunction(unsigned n,
	const double* x, double* grad, void* my_func_data)
{
	if (grad) {
		std::cerr << "Algorithm for optimization requires gradient!" << std::endl;
	}

	std::vector<double>* p = static_cast<std::vector<double>*>(my_func_data);

	int ind = 0;
	int numberOfT = (*p)[ind++];		// number of temperatures

	std::vector<double> T(numberOfT);
	std::vector<double> v(numberOfT);
	for (int i = 0; i < numberOfT; i++) T[i] = (*p)[ind++];		// temperatures
	for (int i = 0; i < numberOfT; i++) v[i] = (*p)[ind++];		// actual reaction velocities

	double R_ = 1.987;

	double func = 0.0;

	for (int i = 0; i < numberOfT; i++)
	{
		if (v[i] != 0.0)
		{
			func += pow((log(v[i]) - (x[0] + x[1] * log(T[i]) - x[2] / R_ / T[i])) / log(v[i]), 2.0);
			//func += abs(abs((log(v[i]) - (x[0] + x[1] * log(T[i]) - x[2] / R_ / T[i])) / log(v[i])));
			//func += abs((log(v[i]) - (x[0] + x[1] * log(T[i]) - x[2] / R_ / T[i])));
		}
		else
		{
			func += abs(x[0] + x[1] * log(T[i]) - x[2] / R_ / T[i]);
		}
	}


	return func;

}

LumpedReaction::LumpedReaction(std::vector<Reaction> reactions, Simulation* sim,
	ThermoOut* thermOut,	ChemkinOut* chemkinout)
{
	std::vector<double> Ts = sim->Temperatures();
	Pressures = sim->Pressures();
	reacs = reactions;
	chemOut = chemkinout;

	targetRates.resize(Pressures.size());
	for (int i = 0; i < Pressures.size(); i++)
		targetRates[i] = std::vector<double>(Ts.size(), 0.0);
	targetTs = Ts;

	if (reactions.size() == 0)
	{
		UTL::error("LumpedReaction constructor called on an empty reactions vector");
		errorOccured = true;
		return;
	}



	// check if reactions are lumpable
	std::vector<Molecola> reactants;
	std::vector<Molecola*> reactantsPointers = reacs[0].reactantList();
	for (auto& reP : reactantsPointers)
	{
		reacNames.push_back(chemOut->molToNameLump(*reP));
		reactants.push_back(*reP);
		sampleReactants.push_back(*reP);
	}
	for (auto& R : reacs)
	{
		std::vector<std::string> reac2Lump;
		std::vector<Molecola*> reacPoin2 = R.reactantList();
		for (auto& reP : reacPoin2)
		{
			reac2Lump.push_back(chemOut->molToNameLump(*reP));
		}
		if (UTL::areEquivalent(&reacNames, &reac2Lump) == false)
		{
			UTL::error("LumpedReaction called on non lumpable reactions");
			errorOccured = true;
			return;
		}
	}


	needsPLOG = false;
	for (auto& R : reacs)
	{
		std::vector<Molecola> reac2;
		std::vector<Molecola*> reacPoin2 = R.reactantList();
		for (auto& reP : reacPoin2)
			reac2.push_back(*reP);
		if (UTL::areEquivalent(&reactants, &reac2) == false)
		{
			needsPLOG = true;
			break;
		}
	}

	std::vector<Molecola> relevReac(reactions.size());
	std::vector<Molecola> commonReactant;
	if (needsPLOG)
	{
		// analyze the reactants and make the list of the relevant isomers
		for (auto& reac : reactants)	// detect a common reactants (if present)
		{
			bool isPresentInAll = true;
			for (auto& r : reactions)
			{
				std::vector<Molecola*> reaPointList = r.reactantList();
				std::vector<Molecola> reaList;
				for (auto& reap : reaPointList)
					reaList.push_back(*reap);
				if (UTL::isPresent(&reaList, reac) == false)
				{
					isPresentInAll = false;
					break;
				}
			}
			if (isPresentInAll == true)
				commonReactant.push_back(reac);
		}
		if (commonReactant.size() > 2)
		{
			UTL::error("In LumpedReaction more than one common reactant has been found");
			errorOccured = true;
			return;
		}
		// build the list of relevant reactants (the one we need the distribution of)
		for (int i = 0; i < reactions.size(); i++)
		{
			std::vector<Molecola*> reaPointList = reactions[i].reactantList();
			for (auto& reaP : reaPointList)
			{
				bool hasBeenFound = false;
				for (auto& commR : commonReactant)
					if (*reaP == commR)
						hasBeenFound = true;
				if (hasBeenFound == false)
				{
					relevReac[i] = *reaP;
					continue;
				}
			}
		}
	}
	else
	{
		std::vector<Molecola*> r = reactions[0].reactantList();
		std::vector<Molecola> HCs;
		for (auto& spec : r)
			if (spec->isSpecialMolecule() == 0)
				HCs.push_back(*spec);
		int maxSize = 0;
		Molecola relR;
		for (auto& spec : HCs)
		{
			if (spec.size() > maxSize)
			{
				maxSize = spec.size();
				relR = spec;
			}
		}
		for (auto& RR : relevReac)
			RR = relR;
	}

	//std::cout << "relev reac" << std::endl;
	//for (auto& rea : relevReac)
	//	std::cout << rea << std::endl;

	// build list of products
	for (auto& reac : reacs)
	{
		std::vector<Molecola*> prodPoint = reac.productList();
		for (auto& PP : prodPoint)
		{
			int returnValue = UTL::addUnique(&prodNames, chemOut->molToNameLump(*PP));
			if (returnValue == -1)
				sampleProducts.push_back(*PP);
		}
	}

	int indexP1 = -1;   // find the pressure value that is closer to 1
	{
		double minDiff = 1000000000000;
		for (int i = 0; i < Pressures.size(); i++)
		{
			if (std::abs(1 - Pressures[i]) < minDiff)
			{
				minDiff = std::abs(1 - Pressures[i]);
				indexP1 = i;
			}
		}
	}
	relevantReactantName = chemOut->molToNameLump(relevReac[0]);
	{
		std::vector<Molecola> uniqueRelevReac;
		UTL::concatenateUnique(&uniqueRelevReac, &relevReac);
		std::vector<std::vector<double>> distrForThermo(uniqueRelevReac.size());
		for (auto& DFT : distrForThermo)
			DFT = std::vector<double>(Ts.size());
		for (int Tind = 0; Tind < Ts.size(); Tind++)
		{
			double tot = 0.0;
			for (int Sind = 0; Sind < uniqueRelevReac.size(); Sind++)
			{
				//double conc = sim->specConcentration(Tind, indexP1,
				//	uniqueRelevReac[Sind]);
				double conc = sim->specDistribution(Tind, indexP1,
					uniqueRelevReac[Sind]);
				distrForThermo[Sind][Tind] = conc;
				tot += conc;
			}
			for (int Sind = 0; Sind < uniqueRelevReac.size(); Sind++)
			{
				distrForThermo[Sind][Tind] = distrForThermo[Sind][Tind]/tot;
			}
		}
		relevantReactantThermo = thermOut->NASAOutputLumped(relevantReactantName,
			uniqueRelevReac, Ts, distrForThermo);
	}
	// compute stoichiometric coefficients
	{
		std::vector<double> totRate(Ts.size());
		std::vector<std::vector<double>> prodRate(Ts.size());
		for (auto& PR : prodRate)
			PR = std::vector<double>(prodNames.size(), 0.0);
		std::vector<std::vector<double>> prodDist(Ts.size());
		for (auto& PD : prodDist)
			PD = std::vector<double>(prodNames.size(), 0.0);
		for (int Tind = 0; Tind < Ts.size(); Tind++)
		{
			for (int reaInd = 0; reaInd < reacs.size(); reaInd++)
			{
				double rate = reacs[reaInd].rateConstant(Ts[Tind]) *
					sim->specDistribution(Tind, indexP1, relevReac[reaInd]);
				totRate[Tind] += rate;
				std::vector<Molecola*> prods = reacs[reaInd].productList();
				for (auto& pr : prods)
				{

					int ind = UTL::addUnique(&prodNames,
						chemOut->molToNameLump(*pr));
					prodRate[Tind][ind] += rate;
				}
			}
			for (int prodInd = 0; prodInd < prodNames.size(); prodInd++)
			{
				prodDist[Tind][prodInd] = prodRate[Tind][prodInd] / totRate[Tind];
			}
		}
		// TO DO: make a proper integration of the rates and distributions over T
		double sumTotRate = 0.0;
		for (auto& tr : totRate)
			sumTotRate += tr;
		stoicCoeff = std::vector<double>(prodNames.size(), 0.0);
		for (int spInd = 0; spInd < prodNames.size(); spInd++)
		{
			for (int Tind = 0; Tind < Ts.size(); Tind++)
			{
				stoicCoeff[spInd] += prodDist[Tind][spInd] * totRate[Tind] / sumTotRate;
			}
		}
	}

	//std::cout << "STOICHIOMETRIC COEFFICIENT" << std::endl;
	//for (int i = 0; i < prodNames.size(); i++)
	//	std::cout << prodNames[i] << "   " << stoicCoeff[i] << std::endl;


	if (!needsPLOG)
	{
		Pressures = { 1.0 };
	}


	for (int Pind = 0; Pind < Pressures.size(); Pind++)
	{
		std::vector<double> targetKs(Ts.size());
		for (int Tind = 0; Tind < Ts.size(); Tind++)
		{
			double k = 0.0;
			//std::cout << "reactoins" << reactions.size() << std::endl;
			//std::cout << "reacs" << reacs.size() << std::endl;
			for (int j = 0; j < reactions.size(); j++)
			{
				//std::cout << j;
				//std::cout << "5 - " << relevReac[j] << std::endl;
				k += sim->specDistribution(Tind, Pind, relevReac[j])
					* reactions[j].rateConstant(Ts[Tind]);
			}
			targetKs[Tind] = k;
		}
		targetRates[Pind] = targetKs;
		nlopt_opt opt_loc;
		opt_loc = nlopt_create(NLOPT_LN_BOBYQA, 3);
		std::vector<double> dataToFunction(Ts.size() * 2 + 1);	// vector of the data to feed to the funciont to optimize
		// populate data to function
		//		first element is the number of temperatures
		//		then there are all the temperatures
		//		then there are all the reaction rates
		dataToFunction[0] = Ts.size();
		int ind = 1;
		for (auto& t : Ts)
		{
			dataToFunction[ind] = t;
			ind++;
		}
		for (auto& k : targetKs)
		{
			dataToFunction[ind] = k;
			ind++;
		}

		nlopt_set_min_objective(opt_loc, optimizationFunction, &dataToFunction);

		double absTol = 1.0e-15;
		//nlopt_set_xtol_abs(opt_loc, &absTol);
		nlopt_set_xtol_rel(opt_loc, 1.0e-15);
		nlopt_set_maxtime(opt_loc, 10.);

		// convert the pre exponential factor from cm3 to m3
		double kinConst[3];
		//kinConst[0] = std::log(reacs[0].A());
		//kinConst[1] = reacs[0].n();
		//kinConst[2] = reacs[0].E();
		computeFirstGuess(kinConst, Ts, targetKs);
		//std::cout << kinConst[0] << "  " << kinConst[1] << "  " << kinConst[2] << std::endl;
		double minf = 0.;
		int exitFlag = nlopt_optimize(opt_loc, kinConst, &minf);
		//std::cout << "exit flag " << exitFlag << std::endl;
		if (exitFlag < 0 && exitFlag != -4)
			UTL::error("Failed to fit kinetic parameters at " +
				std::to_string(Pressures[Pind]) + " atm (exit flag " +
				std::to_string(exitFlag) + ")");
		//std::cout << kinConst[0] << "  " << kinConst[1] << "  " << kinConst[2] << std::endl;
		As.push_back(std::exp(kinConst[0]));
		ns.push_back(kinConst[1]);
		Es.push_back(kinConst[2]);
		nlopt_destroy(opt_loc);
	}



}


LumpedReaction::LumpedReaction(std::vector<Reaction> reactions, ThermoOut* thermOut,
	std::vector<double> temperatures, ChemkinOut* chemkinout)
{
	std::vector<double> Ts = temperatures;
	Pressures = { 1 };
	reacs = reactions;
	chemOut = chemkinout;
	needsPLOG = false;

	targetRates.resize(Pressures.size());
	for (int i = 0; i < Pressures.size(); i++)
		targetRates[i] = std::vector<double>(Ts.size(), 0.0);
	targetTs = Ts;

	if (reactions.size() == 0)
	{
		UTL::error("LumpedReaction constructor called on an empty reactions vector");
		errorOccured = true;
		return;
	}

	// check if reactions are lumpable
	std::vector<Molecola> reactants;
	std::vector<Molecola*> reactantsPointers = reacs[0].reactantList();
	for (auto& reP : reactantsPointers)
	{
		reacNames.push_back(chemOut->molToNameLump(*reP));
		reactants.push_back(*reP);
	}
	for (auto& R : reacs)
	{
		std::vector<std::string> reac2Lump;
		std::vector<Molecola*> reacPoin2 = R.reactantList();
		for (auto& reP : reacPoin2)
		{
			reac2Lump.push_back(chemOut->molToNameLump(*reP));
		}
		if (UTL::areEquivalent(&reacNames, &reac2Lump) == false)
		{
			UTL::error("LumpedReaction called on non lumpable reactions");
			errorOccured = true;
			return;
		}
	}

	// detect if the relevant reactant is the same for all reactions
	bool sameRelReac = true;
	for (auto& R : reacs)
	{
		std::vector<Molecola> reac2;
		std::vector<Molecola*> reacPoin2 = R.reactantList();
		for (auto& reP : reacPoin2)
			reac2.push_back(*reP);
		if (UTL::areEquivalent(&reactants, &reac2) == false)
		{
			sameRelReac = false;
			break;
		}
	}
	std::vector<Molecola> relevReac(reactions.size());
	std::vector<Molecola> commonReactant;
	if(sameRelReac)
	{
		std::vector<Molecola*> r = reactions[0].reactantList();
		std::vector<Molecola> HCs;
		for (auto& spec : r)
			if (spec->isSpecialMolecule() == 0)
				HCs.push_back(*spec);
		int maxSize = 0;
		Molecola relR;
		for (auto& spec : HCs)
		{
			if (spec.size() > maxSize)
			{
				maxSize = spec.size();
				relR = spec;
			}
		}
		for (auto& RR : relevReac)
			RR = relR;
	}
	else
	{
		// analyze the reactants and make the list of the relevant isomers
		for (auto& reac : reactants)	// detect a common reactants (if present)
		{
			bool isPresentInAll = true;
			for (auto& r : reactions)
			{
				std::vector<Molecola*> reaPointList = r.reactantList();
				std::vector<Molecola> reaList;
				for (auto& reap : reaPointList)
					reaList.push_back(*reap);
				if (UTL::isPresent(&reaList, reac) == false)
				{
					isPresentInAll = false;
					break;
				}
			}
			if (isPresentInAll == true)
				commonReactant.push_back(reac);
		}
		if (commonReactant.size() > 2)
		{
			UTL::error("In LumpedReaction more than one common reactant has been found");
			errorOccured = true;
			return;
		}
		// build the list of relevant reactants (the one we need the distribution of)
		for (int i = 0; i < reactions.size(); i++)
		{
			std::vector<Molecola*> reaPointList = reactions[i].reactantList();
			for (auto& reaP : reaPointList)
			{
				bool hasBeenFound = false;
				for (auto& commR : commonReactant)
					if (*reaP == commR)
						hasBeenFound = true;
				if (hasBeenFound == false)
				{
					relevReac[i] = *reaP;
					continue;
				}
			}
		}
	}

	std::vector<std::vector<double>> relevReacDist(temperatures.size());
	for (auto& RRD : relevReacDist)
		RRD = std::vector<double>(relevReac.size(),0.0);

	{
		std::vector<Molecola> uniqueRelevReac;
		UTL::concatenateUnique(&uniqueRelevReac, &relevReac);
		std::vector<double> DGri(uniqueRelevReac.size()-1);
		std::vector<double> Keqi(uniqueRelevReac.size()-1);
		//std::vector<std::vector<double>> uniqReacDist(uniqueRelevReac.size(),0.0);
		std::vector<std::vector<double>> uniqReacDist(temperatures.size());
		for (auto& URD : uniqReacDist)
			URD = std::vector<double>(uniqueRelevReac.size(), 0.0);
		for (int Tind = 0; Tind < temperatures.size(); Tind++)
		{
			double T = temperatures[Tind];
			for (int i = 0; i < uniqueRelevReac.size() - 1; i++)
			{
				double Hfip1 = thermOut->Hf(uniqueRelevReac[i + 1], T);
				double Hfi = thermOut->Hf(uniqueRelevReac[i], T);
				double Sip1 = thermOut->S(uniqueRelevReac[i + 1], T);
				double Si = thermOut->S(uniqueRelevReac[i], T);
				DGri[i] = (Hfip1 - Hfi) - T * (Sip1 - Si);
				Keqi[i] = std::exp(-DGri[i] / 1.987 / T);
			}
			uniqReacDist[Tind][0] = 0.0;
			for (int i = 0; i < uniqReacDist[Tind].size(); i++)
			{
				double productory = 1;
				for (int j = 0; j < i; j++)
				{
					productory *= Keqi[j];
				}
				uniqReacDist[Tind][0] += productory;
			}
			uniqReacDist[Tind][0] = 1.0/ uniqReacDist[Tind][0];
			for (int i = 0; i < uniqReacDist[Tind].size()-1; i++)
			{
				uniqReacDist[Tind][i + 1] = uniqReacDist[Tind][i] * Keqi[i];
			}
			for (int specInd = 0; specInd < relevReacDist[Tind].size(); specInd++)
			{
				for (int specInd2 = 0; specInd2 < uniqueRelevReac.size(); specInd2++)
				{
					if (relevReac[specInd] == uniqueRelevReac[specInd2])
					{
						relevReacDist[Tind][specInd] = uniqReacDist[Tind][specInd2];
						break;
					}
				}
			}
		}
	}

	//// debug ->
	//std::cout << "DEBUG START" << std::endl;
	//for (int i = 0; i < relevReacDist[0].size(); i++)
	//	std::cout << relevReacDist[0][i] << std::endl;
	//std::cout << "DEBUG FINISH" << std::endl;
	//// <- debug

	// build list of products
	for (auto& reac : reacs)
	{
		std::vector<Molecola*> prodPoint = reac.productList();
		for (auto& PP : prodPoint)
			UTL::addUnique(&prodNames, chemOut->molToNameLump(*PP));
	}

	int indexP1 = -1;   // find the pressure value that is closer to 1
	{
		double minDiff = 1000000000000;
		for (int i = 0; i < Pressures.size(); i++)
		{
			if (std::abs(1 - Pressures[i]) < minDiff)
			{
				minDiff = std::abs(1 - Pressures[i]);
				indexP1 = i;
			}
		}
	}
	// compute stoichiometric coefficients
	{
		std::vector<double> totRate(Ts.size());
		std::vector<std::vector<double>> prodRate(Ts.size());
		for (auto& PR : prodRate)
			PR = std::vector<double>(prodNames.size(), 0.0);
		std::vector<std::vector<double>> prodDist(Ts.size());
		for (auto& PD : prodDist)
			PD = std::vector<double>(prodNames.size(), 0.0);
		for (int Tind = 0; Tind < Ts.size(); Tind++)
		{
			for (int reaInd = 0; reaInd < reacs.size(); reaInd++)
			{
				double rate = reacs[reaInd].rateConstant(Ts[Tind]) *
					relevReacDist[Tind][reaInd];
				totRate[Tind] += rate;
				std::vector<Molecola*> prods = reacs[reaInd].productList();
				for (auto& pr : prods)
				{

					int ind = UTL::addUnique(&prodNames,
						chemOut->molToNameLump(*pr));
					prodRate[Tind][ind] += rate;
				}
			}
			for (int prodInd = 0; prodInd < prodNames.size(); prodInd++)
			{
				prodDist[Tind][prodInd] = prodRate[Tind][prodInd] / totRate[Tind];
			}
		}
		// TO DO: make a proper integration of the rates and distributions over T
		double sumTotRate = 0.0;
		for (auto& tr : totRate)
			sumTotRate += tr;
		stoicCoeff = std::vector<double>(prodNames.size(), 0.0);
		for (int spInd = 0; spInd < prodNames.size(); spInd++)
		{
			for (int Tind = 0; Tind < Ts.size(); Tind++)
			{
				stoicCoeff[spInd] += prodDist[Tind][spInd] * totRate[Tind] / sumTotRate;
			}
		}
	}

	//std::cout << "STOICHIOMETRIC COEFFICIENT" << std::endl;
	//for (int i = 0; i < prodNames.size(); i++)
	//	std::cout << prodNames[i] << "   " << stoicCoeff[i] << std::endl;
	

	if (!needsPLOG)
	{
		Pressures = { 1.0 };
	}


	for (int Pind = 0; Pind < Pressures.size(); Pind++)
	{
		std::vector<double> targetKs(Ts.size());
		for (int Tind = 0; Tind < Ts.size(); Tind++)
		{
			double k = 0.0;
			//std::cout << "reactoins" << reactions.size() << std::endl;
			//std::cout << "reacs" << reacs.size() << std::endl;
			for (int j = 0; j < reactions.size(); j++)
			{
				//std::cout << j;
				//std::cout << "5 - " << relevReac[j] << std::endl;
				k += relevReacDist[Tind][j]
					* reactions[j].rateConstant(Ts[Tind]);
			}
			targetKs[Tind] = k;
		}
		targetRates[Pind] = targetKs;
		nlopt_opt opt_loc;
		opt_loc = nlopt_create(NLOPT_LN_BOBYQA, 3);
		std::vector<double> dataToFunction(Ts.size() * 2 + 1);	// vector of the data to feed to the funciont to optimize
		
		// populate data to function
		//		first element is the number of temperatures
		//		then there are all the temperatures
		//		then there are all the reaction rates
		dataToFunction[0] = Ts.size();
		int ind = 1;
		for (auto& t : Ts)
		{
			dataToFunction[ind] = t;
			ind++;
		}
		for (auto& k : targetKs)
		{
			dataToFunction[ind] = k;
			ind++;
		}

		nlopt_set_min_objective(opt_loc, optimizationFunction, &dataToFunction);

		double absTol = 1.0e-15;
		//nlopt_set_xtol_abs(opt_loc, &absTol);
		nlopt_set_xtol_rel(opt_loc, 1.0e-15);
		nlopt_set_maxtime(opt_loc, 10.);

		// convert the pre exponential factor from cm3 to m3
		double kinConst[3];
		//kinConst[0] = std::log(reacs[0].A());
		//kinConst[1] = reacs[0].n();
		//kinConst[2] = reacs[0].E();
		computeFirstGuess(kinConst, Ts, targetKs);

		//std::cout << kinConst[0] << "  " << kinConst[1] << "  " << kinConst[2] << std::endl;
		double minf = 0.;
		int exitFlag = nlopt_optimize(opt_loc, kinConst, &minf);
		//std::cout << "exit flag " << exitFlag << std::endl;
		if (exitFlag < 0 && exitFlag != -4)
			UTL::error("Failed to fit kinetic parameters at " +
				std::to_string(Pressures[Pind]) + " atm (exit flag " + 
				std::to_string(exitFlag) +")");
		if (exitFlag == -4)
			UTL::warning("Fit of kinetic paramenters at " +
				std::to_string(Pressures[Pind]) + " atm stopped because of roundoff errors limited progress");
		//std::cout << kinConst[0] << "  " << kinConst[1] << "  " << kinConst[2] << std::endl;
		As.push_back(std::exp(kinConst[0]));
		ns.push_back(kinConst[1]);
		Es.push_back(kinConst[2]);
		nlopt_destroy(opt_loc);
	}
}

bool compareFunc(std::pair<double, int>& a, std::pair<double, int>& b)
{
	return a.first < b.first;
}

std::string LumpedReaction::print()
{
	if (errorOccured)
	{
		return "!Error occured in lumping this reaction\n";
	}
	std::stringstream out;
	for (int i = 0; i < reacNames.size(); i++)
	{
		out << reacNames[i];
		if (i < reacNames.size() - 1)
			out << "+";
	}
	out << "=>";
	bool firstPrinted = false;

	// rounding algorithm taken from doi.org/10.48550/arXiv.1501.00014
	
	std::vector<double> roundedStoicCoeff(stoicCoeff.size(), 0.0);
	std::vector<double> shortFalls(stoicCoeff.size(), 0.0);
	double digitNumber = 0.0001;
	double integerSum = 0;
	for (int i = 0; i < stoicCoeff.size(); i++)
	{
		roundedStoicCoeff[i] = stoicCoeff[i] / digitNumber;
		integerSum += stoicCoeff[i];
	}
	//std::cout << "integer sum = " << integerSum << std::endl;
	integerSum = std::round(integerSum)/digitNumber;
	double shortFallSum = 0;
	for (int i = 0; i < stoicCoeff.size(); i++)
	{
		shortFalls[i] = roundedStoicCoeff[i] - std::floor(roundedStoicCoeff[i]);
		shortFallSum += shortFalls[i];
	}

	std::vector<std::pair<double, int>> coupledVectors(shortFalls.size());
	for (int i = 0; i < shortFalls.size(); i++)
		coupledVectors[i] = std::pair<double, int>(shortFalls[i], i);
	//std::cout << "unsorted" << std::endl;
	//for (int i = 0; i < coupledVectors.size(); i++)
	//	std::cout << coupledVectors[i].first << "   " << coupledVectors[i].second << std::endl;
	std::sort(coupledVectors.begin(), coupledVectors.end(), compareFunc);
	//std::cout << "sorted" << std::endl;
	//for (int i = 0; i < coupledVectors.size(); i++)
	//	std::cout << coupledVectors[i].first << "   " << coupledVectors[i].second << std::endl;

	for (int i = 0; i < std::round(shortFallSum); i++)
		roundedStoicCoeff[coupledVectors[coupledVectors.size() - 1 - i].second] = std::ceil(roundedStoicCoeff[coupledVectors[coupledVectors.size() - 1 - i].second]);
	for (int i = std::round(shortFallSum); i < roundedStoicCoeff.size(); i++)
		roundedStoicCoeff[coupledVectors[coupledVectors.size() - 1 - i].second] = std::floor(roundedStoicCoeff[coupledVectors[coupledVectors.size() - 1 - i].second]);

	for (int i = 0; i < roundedStoicCoeff.size(); i++)
		roundedStoicCoeff[i] = roundedStoicCoeff[i] * digitNumber;

	//std::cout << "start" << std::endl;
	//for (int i = 0; i < stoicCoeff.size(); i++)
	//{
	//	std::cout << std::setw(20) << stoicCoeff[i] << std::setw(20) << roundedStoicCoeff[i] << std::setw(5) << sampleProducts[i].numberOfC() <<
	//		std::setw(5) << sampleProducts[i].numberOfH()<< std::setw(5) << sampleProducts[i].numberOfO() << std::endl;
	//}
	//std::cout << "end" << std::endl;

	for (int i = 0; i < prodNames.size(); i++)
	{
		char buff[10];
		sprintf_s(buff, 10, "%6.4f", roundedStoicCoeff[i]);
		if (std::string(buff) == "0.0000")
			continue;
		if (firstPrinted)
			out << "+";
		if (std::string(buff) != "1.0000")
			out << buff;
		firstPrinted = true;
		out << prodNames[i];
	}
	{
		char buff[100];
		sprintf_s(buff, 100, " %9.3e %8.3f %8.0f\n", As[0], ns[0], Es[0]);
		out << buff;
	}
	if (Pressures.size() > 1)
	{
		char buff[100];
		for (int i = 0; i < Pressures.size(); i++)
		{
			sprintf_s(buff, 100, "PLOG / %7.3f       %9.3e %8.3f %8.0f  /\n",
				Pressures[i], As[i], ns[i], Es[i]);
			out << buff;
		}
	}
	return out.str();
}


std::string LumpedReaction::printTargetRates()
{
	std::stringstream outSS;
	outSS << print();
	for (int j = 0; j < Pressures.size(); j++)
	{
		outSS << "Pressure " << Pressures[j] << " atm:" << std::endl;
		for (int i = 0; i < targetTs.size(); i++)
			outSS << "   " << std::setw(20) << targetTs[i] << std::setw(20) << 
			targetRates[j][i] << std::endl;
		outSS << std::endl;
	}
	return outSS.str();
}

void LumpedReaction::computeFirstGuess(double x[], std::vector<double> temps, std::vector<double> ks)
{
	if (temps.size() < 3)
		UTL::fatalError("LumpedReaction::computeFirstGuess called with less than 3 target points");
	if (temps.size() != ks.size())
		UTL::fatalError("LumpedReaction::computeFirstGuess called with temp vector of different size than ks vector");

	int middleIndex = std::floor(temps.size()/2);
	double T1 = temps[0];
	double T2 = temps[middleIndex];
	double T3 = temps[temps.size() - 1];
	double lnk1 = std::log(ks[0]);
	double lnk2 = std::log(ks[middleIndex]);
	double lnk3 = std::log(ks[ks.size()-1]);
	// debug ->
	//std::cout << "T1 = " << T1 << ", T2 = " << T2 << ", T3 = " << T3 << std::endl;
	//std::cout << "lnK1 = " << lnk1 << ", lnK2 = " << lnk2 << ", lnK3 = " << lnk3 << std::endl;
	// <- debug
	double R = 1.987;
	double alpha = (1.0 / R) * (1 / T2 - 1 / T1) * (std::log(T3) - std::log(T1)) /
		(std::log(T2) - std::log(T1)) + (1 / R / T1 - 1 / R / T3);
	double X3 = (lnk3 - lnk1 - (std::log(T3) - std::log(T1)) /
		(std::log(T2) - std::log(T1)) * (lnk2 - lnk1)) / alpha;
	double X2 = (lnk2 - lnk1 + X3 / R * (1 / T2 - 1 / T1)) / (std::log(T2) - std::log(T1));
	double X1 = lnk1 - X2 * std::log(T1) + X3 / R / T1;

	// debug ->
	//std::cout << "x1 = " << X1 << ", x2 = " << X2 << ", x3 = " << X3 << std::endl;
	// <- debug
	x[0] = X1;
	x[1] = X2;
	x[2] = X3;
}