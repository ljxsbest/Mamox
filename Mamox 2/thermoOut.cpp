#include "thermoOut.h"

double optimFunc(unsigned n, const double* x, double* grad, void* my_func_data)
{
	// x[0] = A_L; x[1] = n_L; x[3] = Eatt_L
	if (grad) {
		std::cerr << "Algorithm for optimization requires gradient!" << std::endl;
	}

	
	std::vector<double>* data = static_cast<std::vector<double>*>(my_func_data);
	if (data->size() % 2 != 0)
		std::cerr << "ERROR: In thermoOut optimFunction, data vector has odd number of elements!" << std::endl;
	int size = data->size() / 2;
	std::vector<double> Cp(size);
	std::vector<double> T(size);
	for (int i = 0; i < size; i++)
	{
		Cp[i] = (*data)[i];
		T[i]  = (*data)[i+size];
	}


	double R_ = 1.987;

	double func = 0.0;

	for (int i = 0; i < size; i++)
	{
		func += std::pow((Cp[i] / R_ - (x[0] +x[1]*T[i] +x[2]*T[i]*T[i] +x[3]*T[i]*T[i]*T[i]+ x[4]*T[i]*T[i]*T[i]*T[i]))/(Cp[i] / R_), 2);
	}

	return func;
}

ThermoOut::ThermoOut(std::string groupsPath, std::string radicalCorrectionsPath,
	std::string knownValuesPath, ChemkinOut* chemout_)
{
	chemOut = chemout_;

	groupsFile.open(groupsPath);
	if (!groupsFile)
		std::cerr << "ERROR: ThermoOut was not able to open file containing the group contributions!" << std::endl;
	
	correctionsFile.open(radicalCorrectionsPath);
	if (!correctionsFile)
		std::cerr << "ERROR: ThermoOut was not able to open file containing the radicals corrections!" << std::endl;

	knownFile.open(knownValuesPath);
	if (!knownFile)
		std::cerr << "ERROR: ThermoOut was not able to open file containing the known values!" << std::endl;
	else
	{
		knownFileProvided = true;
		
		std::string line;

		//bool fileStarted = false;
		//while (std::getline(knownFile, line) && !fileStarted) // search the beginning of the file
		//{
		//	if (line.substr(0, 6) == "THERMO")
		//	{
		//		fileStarted = true;
		//		std::getline(knownFile, line);	//flush the line with the temperatures
		//	}
		//}
		//
		//if (!fileStarted)
		//{
		//	std::cerr << "WARNING: the known values thermo data file is formatted improperly.\nIt has been ingored." << std::endl;
		//	return;
		//}

		while (std::getline(knownFile, line))
		{
			if (line[0] == '!') // if it is a comment skip
				continue;
			if (line.size() < 80)
				continue;
			if (line[79] == '1')
			{
				std::string data;
				std::stringstream ss(line);
				std::string speciesName;
				std::getline(ss, speciesName, ' ');

				data.append(line);
				data.append("\n");
				std::getline(knownFile, line);
				if (line[79] != '2')
					continue;
				data.append(line);
				data.append("\n");
				std::getline(knownFile, line);
				if (line[79] != '3')
					continue;
				data.append(line);
				data.append("\n");
				std::getline(knownFile, line);
				if (line[79] != '4')
					continue;
				data.append(line);
				data.append("\n");

				knownMoleculesNames.push_back(speciesName);
				knownMoleculesData.push_back(data);
			}
		}
	}

	importGroups();
	importCorrections();

}

ThermoOut::ThermoOut(std::string groupsPath, std::string radicalCorrectionsPath, ChemkinOut* chemout_)
{
	chemOut = chemout_;

	groupsFile.open(groupsPath);
	if (!groupsFile)
		std::cerr << "ERROR: ThermoOut was not able to open file containing the group contributions!" << std::endl;

	correctionsFile.open(radicalCorrectionsPath);
	if (!correctionsFile)
		std::cerr << "ERROR: ThermoOut was not able to open file containing the radicals corrections!" << std::endl;

	knownFileProvided = false;

	importGroups();
	importCorrections();
}

void ThermoOut::NASAParameters(Molecola* mol, std::vector<double>* TLimsVec,
	std::vector<double>* paramLowT, std::vector<double>* paramHighT)
{
	Eigen::VectorXd pHighT;
	Eigen::VectorXd pLowT;

	std::vector<double> TLVec;

	int posInKnown = -1;
	if (knownFileProvided)
	{
		// check if the molecule is in the list of known molecules, if yes return the known values
		posInKnown = isInKnownList(mol);
	}
	if (posInKnown != -1)
	{
		knownMoleculeData(mol, &pLowT, &pHighT, &TLVec);
	}
	else
	{
		// if the molecule is not in the known molecules list compute the properties
		Molecola molToAnalyze;

		// if mol is a radical analyze its base molecule, otherwise analyze the mol itsef
		if (mol->posCrad() != 0 || mol->posCOOrad() != 0 || mol->posCOrad() != 0)
		{
			molToAnalyze = mol->noRadicalMolecule();
			//std::cout << *mol << "    =>    " << molToAnalyze << std::endl;
		}
		else
			molToAnalyze = *mol;


		//std::cout << *mol << "   =>   " << molToAnalyze << std::endl;

		std::vector<std::string> groups;
		std::vector<int> frequencies;
		//std::cout << *mol << "    =>    " << molToAnalyze << std::endl;
		for (int i = 1; i < molToAnalyze.size() + 1; i++)
		{
			std::string groupName;
			if (molToAnalyze.tipo(i) != 8)
			{
				if (molToAnalyze.isole(i) == false && molToAnalyze.ischeto(i) == false)
					groupName.append("C");
				else if (molToAnalyze.isole(i) == true && molToAnalyze.isKetene(i) == 0)
					groupName.append("CD");
				else if (molToAnalyze.ischeto(i) == true && molToAnalyze.isKetene(i) == 0)
					groupName.append("CO");
				else if (molToAnalyze.isKetene(i) == 2)
				{
					groupName.append("KETENE");					// NB need to be corrected with the proper ketene gorup!
					addGroupToList(groupName, &groups, &frequencies);
					continue;
				}
				else if (molToAnalyze.isKetene(i) == 1)
					continue;

				int SbondedC = molToAnalyze.numSingleBondedCarbons(i);
				if (SbondedC > 0)
				{
					groupName.append("/C");
					if (SbondedC > 1)
						groupName.append(std::to_string(SbondedC));
				}
				int DbondedC = molToAnalyze.numDoubleBondedCarbons(i);
				if (DbondedC > 0 && molToAnalyze.isole(i) == false)
				{
					groupName.append("/CD");
					if (DbondedC > 1)
						groupName.append(std::to_string(DbondedC));
				}
				int bondedKeto = molToAnalyze.numBondedKetoCarbons(i);
				if (bondedKeto > 0)
				{
					groupName.append("/CO");
					if (bondedKeto > 1)
						groupName.append(std::to_string(bondedKeto));
				}
				int bondedH = molToAnalyze.numAbstractableH(i);
				if (bondedH > 0)
				{
					groupName.append("/H");
					if (bondedH > 1)
						groupName.append(std::to_string(bondedH));
				}
				if (molToAnalyze.tipo(i) == 4)
					groupName.append("/OO");
				if (molToAnalyze.isetero(i))
					groupName.append("/O");
			}
			if (groupName == "CO/C/CD")
				groupName = "CO/CD/C";
			if (groupName == "C/C2/CO/OO")
				groupName = "C/C2/CO/O";
			if (groupName == "C/C2/CD/OO")
				groupName = "C/CD/C2/O";

			addGroupToList(groupName, &groups, &frequencies);
		}
		for (int i = 0; i < molToAnalyze.numCOOH(); i++)
			addGroupToList("OO/C/H", &groups, &frequencies);
		if (molToAnalyze.numberEthero() == 1)
		{
			addGroupToList("O/C2", &groups, &frequencies);
			std::string groupName = "CY/C";
			groupName.append(std::to_string(molToAnalyze.dist(molToAnalyze.posEthero()[0],
				molToAnalyze.posEthero()[1]) + 1));
			groupName.append("O");
			addGroupToList(groupName, &groups, &frequencies);
		}


		double Hf = 0.0;
		double S = 0.0;
		double Cp300 = 0.0;
		double Cp400 = 0.0;
		double Cp500 = 0.0;
		double Cp600 = 0.0;
		double Cp800 = 0.0;
		double Cp1000 = 0.0;
		double Cp1500 = 0.0;

		for (int i = 0; i < groups.size(); i++)			// TODO: consider the possibility of other units of measure
		{
			Group currGroup = findGroup(groups[i], &groupList);
			Hf += frequencies[i] * currGroup.Hf * 1000.0;			// [cal/mol]
			S += frequencies[i] * currGroup.S;					// [cal/mol/K]
			Cp300 += frequencies[i] * currGroup.Cp300;				// [cal/mol/K]
			Cp400 += frequencies[i] * currGroup.Cp400;				// [cal/mol/K]
			Cp500 += frequencies[i] * currGroup.Cp500;				// [cal/mol/K]
			Cp600 += frequencies[i] * currGroup.Cp600;				// [cal/mol/K]
			Cp800 += frequencies[i] * currGroup.Cp800;				// [cal/mol/K]
			Cp1000 += frequencies[i] * currGroup.Cp1000;				// [cal/mol/K]
			Cp1500 += frequencies[i] * currGroup.Cp1500;				// [cal/mol/K]
		}

		// reference http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/thermo.html#yu

		// symmetry correction
		S -= 1.987 * log(double(molToAnalyze.numberOfSymmetries()));

		// DEBUG
		//if(chemOut->molToName(*mol) == "ISOC16")
		//	std::cout << "num sym = " << molToAnalyze.numberOfSymmetries() << ",  corr:"
		//		<< 1.987 * log(double(molToAnalyze.numberOfSymmetries())) << std::endl;


		// chirality correction
		if (molToAnalyze.isChiral())
			S += 1.987 * log(2);

		// add Hydrogen Bond Increment if the molecule is a radical

		std::vector<std::string> radGroups;
		std::vector<int> radFrequencies;
		if (mol->posCrad() != 0 && mol->posCOOrad() == 0 && mol->posCOrad() == 0)	// if the molecule is a carbon centered radical
		{
			//Rpmet, Rpet, Rp, Rs, Rt
			switch (mol->tipoR(mol->posCrad()))
			{
			case Rpmet:
				addGroupToList("P", &radGroups, &radFrequencies);
				break;
			case Rpet:
				addGroupToList("P", &radGroups, &radFrequencies);
				break;
			case Rp:
				if (mol->isCradAllylic(mol->posCrad()))
					addGroupToList("ALLYLP", &radGroups, &radFrequencies);
				else
					addGroupToList("P", &radGroups, &radFrequencies);
				break;
			case Rs:
				if (mol->isCradAllylic(mol->posCrad()))
					addGroupToList("ALLYLS", &radGroups, &radFrequencies);
				else
					addGroupToList("S", &radGroups, &radFrequencies);
				break;
			case Rt:
				if (mol->isCradAllylic(mol->posCrad()))
					addGroupToList("ALLYLT", &radGroups, &radFrequencies);
				else
					addGroupToList("T", &radGroups, &radFrequencies);
				break;
			}

		}
		else if (mol->posCrad() == 0 && mol->posCOOrad() != 0 && mol->posCOrad() == 0)
		{
			addGroupToList("ALPEROX", &radGroups, &radFrequencies);
		}
		else if (mol->posCrad() == 0 && mol->posCOOrad() == 0 && mol->posCOrad() != 0)
		{
			addGroupToList("ALKOXY", &radGroups, &radFrequencies);
		}


		for (int i = 0; i < radGroups.size(); i++)			// TODO: consider the possibility of other units of measure
		{
			Group currGroup = findGroup(radGroups[i], &correctionsList);
			Hf += radFrequencies[i] * currGroup.Hf * 1000.0 - 52100;	// [cal/mol/K]
			S += radFrequencies[i] * currGroup.S;						// [cal/mol/K]
			Cp300 += radFrequencies[i] * currGroup.Cp300;				// [cal/mol/K]
			Cp400 += radFrequencies[i] * currGroup.Cp400;				// [cal/mol/K]
			Cp500 += radFrequencies[i] * currGroup.Cp500;				// [cal/mol/K]
			Cp600 += radFrequencies[i] * currGroup.Cp600;				// [cal/mol/K]
			Cp800 += radFrequencies[i] * currGroup.Cp800;				// [cal/mol/K]
			Cp1000 += radFrequencies[i] * currGroup.Cp1000;				// [cal/mol/K]
			Cp1500 += radFrequencies[i] * currGroup.Cp1500;				// [cal/mol/K]
		}

		////DEBUG -START
		//if (chemOut->molToName(*mol) == "I2C8KET1-2")
		//{
		//	std::cout << "DEBUG" << std::endl;
		//	std::cout << Cp300 << "  " << Cp400 << "  " << Cp500 << "  " << Cp600 << "  " <<
		//		Cp800 << "  " << Cp1000 << "  " << Cp1500 << std::endl;
		//	std::cout << Hf << "   " << S << std::endl;
		//
		//}
		////DEBUG - END

		Eigen::VectorXd Temps(7);
		Eigen::VectorXd CpsVec(7);

		// fit the parameters      TODO: make it automatically detect the temperatures if values at different temperatures are provided
		Temps << 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0;
		CpsVec << Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500;
		fitNASAParam(Temps, CpsVec, Hf, 298.0, S, 298.0, &pLowT);
		pHighT = pLowT;
		TLVec = { 300.0, 800.0, 1500.0 };
	}
	*TLimsVec = TLVec;
	paramHighT->resize(pHighT.size());
	paramLowT->resize(pLowT.size());
	for (int i = 0; i < pHighT.size(); i++)
	{
		(*paramHighT)[i] = pHighT(i);
		(*paramLowT)[i] = pLowT(i);
	}

}

std::string ThermoOut::NASAOutput(Molecola* mol, std::string name)
{
	
	int numC = mol->numberOfC();
	int numH = mol->numberOfH();
	int numO = mol->numberOfO();

	std::vector<double> TLimsVec;
	std::vector<double> paramLowT;
	std::vector<double> paramHighT;

	NASAParameters(mol, &TLimsVec, &paramLowT, &paramHighT);
	

	return NASAParametersToString(name, numC, numH, numO, TLimsVec, paramLowT, paramHighT);

}

std::string ThermoOut::NASAOutput(Molecola* mol)
{
	std::string name = chemOut->molToName(*mol);
	return NASAOutput(mol, name);
}

std::string ThermoOut::NASAParametersToString(std::string name, int numC, int numH, int numO, std::vector<double> Tlims, Eigen::VectorXd paramLowT, Eigen::VectorXd paramHighT)
{
	std::stringstream output;

	output << std::left << std::setw(24) << name;

	if (numC > 0)
		output << "C" << std::right << std::setw(4) << numC;
	else
		output << "     ";

	if (numH > 0)
		output << "H" << std::right << std::setw(4) << numH;
	else
		output << "     ";

	if (numO > 0)
		output << "O" << std::right << std::setw(4) << numO;
	else
		output << "     ";

	char buffer[120];
	//output << "     G    300.00   1500.00  800.00      1" << std::endl;
	sprintf_s(buffer, 120, "     G   %7.2f  %7.2f  %7.2f      1", Tlims[0], Tlims[2], Tlims[1]);
	output << buffer << std::endl;
	sprintf_s(buffer, 120, "%+.8E%+.8E%+.8E%+.8E%+.8E    2", paramHighT(0), paramHighT(1), paramHighT(2), paramHighT(3), paramHighT(4));
	output << buffer << std::endl;
	sprintf_s(buffer, 120, "%+.8E%+.8E%+.8E%+.8E%+.8E    3", paramHighT(5), paramHighT(6), paramLowT(0), paramLowT(1), paramLowT(2));
	output << buffer << std::endl;
	sprintf_s(buffer, 120, "%+.8E%+.8E%+.8E%+.8E                   4", paramLowT(3), paramLowT(4), paramLowT(5), paramLowT(6));
	output << buffer << std::endl;

	return output.str();
}

std::string ThermoOut::NASAParametersToString(std::string name, int numC, int numH, int numO, std::vector<double> Tlims, std::vector<double> lowTParams, std::vector<double> highTParams)
{
	//std::cout << "nasa param" << std::endl;
	Eigen::VectorXd lowTParamsEigen(lowTParams.size());
	for (int i = 0; i < lowTParams.size(); i++)
		lowTParamsEigen(i) = lowTParams[i];

	Eigen::VectorXd highTParamsEigen(highTParams.size());
	for (int i = 0; i < highTParams.size(); i++)
		highTParamsEigen(i) = highTParams[i];

	return NASAParametersToString(name, numC, numH, numO, Tlims, lowTParamsEigen, highTParamsEigen);
}

void ThermoOut::fitNASAParam(Eigen::VectorXd temps, Eigen::VectorXd Cps, double Hf, double T_Hf, double S, double T_S, Eigen::VectorXd* params)
{
	if (temps.size() != Cps.size())
	{
		std::cerr << "ERROR: in ThermoOut::fitCpParam temp and Cps are not the same size!" << std::endl;
		return;
	}
	int size = 5;
	
	Eigen::VectorXd Cps5(size, 1);
	Eigen::MatrixXd matrix(size, size);
	for (int i = 0; i < size; i++)
	{
		Cps5(i) = Cps(i)/1.987;
		for (int j = 0; j < size; j++)
		{
			matrix(i, j) = std::pow(temps(i), j);
		}
	}
	Eigen::MatrixXd invMatr = matrix.inverse();

	Eigen::VectorXd CpParamFG = invMatr * Cps5;	// compute first guess parameters 

	nlopt_opt opt_loc;
	opt_loc = nlopt_create(NLOPT_LN_BOBYQA, 5);
	std::vector<double> inData(Cps.size() * 2);
	for (int i = 0; i < Cps.size(); i++)
	{
		inData[i] = Cps[i];
		inData[i+Cps.size()] = temps[i];
	}
	//std::cout << "debug indata ";
	//for (int i = 0; i < inData.size(); i++)
	//	std::cout << inData[i] << "  ";
	//std::cout << std::endl;

	nlopt_set_min_objective(opt_loc, optimFunc, &inData);

	nlopt_set_xtol_rel(opt_loc, 1.0e-15);
	nlopt_set_maxtime(opt_loc, 10.);

	double CpParam[5];
	CpParam[0] = CpParamFG[0];	// set first guess parameters
	CpParam[1] = CpParamFG[1];
	CpParam[2] = CpParamFG[2];
	CpParam[3] = CpParamFG[3];
	CpParam[4] = CpParamFG[4];
	double minf = 0.;
	int exitFlag = nlopt_optimize(opt_loc, CpParam, &minf);
	//std::cout << "DEBUG: minf = " << minf << std::endl;
	
	double par5 = T_Hf * (Hf/1.987/ T_Hf - CpParam[0] - T_Hf *CpParam[1]/2.0 
		- std::pow(T_Hf, 2) * CpParam[2] / 3.0 - std::pow(T_Hf, 3) * CpParam[3] / 4.0 
		- std::pow(T_Hf, 4) * CpParam[4] / 5.0);

	double par6 = S / 1.987 - CpParam[0] * std::log(T_S) - CpParam[1] * T_S
		- CpParam[2] / 2.0 * std::pow(T_S, 2) - CpParam[3] / 3.0 * std::pow(T_S, 3)
		- CpParam[4] / 4.0 * std::pow(T_S, 4);
	
	params->resize(7);

	for (int i = 0; i < 5; i++)
		(*params)(i) = CpParam[i];
	
	(*params)(5) = par5;
	(*params)(6) = par6;
	return;
}

int ThermoOut::isInKnownList(Molecola* mol)
{
	return isInKnownList(chemOut->molToName(*mol));
}

int ThermoOut::isInKnownList(std::string molName)
{
	for (int i = 0; i < knownMoleculesNames.size(); i++)
	{
		if (molName == knownMoleculesNames[i])
		{
			return i;
		}
	}
	return -1;
}

void ThermoOut::importGroups()
{
	std::string line;
	std::getline(groupsFile, line);		//discard first line
	while (!groupsFile.eof())
	{
		std::getline(groupsFile, line);
		std::stringstream lineSS(line);

		Group currentGroup;

		std::getline(lineSS, line, ',');		// get group name
		currentGroup.name = line;
		std::getline(lineSS, line, ',');		// get group units of measure
		currentGroup.units = line;
		std::getline(lineSS, line, ',');		// get group heat of formation
		if (line != "")
			currentGroup.Hf = std::stod(line);
		else
			currentGroup.Hf = 0;
		std::getline(lineSS, line, ',');		// get group entropy
		if (line != "")
			currentGroup.S = std::stod(line);
		else
			currentGroup.S = 0;
		std::getline(lineSS, line, ',');		// get group Cps
		if (line != "")
			currentGroup.Cp300 = std::stod(line);
		else
			currentGroup.Cp300 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp400 = std::stod(line);
		else
			currentGroup.Cp400 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp500 = std::stod(line);
		else
			currentGroup.Cp500 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp600 = std::stod(line);
		else
			currentGroup.Cp600 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp800 = std::stod(line);
		else
			currentGroup.Cp800 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp1000 = std::stod(line);
		else
			currentGroup.Cp1000 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp1500 = std::stod(line);
		else
			currentGroup.Cp1500 = 0;

		groupList.push_back(currentGroup);
	}
}

void ThermoOut::importCorrections()
{
	std::string line;
	std::getline(correctionsFile, line);		//discard first line
	while (!correctionsFile.eof())
	{
		std::getline(correctionsFile, line);
		std::stringstream lineSS(line);

		Group currentGroup;

		std::getline(lineSS, line, ',');		// get group name
		currentGroup.name = line;
		std::getline(lineSS, line, ',');		// get group units of measure
		currentGroup.units = line;
		std::getline(lineSS, line, ',');		// get group heat of formation
		if (line != "")
			currentGroup.Hf = std::stod(line);
		else
			currentGroup.Hf = 0;
		std::getline(lineSS, line, ',');		// get group entropy
		if (line != "")
			currentGroup.S = std::stod(line);
		else
			currentGroup.S = 0;
		std::getline(lineSS, line, ',');		// get group Cps
		if (line != "")
			currentGroup.Cp300 = std::stod(line);
		else
			currentGroup.Cp300 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp400 = std::stod(line);
		else
			currentGroup.Cp400 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp500 = std::stod(line);
		else
			currentGroup.Cp500 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp600 = std::stod(line);
		else
			currentGroup.Cp600 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp800 = std::stod(line);
		else
			currentGroup.Cp800 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp1000 = std::stod(line);
		else
			currentGroup.Cp1000 = 0;
		std::getline(lineSS, line, ',');
		if (line != "")
			currentGroup.Cp1500 = std::stod(line);
		else
			currentGroup.Cp1500 = 0;

		correctionsList.push_back(currentGroup);
	}
}

void ThermoOut::addGroupToList(std::string groupName, std::vector<std::string>* groupVec,
	std::vector<int>* freqVec)
{
	for (int i = 0; i < groupVec->size(); i++)
	{
		if (groupName == (*groupVec)[i])	// if present increase
		{
			(*freqVec)[i]++;
			return;
		}
	}
	// if not present add it
	groupVec->push_back(groupName);
	freqVec->push_back(1);

	return;
}

Group ThermoOut::findGroup(std::string groupName, std::vector<Group>* groupVector)
{
	for (int i = 0; i < groupVector->size(); i++)
	{
		if (groupName == (*groupVector)[i].name)
			return (*groupVector)[i];
	}

	// JIAXIN std::cout << "WARNING: in ThermoOut::findGroup " << groupName << " was not found!" << std::endl;
	Group errorGroup;
	errorGroup.name = "ERR";
	errorGroup.Hf = 0.0;
	errorGroup.S = 0.0;
	errorGroup.Cp300 = 0.0;
	errorGroup.Cp400 = 0.0;
	errorGroup.Cp500 = 0.0;
	errorGroup.Cp600 = 0.0;
	errorGroup.Cp800 = 0.0;
	errorGroup.Cp1000 = 0.0;
	errorGroup.Cp1500 = 0.0;
	return  errorGroup;
}

double ThermoOut::Hf_T(double T, Eigen::VectorXd params)	// return the heat of formation at the temperature T in [cal/mol]
{
	return 1.987 * T * (params(0) + params(1) * std::pow(T, 1) / 2 + params(2) * std::pow(T, 2) / 3 
		+ params(3) * std::pow(T, 3) / 4 + params(4) * std::pow(T, 4) / 5 + params(5) / T);
}

double ThermoOut::S_T(double T, Eigen::VectorXd params)	// return the entropy at the temperature T in [cal/mol/K]
{
	return 1.987 * (params(0) * std::log(T) + params(1) * T + params(2) * std::pow(T, 2) / 2 
		+ params(3) * std::pow(T, 3) / 3 + params(4) * std::pow(T, 4) / 4 + params(6));
}

double ThermoOut::Cp_T(double T, Eigen::VectorXd params)	// return the specific heat at the temperature T in [cal/mol/K]
{
	return 1.987 * (params(0) + params(1) * T + params(2) * std::pow(T, 2)
		+ params(3) * std::pow(T, 3) + params(4) * std::pow(T, 4));
}

double ThermoOut::Hf_T(double T, std::vector<double> params)	// return the heat of formation at the temperature T in [cal/mol]
{
	Eigen::VectorXd paramEigen(params.size());
	for (int i = 0; i < params.size(); i++)
		paramEigen(i) = params[i];
	return Hf_T(T, paramEigen);
}

double ThermoOut::S_T(double T, std::vector<double> params)		// return the entropy at the temperature T in [cal/mol/K]
{
	Eigen::VectorXd paramEigen(params.size());
	for (int i = 0; i < params.size(); i++)
		paramEigen(i) = params[i];
	return S_T(T, paramEigen);
}

double ThermoOut::Cp_T(double T, std::vector<double> params)	// return the specific heat at the temperature T in [cal/mol/K]	
{
	Eigen::VectorXd paramEigen(params.size());
	for (int i = 0; i < params.size(); i++)
		paramEigen(i) = params[i];
	return Cp_T(T, paramEigen);
}


std::string ThermoOut::NASAOutputLumped(std::string name, std::vector<Molecola> mols,
	std::vector<double> temperatures, std::vector<std::vector<double>> distr)
{
	if (mols.size() == 1)
	{
		//std::cout << name << std::endl;
		return NASAOutput(&(mols[0]), name);
	}
	
	std::vector<std::vector<double>> paramsLowT(mols.size());
	std::vector<std::vector<double>> paramsHighT(mols.size());
	std::vector<std::vector<double>> tempsLimits(mols.size());
	
	for (int i = 0; i < mols.size(); i++)
		NASAParameters(&(mols[i]), &(tempsLimits[i]), &(paramsLowT[i]), &(paramsHighT[i]));
		
	double lowTLim = 0;
	double midTLim = 0;
	double highTLim = 1.0e20;

	// find low limit
	for (int i = 0; i < tempsLimits.size(); i++)
		if (tempsLimits[i][0] > lowTLim)
			lowTLim = tempsLimits[i][0];
	// find high T limit
	for (int i = 0; i < tempsLimits.size(); i++)
		if (tempsLimits[i][2] < highTLim)
			highTLim = tempsLimits[i][2];
	// find midT
	for(int i = 0; i < tempsLimits.size(); i++)
		midTLim += tempsLimits[i][1];
	midTLim = midTLim / tempsLimits.size();

	if (lowTLim < temperatures[0])
		lowTLim = temperatures[0];
	if (highTLim > temperatures[temperatures.size() - 1])
		highTLim = temperatures[temperatures.size() - 1];
	if (midTLim < lowTLim || midTLim > highTLim)
		midTLim = (lowTLim + highTLim) / 2.0;

	// compute the Cp points in the low temperature range
	std::vector<double> TsLowT(7);
	for (int i = 0; i < 7; i++)
		TsLowT[i] = lowTLim + (midTLim - lowTLim) / 6 * i;
	
	std::vector<double> CpsLowT(7, 0.0);
	std::vector<std::vector<double>> distrLowT(mols.size());
	for (int i = 0; i < mols.size(); i++)
	{
		distrLowT[i] = UTL::linInterpolation(temperatures, distr[i], TsLowT);
	}
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < mols.size(); j++)
		{
			if(TsLowT[i] <= tempsLimits[j][1])
				CpsLowT[i] += Cp_T(TsLowT[i], paramsLowT[j]) * distrLowT[j][i];
			else
				CpsLowT[i] += Cp_T(TsLowT[i], paramsHighT[j]) * distrLowT[j][i];
		}
	}

	// compute the Cp points in the high temperature range
	std::vector<double> TsHighT(7);
	for (int i = 0; i < 7; i++)
		TsHighT[i] = midTLim + (highTLim - midTLim) / 6 * i;

	std::vector<double> CpsHighT(7, 0.0);
	std::vector<std::vector<double>> distrHighT(mols.size());
	for (int i = 0; i < mols.size(); i++)
	{
		distrHighT[i] = UTL::linInterpolation(temperatures, distr[i], TsHighT);
	}
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < mols.size(); j++)
		{
			if (TsHighT[i] <= tempsLimits[j][1])
				CpsHighT[i] += Cp_T(TsHighT[i], paramsLowT[j]) * distrHighT[j][i];
			else
				CpsHighT[i] += Cp_T(TsHighT[i], paramsHighT[j]) * distrHighT[j][i];
		}
	}

	// compute Hf for mid T range
	double HfMidT = 0.0;
	for (int i = 0; i < mols.size(); i++)
	{
		if (midTLim <= tempsLimits[i][1])
			HfMidT += Hf_T(midTLim, paramsLowT[i]) * distrHighT[i][0];
		else
			HfMidT += Hf_T(midTLim, paramsHighT[i]) * distrHighT[i][0];
	}

	// compute S for mid T range				// TODO: check, maybe it is better to use the mixing rule of mixtures (cosidering entropy of mixing)
	double SMidT = 0.0;
	for (int i = 0; i < mols.size(); i++)
	{
		if (midTLim <= tempsLimits[i][1])
			SMidT += S_T(midTLim, paramsLowT[i]) * distrHighT[i][0];
		else
			SMidT += S_T(midTLim, paramsHighT[i]) * distrHighT[i][0];
	}
	//std::cout << "Hf = " << HfMidT << ",     S = " << SMidT << std::endl;
	std::vector<double> lowTLumpParam;
	std::vector<double> highTLumpParam;
	fitNASAParam(TsLowT, CpsLowT, HfMidT, midTLim, SMidT, midTLim, &lowTLumpParam);
	fitNASAParam(TsHighT, CpsHighT, HfMidT, midTLim, SMidT, midTLim, &highTLumpParam);

	std::vector<double> TLimVec = { lowTLim, midTLim, highTLim };
	return NASAParametersToString(name, mols[0].numberOfC(), mols[0].numberOfH(), 
		mols[0].numberOfO(), TLimVec, lowTLumpParam, highTLumpParam);
}

bool isValidNumber(std::string str)
{
	bool isThereADigit = false;
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] != ' ')
		{
			if (std::isdigit(str[i]))
				isThereADigit = true;
			else
				return false;
		}
	}
	if (isThereADigit)
		return true;
	else
		return false;
}


void ThermoOut::knownMoleculeData(std::string molName, std::vector<double>* parLowT,
	std::vector<double>* parHighT, std::vector<double>* Temps)
{
	Temps->resize(3);
	parLowT->resize(7);
	parHighT->resize(7);

	int index = isInKnownList(molName);
	if (index == -1)
	{
		std::cerr << "WARNING: knownMoleculeData called on a non known molecule ("<< molName << ")" << std::endl;
		return;
	}

	std::string line;
	std::string buffer;
	std::stringstream dataSS(knownMoleculesData[index]);

	std::getline(dataSS, line);

	std::stringstream ss(line);
	species spec;

	// read species name
	std::getline(ss, buffer, ' ');
	//count the number of atoms
	int numAt = 0;
	std::vector<int> posNumbers = { 25, 30, 35 };
	int numC, numH, numO;
	for (int i = 0; i < 3; i++)
	{
		std::string numStr = line.substr(posNumbers[i], 4);
		if (isValidNumber(numStr))
		{
			numAt = std::stoi(numStr);
			if (numAt != 0)
			{
				char atomLetter = line[posNumbers[i] - 1];
				switch (atomLetter)
				{
				case 'C':
					numC = numAt;
					break;
				case 'H':
					numH = numAt;
					break;
				case 'O':
					numO = numAt;
					break;
				default:
					std::cerr << "WARNING: atom not recognized while reading species " << molName << "." << std::endl;
					break;
				}
			}
		}
	}

	// get the temperatures limit of the polinomials
	std::string temperaturesString = line.substr(45, 29);
	std::stringstream temperaturesStrStream(temperaturesString);
	std::string temp;

	bool tempFound = false;
	while (tempFound == false)
	{
		std::getline(temperaturesStrStream, temp, ' ');
		if (temp != "")
			tempFound = true;
	}
	(*Temps)[0] = std::stod(temp);

	tempFound = false;
	while (tempFound == false)
	{
		std::getline(temperaturesStrStream, temp, ' ');
		if (temp != "")
			tempFound = true;
	}
	(*Temps)[2] = std::stod(temp);

	tempFound = false;
	while (tempFound == false)
	{
		std::getline(temperaturesStrStream, temp, ' ');
		if (temp != "")
			tempFound = true;
	}
	(*Temps)[1] = std::stod(temp);

	// get the parameters of the polynomials
	
	std::getline(dataSS, line);
	std::string paramStr;
	double paramDoub;
	paramStr = line.substr(0, 15);
	paramDoub = std::stod(paramStr);
	(*parHighT)[0] = paramDoub;

	paramStr = line.substr(15, 15);
	paramDoub = std::stod(paramStr);
	(*parHighT)[1] = paramDoub;

	paramStr = line.substr(30, 15);
	paramDoub = std::stod(paramStr);
	(*parHighT)[2] = paramDoub;

	paramStr = line.substr(45, 15);
	paramDoub = std::stod(paramStr);
	(*parHighT)[3] = paramDoub;

	paramStr = line.substr(60, 15);
	paramDoub = std::stod(paramStr);
	(*parHighT)[4] = paramDoub;

	std::getline(dataSS, line);

	paramStr = line.substr(0, 15);
	paramDoub = std::stod(paramStr);
	(*parHighT)[5] = paramDoub;

	paramStr = line.substr(15, 15);
	paramDoub = std::stod(paramStr);
	(*parHighT)[6] = paramDoub;

	paramStr = line.substr(30, 15);
	paramDoub = std::stod(paramStr);
	(*parLowT)[0] = paramDoub;

	paramStr = line.substr(45, 15);
	paramDoub = std::stod(paramStr);
	(*parLowT)[1] = paramDoub;

	paramStr = line.substr(60, 15);
	paramDoub = std::stod(paramStr);
	(*parLowT)[2] = paramDoub;

	std::getline(dataSS, line);

	paramStr = line.substr(0, 15);
	paramDoub = std::stod(paramStr);
	(*parLowT)[3] = paramDoub;

	paramStr = line.substr(15, 15);
	paramDoub = std::stod(paramStr);
	(*parLowT)[4] = paramDoub;

	paramStr = line.substr(30, 15);
	paramDoub = std::stod(paramStr);
	(*parLowT)[5] = paramDoub;

	paramStr = line.substr(45, 15);
	paramDoub = std::stod(paramStr);
	(*parLowT)[6] = paramDoub;

}

void ThermoOut::knownMoleculeData(Molecola* mol, std::vector<double>* parLowT,
	std::vector<double>* parHighT, std::vector<double>* Temps)
{
	knownMoleculeData(chemOut->molToName(*mol), parLowT, parHighT, Temps);
}

void ThermoOut::fitNASAParam(std::vector<double> temps, std::vector<double> Cps, double Hf, double T_Hf, double S, double T_S, std::vector<double>* params)
{
	Eigen::VectorXd tempsEigen(temps.size());
	for (int i = 0; i < temps.size(); i++)
		tempsEigen(i) = temps[i];

	Eigen::VectorXd CpsEigen(Cps.size());
	for (int i = 0; i < Cps.size(); i++)
		CpsEigen(i) = Cps[i];

	Eigen::VectorXd paramsEigen;

	fitNASAParam(tempsEigen, CpsEigen, Hf, T_Hf, S, T_S, &paramsEigen);

	params->resize(paramsEigen.size());
	for (int i = 0; i < paramsEigen.size(); i++)
		(*params)[i] = paramsEigen(i);
}

void ThermoOut::knownMoleculeData(std::string molName, Eigen::VectorXd* parLowT,
	Eigen::VectorXd* parHighT, std::vector<double>* Temps)
{
	std::vector<double> parLowTSTD;
	std::vector<double> parHighTSTD;

	knownMoleculeData(molName, &parLowTSTD, &parHighTSTD, Temps);

	parLowT->resize(parLowTSTD.size());
	parHighT->resize(parHighTSTD.size());
	for (int i = 0; i < parLowTSTD.size(); i++)
	{
		(*parLowT)(i) = parLowTSTD[i];
		(*parHighT)(i) = parHighTSTD[i];
	}
}

void ThermoOut::knownMoleculeData(Molecola* mol, Eigen::VectorXd* parLowT,
	Eigen::VectorXd* parHighT, std::vector<double>* Temps)
{
	knownMoleculeData(chemOut->molToName(*mol), parLowT, parHighT, Temps);
}


double ThermoOut::Hf(Molecola mol, double T)	// return the heat of formation at the temperature T in [cal/mol]
{

	std::vector<double> TLimsVec;
	std::vector<double> paramLowT;
	std::vector<double> paramHighT;
	NASAParameters(&mol, &TLimsVec, &paramLowT, &paramHighT);
	if(T < TLimsVec[1])
		return Hf_T(T, paramLowT);
	else
		return Hf_T(T, paramHighT);
}

double ThermoOut::S(Molecola mol, double T)	// return the entropy at the temperature T in [cal/mol/K]
{
	std::vector<double> TLimsVec;
	std::vector<double> paramLowT;
	std::vector<double> paramHighT;
	NASAParameters(&mol, &TLimsVec, &paramLowT, &paramHighT);
	if (T < TLimsVec[1])
		return S_T(T, paramLowT);
	else
		return S_T(T, paramHighT);
}
