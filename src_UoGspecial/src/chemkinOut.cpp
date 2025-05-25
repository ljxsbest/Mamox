#include "chemkinOut.h"
#include "XML.h"

int ChemkinOut::Add(Molecola newMol, std::vector<Molecola>* vec)   // add newMol to vec
{
	for (int i = 0; i < vec->size(); i++)
		if (newMol == (*vec)[i]) return i;	// if it's present return it's index

	vec->push_back(newMol);					// else add it and return it's index
	return vec->size() - 1;
};

ChemkinOut::ChemkinOut(std::string fileNameDetailed)
	:outFileDetailed(fileNameDetailed, std::ofstream::out)
{
	XML xmlFile("resources/namesList/molList_pengzhi.xml", 'r');
	std::cout << std::endl << "Loading named species list ..." << std::endl;
	xmlFile.getMolecules(&namedSpecies, &namesList);
	std::cout << "		... named species loaded." << std::endl << std::endl;
}

int ChemkinOut::findPosition(std::vector<int>* vect, int num)		// if num is present in vect return its position, otherwise return -1
{
	for (int i = 0; i < (*vect).size(); i++)
	{
		if ((*vect)[i] == num) return i;
	}
	return -1;
}


bool ChemkinOut::addUniqueName(std::vector<std::string>* vect, std::string name)
{
	if (std::find(vect->begin(), vect->end(), name) == vect->end()) {
		// name not in vect, add it
		vect->push_back(name);
		return true;
	}
	else
	{
		return false;
	}
}

std::vector<int> ChemkinOut::searchDuplicates(std::vector<std::vector<std::string>>* vect, std::vector<std::string> reac)
{
	std::vector<int> indexes;
	for (int i = 0; i < (*vect).size(); i++)
	{
		if ((*vect)[i] == reac) indexes.push_back(i);
	}
	return indexes;
}


void printVector(std::vector<std::string>* vec)
{
	for (int i = 0; i < vec->size(); i++)
		std::cout << (*vec)[i] << std::endl;
	std::cout << std::endl;
}


void ChemkinOut::writeHeadingDetailed(std::stringstream* sstr)
{
	detailedHeading << sstr->rdbuf();
}

void ChemkinOut::writeHeadingLumped(std::stringstream* sstr)
{
	lumpedHeading << sstr->rdbuf();
}

void removeH(std::string* name) // if it is present in the name, removes one hydrogen to a alkane name in order to make it a substitute group
{
	int Hpos = -1;
	for (int i = 0; i < name->size(); i++)
	{
		if ((*name)[i] == 'H')
			Hpos = i;
	}
	if (Hpos == name->size() - 1)
	{
		name->erase(Hpos, 1);
	}
	else if (Hpos == -1)
	{

	}
	else 
	{
		std::string num = name->substr(Hpos + 1, name->size() - 1 - Hpos);
		bool isNumber = true;
		for (int i = 0; i < num.size(); i++)
			if (std::isdigit(num[i]) == false)
				isNumber = false;
		if (isNumber)
		{
			int val = std::stoi(num);
			val = val - 1;
			if (val == 0)
			{
				(*name) = name->substr(0, Hpos);
			}
			else
			{
				(*name) = name->substr(0, Hpos+1);
				name->append(std::to_string(val));
			}
		}
	}
}

std::string ChemkinOut::alkaneStructurePrefix(Molecola* mol)
{
	std::vector<int> mainChain(0);
	std::vector<int> mePos(0);
	std::vector<int> etPos(0);
	std::vector<int> map(0);
	mol->orderedListOfMeAndEt(&mainChain, &mePos, &etPos, &map);
	int numMe = mePos.size();
	int numEt = etPos.size();

	std::string prefix;
	if (mol->isLinear())
		prefix.append("N");
	else
	{
		bool alreadyNamed = false;
		if (numMe == 2 && numEt == 0)
		{
			if ((mePos[0] == mainChain[1] && mePos[1] == mainChain[1])
				|| (mePos[0] == mainChain[mainChain.size() - 2] && mePos[1] == mainChain[mainChain.size() - 2]))
			{
				prefix.append("NE");
				alreadyNamed = true;
			}
			else if ((mePos[0] == mainChain[mainChain.size() - 2] && mePos[1] == mainChain[1])
				|| (mePos[0] == mainChain[1] && mePos[1] == mainChain[mainChain.size() - 2]))
			{
				prefix.append("X");
				alreadyNamed = true;
			}
		}
		if (alreadyNamed == false)
		{
			prefix.append("I");
			for (int i = 0; i < numMe; i++)
				prefix.append(std::to_string(mePos[i]));
			if (numEt > 0)
			{
				for (int i = 0; i < numEt; i++)
					prefix.append(std::to_string(etPos[i]));
				prefix.append("ET");
			}
		}

	}

	// TODO : make the following part of the code more generic, like adding another list of known molecule
	//        and attributing the proper prefix to the known molecules.
	if ((mol->numberOfC() == 16) && (prefix == "I2244688"))
	{
		prefix = "ISO";
	}
	if ((mol->numberOfC() == 15) && (prefix == "I224468"))
	{
		prefix = "ISO";
	}
	if ((mol->numberOfC() == 15) && (prefix == "I224488"))
	{
		prefix = "ISO2";
	}
	if ((mol->numberOfC() == 15) && (prefix == "I224668"))
	{
		prefix = "ISO3";
	}
	if ((mol->numberOfC() == 15) && (prefix == "I224688"))
	{
		prefix = "ISO4";
	}
	if ((mol->numberOfC() == 12) && (prefix == "I22466"))
	{
		prefix = "ISO";
	}

	return prefix;
}

std::string ChemkinOut::molToName(Molecola mol)
{

	std::string name;
	species type = mol.kindOfSPecies();

	if (type == unidentified_)
	{
		std::cerr << "Error in naming the molecule " << mol << ", it doesn't match any expected molecule type." << std::endl;
		mol.printStruttura();
		std::cout << std::endl;
		mol.printDoppioO();
		std::cout << std::endl;
		std::cout  
			<< "numberR     "<< mol.numCrad()			  << std::endl
			<< "numberCOO   "<< mol.numCOORad()		  << std::endl
			<< "numberCOOH  "<< mol.numCOOH()			  << std::endl
			<< "numberEther "<< mol.numberEthero()	  << std::endl
			<< "numberOLE   "<< mol.numberOle()		  << std::endl
			<< "numberCO    "<< mol.numKeto()			  << std::endl
			<< "specMol     "<< mol.isSpecialMolecule() << std::endl
			<< "numLinEth   "<< mol.numLinEther()       << std::endl
			<< "numCOR      "<< mol.numCORad() << std::endl << std::endl;
		exit(0);
	}

	std::string n = checkIfNamed(&mol);
	//n = "not found";
	if (n != "not found")
		name = n;
	else
	{
		if (type == special_)
		{
			switch (mol.isSpecialMolecule())
			{
			case 1:
				name.append("N2");
				break;
			case 2:
				name.append("O2");
				break;
			case 3:
				name.append("O");
				break;
			case 4:
				name.append("OH");
				break;
			case 5:
				name.append("HO2");
				break;
			case 6:
				name.append("H2O");
				break;
			case 7:
				name.append("H2O2");
				break;
			case 8:
				name.append("H");
				break;
			case 9:
				name.append("H2");
				break;
			case 10:
				name.append("CO");
				break;
			}
		}
		else
		{
			// name the base alkane
			std::vector<int> mainChain(0);
			std::vector<int> mePos(0);
			std::vector<int> etPos(0);
			std::vector<int> map(0);
			mol.orderedListOfMeAndEt(&mainChain, &mePos, &etPos, &map);
			int numMe = mePos.size();
			int numEt = etPos.size();

			//if (mol.isLinear())
			//	name.append("N");
			//else
			//{
			//	bool alreadyNamed = false;
			//	if (numMe == 2 && numEt == 0)
			//	{
			//		if ((mePos[0] == mainChain[1] && mePos[1] == mainChain[1])
			//			|| (mePos[0] == mainChain[mainChain.size() - 2] && mePos[1] == mainChain[mainChain.size() - 2]))
			//		{
			//			name.append("NE");
			//			alreadyNamed = true;
			//		}
			//		else if ((mePos[0] == mainChain[mainChain.size() - 2] && mePos[1] == mainChain[1])
			//			|| (mePos[0] == mainChain[1] && mePos[1] == mainChain[mainChain.size() - 2]))
			//		{
			//			name.append("X");
			//			alreadyNamed = true;
			//		}
			//	}
			//	if (alreadyNamed == false)
			//	{
			//		name.append("I");
			//		for (int i = 0; i < numMe; i++)
			//			name.append(std::to_string(mePos[i]));
			//		if (numEt > 0)
			//		{
			//			for (int i = 0; i < numEt; i++)
			//				name.append(std::to_string(etPos[i]));
			//			name.append("ET");
			//		}
			//	}
			//
			//}
			name.append(alkaneStructurePrefix(&mol));

			name.append("C");
			name.append(std::to_string(mol.numberOfC()));

			switch (type)
			{
			case R_: // molecule is R
			{
				name.append("-");
				name.append(std::to_string(map[mol.posCrad()-1]));
				break;
			}

			case ROO_:	// molecule is ROO
			{
				name.append("-");
				name.append(std::to_string(map[mol.posCOOrad()-1]));
				name.append("O2");
				break;
			}

			case QOOH_:	// molecule is QOOH
			{
				name.append("OOH");
				name.append(std::to_string(map[mol.posCOOH()-1]));
				name.append("-");
				name.append(std::to_string(map[mol.posCrad()-1]));
				break;
			}

			case OOQOOH_:	// molecule is OOQOOH
			{
				name.append("OOH");
				name.append(std::to_string(map[mol.posCOOH() - 1]));
				name.append("-");
				name.append(std::to_string(map[mol.posCOOrad()-1]));
				name.append("O2");
				break;
			}

			case OLE_:	// molecule is olefin
			{
				// this must be double checked and modified if the numbering of branching is changed
				// for now it should prevent two different molecules to have the same name
				std::vector<int> posOle = mol.posOle(1);
				int pos1 = map[posOle[0]-1];
				int pos2 = map[posOle[1]-1];
				int pos = 0;
				if (pos1 <= mainChain.size() && pos2 <= mainChain.size())
					pos = std::min(pos1, pos2);
				else if (pos1 <= mainChain.size() && pos2 > mainChain.size())
					pos = pos2;
				else if (pos1 > mainChain.size() && pos2 <= mainChain.size())
					pos = pos1;
				else if (pos1 > mainChain.size() && pos2 > mainChain.size())
					pos = std::max(pos1, pos2);
				name.append("D");
				name.append(std::to_string(pos));
				break;
			}

			case CO_:	// molecule is CO
			{

				//name.append("-");
				//name.append(std::to_string(map[mol.posKeto() - 1]));
				//name.append("O");
				//std::cout << mol << std::endl;
				
				name = "";		// erase the base alkane name
				Molecola m1, m2;
				int posKeto = mol.posKeto();
				
				if (posKeto == mainChain[0])
				{
					mol.spezza(mainChain[0], mainChain[1], &m1, &m2);
					m2 = m2.parentFuel();
					//std::cout << m1 << std::endl;
					//std::cout << m2 << std::endl;
					if (m2.isLinear())
					{
						name.append("C");
						if (m2.numberOfC() > 1)
							name.append(std::to_string(m2.numberOfC()));
						if (m2.numberOfH() - 1 > 0)
						{
							name.append("H");
							if (m2.numberOfH()-1 > 1)
								name.append(std::to_string(m2.numberOfH() - 1));
						}
					}
					else
					{
						std::string na = molToName(m2);
						removeH(&na);
						name.append(na);
					}
					name.append("C");
					if (mol.numH(posKeto) > 0)
					{
						name.append("H");
						name.append(std::to_string(mol.numH(posKeto)));
					}
					name.append("O");
				}
				else if (posKeto == mainChain[mainChain.size() - 1])
				{
					mol.spezza(mainChain[mainChain.size()-2], mainChain[mainChain.size()-1], &m1, &m2);
					m1 = m1.parentFuel();
					//std::cout << m1 << std::endl;
					//std::cout << m2 << std::endl;
					if (m1.isLinear())
					{
						name.append("C");
						if (m1.numberOfC() > 1)
							name.append(std::to_string(m1.numberOfC()));
						if (m1.numberOfH() - 1 > 0)
						{
							name.append("H");
							if (m1.numberOfH() - 1 > 1)
								name.append(std::to_string(m1.numberOfH() - 1));
						}
					}
					else
					{
						std::string na = molToName(m1);
						removeH(&na);
						name.append(na);
					}
					name.append("C");
					if (mol.numH(posKeto) > 0)
					{
						name.append("H");
						if(mol.numH(posKeto) > 1)
							name.append(std::to_string(mol.numH(posKeto)));
					}
					name.append("O");
				}
				else
				{
					// find index of main chain to which the =O is attached
					int ind = 0;
					for (int i = 0; i < mainChain.size(); i++)
						if (mainChain[i] == posKeto)
							ind = i;
					Molecola m3;
					//std::cout << mol << std::endl;
					//std::cout << ind << std::endl;
					mol.spezza(mainChain[ind], mainChain[ind - 1], &m3, &m1);
					mol.spezza(mainChain[ind], mainChain[ind +1], &m3, &m2);
					m1 = m1.parentFuel();
					m2 = m2.parentFuel();
					//std::cout << m1 << std::endl;
					//std::cout << m2 << std::endl;
					//std::cout << m3 << std::endl;
					if (m1.isLinear())
					{
						name.append("C");
						if (m1.numberOfC() > 1)
							name.append(std::to_string(m1.numberOfC()));
						if (m1.numberOfH() - 1 > 0)
						{
							name.append("H");
							if (m1.numberOfH() - 1 > 1)
								name.append(std::to_string(m1.numberOfH() - 1));
						}
					}
					else
					{
						std::string na = molToName(m1);
						removeH(&na);
						name.append(na);
					}
					name.append("C");
					if (mol.numH(posKeto) > 0)
					{
						name.append("H");
						if (mol.numH(posKeto) > 1)
							name.append(std::to_string(mol.numH(posKeto)));
					}
					name.append("O");
					if (m2.isLinear())
					{
						name.append("C");
						if (m2.numberOfC() > 1)
							name.append(std::to_string(m2.numberOfC()));
						if (m2.numberOfH() - 1 > 0)
						{
							name.append("H");
							if (m2.numberOfH() - 1 > 1)
								name.append(std::to_string(m2.numberOfH() - 1));
						}
					}
					else
					{
						std::string na = molToName(m2);
						removeH(&na);
						name.append(na);
					}
				}

				break;
			}

			case cEth_:	// molecule is an ether
			{
				name.append("O");
				std::vector<int> pos = mol.posEthero();
				int pos1 = std::min(map[pos[0] - 1], map[pos[1] - 1]);
				int pos2 = std::max(map[pos[0] - 1], map[pos[1] - 1]);
				name.append(std::to_string(pos1));
				name.append("-");
				name.append(std::to_string(pos2));
				//std::cout << std::setw(20) << mol << std::setw(5) << pos[0] << std::setw(5) << pos[1] << "  ->  " << std::setw(5) << map[pos[0]-1] << std::setw(5) << map[pos[1]-1] << "    " << name << std::endl;
				break;
			}

			case RO_:	// molecule is an RO
			{
				name.append("-");
				name.append(std::to_string(map[mol.posCrad()-1]));
				name.append("-");
				name.append(std::to_string(map[mol.posKeto()-1]));
				name.append("O");
				break;
			}

			case KHP_:	// molecule is KHP
			{
				name.append("KET");
				name.append(std::to_string(map[mol.posKeto()-1]));
				name.append("-");
				name.append(std::to_string(map[mol.posCOOH()-1]));
				break;
			}

			case ROOH_:	// molecule is ROOH
			{
				name.append("-");
				name.append(std::to_string(mol.posCOOH()));
				name.append("OOH");
				break;
			}

			case POOH2_:	// molecule is P(OOH)2
			{
				std::vector<int> posOOH = mol.posOOHinPOOH2();
				int pos1 = std::min(map[posOOH[0] - 1], map[posOOH[1] - 1]);
				int pos2 = std::max(map[posOOH[0] - 1], map[posOOH[1] - 1]);
				name.append("Q");
				name.append(std::to_string(pos1));
				name.append(std::to_string(pos2));
				name.append("-");
				name.append(std::to_string(map[mol.posCrad()-1]));
				break;
			}

			case ROOH2_:	// molecule is R(OOH)2
			{
				std::vector<int> posOOH = mol.posOOHinROOH2();
				int pos1 = std::min(map[posOOH[0] - 1], map[posOOH[1] - 1]);
				int pos2 = std::max(map[posOOH[0] - 1], map[posOOH[1] - 1]);
				name.append("Q");
				name.append(std::to_string(pos1));
				name.append(std::to_string(pos2));
				break;
			}

			case oleOOH_:	// molecule is ole-OOH
			{
				// this must be double checked and modified if the numbering of branching is changed
				// for now it should prevent two different molecules to have the same name
				std::vector<int> posOle = mol.posOle(1);
				int pos1 = map[posOle[0] - 1];
				int pos2 = map[posOle[1] - 1];
				int pos = 0;
				if (pos1 <= mainChain.size() && pos2 <= mainChain.size())
					pos = std::min(pos1, pos2);
				else if (pos1 <= mainChain.size() && pos2 > mainChain.size())
					pos = pos2;
				else if (pos1 > mainChain.size() && pos2 <= mainChain.size())
					pos = pos1;
				else if (pos1 > mainChain.size() && pos2 > mainChain.size())
					pos = std::max(pos1, pos2);
				name.append("D");
				name.append(std::to_string(pos));
				name.append("-");
				name.append(std::to_string(map[mol.posCOOH()-1]));
				name.append("OOH");
				break;
			}

			case cEthOOH_:	// molecule is an ether-OOH
			{
				name.append("O");
				std::vector<int> pos = mol.posEthero();
				int pos1 = std::min(map[pos[0] - 1], map[pos[1] - 1]);
				int pos2 = std::max(map[pos[0] - 1], map[pos[1] - 1]);
				name.append(std::to_string(pos1));
				name.append(std::to_string(pos2));
				name.append("-");
				name.append(std::to_string(map[mol.posCOOH()-1]));
				name.append("OOH");
				break;
			}

			case oleR_: // molecule is oleR
			{
				// this must be double checked and modified if the numbering of branching is changed
				// for now it should prevent two different molecules to have the same name
				std::vector<int> posOle = mol.posOle(1);
				int pos1 = map[posOle[0] - 1];
				int pos2 = map[posOle[1] - 1];
				int pos = 0;
				if (pos1 <= mainChain.size() && pos2 <= mainChain.size())
					pos = std::min(pos1, pos2);
				else if (pos1 <= mainChain.size() && pos2 > mainChain.size())
					pos = pos2;
				else if (pos1 > mainChain.size() && pos2 <= mainChain.size())
					pos = pos1;
				else if (pos1 > mainChain.size() && pos2 > mainChain.size())
					pos = std::max(pos1, pos2);
				name.append("D");
				name.append(std::to_string(pos));
				name.append("-");
				name.append(std::to_string(map[mol.posCrad()-1]));
				break;
			}

			case oleCO_:	// molecule is oleCO
			{
				//std::stringstream nameMol;
				//nameMol << mol;
				//int asdf = 0;
				//if (nameMol.str() == "C=C(-C)-C(-C)-CO")
				//	asdf = 1;
				//if (nameMol.str() == "C-C(-C)-C(=C)-CO")
				//	asdf = 1;

				// this must be double checked and modified if the numbering of branching is changed
				// for now it should prevent two different molecules to have the same name
				bool isI3C6 = false;
				if (name == "I3C6")
					isI3C6 = true;

				std::vector<int> posOle = mol.posOle(1);
				int pos1 = map[posOle[0] - 1];
				int pos2 = map[posOle[1] - 1];
				int pos = 0;
				if (pos1 <= mainChain.size() && pos2 <= mainChain.size())
					pos = std::min(pos1, pos2);
				else if (pos1 <= mainChain.size() && pos2 > mainChain.size())
					pos = pos2;
				else if (pos1 > mainChain.size() && pos2 <= mainChain.size())
					pos = pos1;
				else if (pos1 > mainChain.size() && pos2 > mainChain.size())
					pos = std::max(pos1, pos2);
				name.append("D");
				name.append(std::to_string(pos));
				name.append("-");
				name.append(std::to_string(map[mol.posKeto() - 1]));
				name.append("O");
				
				if (isI3C6)
				{
					name.append("-CO");
				}

				//if (asdf == 1)
				//{
				//	std::cout << mol << std::endl;
				//	std::cout << posOle[0] << "   " << posOle[1] << std::endl;
				//	std::cout << pos1 << "   " << pos2 << std::endl;
				//	for (int as = 0; as < mainChain.size(); as++)
				//		std::cout << as+1 << "    ";
				//	std::cout << std::endl;
				//	for (int as = 0; as < mainChain.size(); as++)
				//		std::cout << mainChain[as] << "    ";
				//	std::cout << std::endl;
				//	std::cout << name << std::endl;
				//	std::cout << std::endl;
				//}
				
				break;
			}

			case cEthR_:	// molecule is an etherR
			{
				name.append("O");
				std::vector<int> pos = mol.posEthero();
				int pos1 = std::min(map[pos[0] - 1], map[pos[1] - 1]);
				int pos2 = std::max(map[pos[0] - 1], map[pos[1] - 1]);
				name.append(std::to_string(pos1));
				name.append(std::to_string(pos2));
				name.append("-");
				name.append(std::to_string(map[mol.posCrad() - 1]));
				break;
			}

			case cEthCO_:	// molecule is etherCO
			{
				name.append("O");
				std::vector<int> pos = mol.posEthero();
				int pos1 = std::min(map[pos[0] - 1], map[pos[1] - 1]);
				int pos2 = std::max(map[pos[0] - 1], map[pos[1] - 1]);
				name.append(std::to_string(pos1));
				name.append(std::to_string(pos2));
				name.append("-");
				name.append(std::to_string(map[mol.posKeto() - 1]));
				name.append("O");
				break;
			}

			case lEthRO_:	// molecule is linEtherRO
			{
				name.append("C");
				name.append(std::to_string(mol.numberOfC()));
				name.append("linEtherRO-");
				name.append(std::to_string(Add(mol, &linEtheROIsomersList) + 1));
				break;
			}

			case alkRO_:	// molecule is alkenyl RO
			{
				// this must be double checked and modified if the numbering of branching is changed
				// for now it should prevent two different molecules to have the same name
				std::vector<int> posOle = mol.posOle(1);
				int pos1 = map[posOle[0] - 1];
				int pos2 = map[posOle[1] - 1];
				int pos = 0;
				if (pos1 <= mainChain.size() && pos2 <= mainChain.size())
					pos = std::min(pos1, pos2);
				else if (pos1 <= mainChain.size() && pos2 > mainChain.size())
					pos = pos2;
				else if (pos1 > mainChain.size() && pos2 <= mainChain.size())
					pos = pos1;
				else if (pos1 > mainChain.size() && pos2 > mainChain.size())
					pos = std::max(pos1, pos2);
				name.append("-");
				name.append(std::to_string(pos));
				name.append("D");
				name.append(std::to_string(map[mol.posCOrad()-1]));
				name.append("O");
				break;
			}
			}

		}

	}


	if (name == "NEOC5H11CONEOC5H11")
		name = "NEOC5CONEOC5";
	if (name == "I2244C9CONEOC5H11")
		name = "I2244C9CONEOC5";
	//addSpeciesToList(type, name);
	return name;
}

void ChemkinOut::addSpeciesToKinFile(Molecola mol)
{
	species type = mol.kindOfSPecies();
	std::string name = molToName(mol);
	addSpeciesToList(type, name);
}

void ChemkinOut::addSpeciesToList(species type, std::string name)
{
	switch (type)
	{
	case special_:

		break;
	case fuel_:
		//addUniqueName(&, name);
		break;
	case R_:
		addUniqueName(&RList, name);
		break;
	case ROO_:
		addUniqueName(&ROOList, name);
		break;
	case QOOH_:
		addUniqueName(&QOOHList, name);
		break;
	case OOQOOH_:
		addUniqueName(&OOQOOHList, name);
		break;
	case OLE_:
		addUniqueName(&OLEList, name);
		break;
	case CO_:
		addUniqueName(&COList, name);
		break;
	case cEth_:
		addUniqueName(&EtherList, name);
		break;
	case RO_:
		addUniqueName(&ROList, name);
		break;
	case KHP_:
		addUniqueName(&KHPList, name);
		break;
	case ROOH_:
		addUniqueName(&ROOHList, name);
		break;
	case POOH2_:
		addUniqueName(&POOH2List, name);
		break;
	case oleOOH_:
		addUniqueName(&oleOOHList, name);
		break;
	case cEthOOH_:
		addUniqueName(&etherOOHList, name);
		break;
	case oleR_:
		addUniqueName(&oleRList, name);
		break;
	case oleCO_:
		addUniqueName(&oleCOList, name);
		break;
	case cEthR_:
		addUniqueName(&etherRList, name);
		break;
	case cEthCO_:
		addUniqueName(&etherCOList, name);
		break;
	case lEthRO_:
		addUniqueName(&linEtheROList, name);
		break;
	case alkRO_:
		addUniqueName(&alkROList, name);
		break;
	case ROOH2_:
		break;
	default:
		std::cerr << "WARNING! Tried to add " << name << " to the list of species. It has not beign recognized!" << std::endl;
		break;
	}
}

std::string ChemkinOut::molToNameLump(Molecola mol)
{
	std::string name;
	int type = -1;
	// search the type of molecule:
	//		1  = fuel
	//		2  = R
	//		3  = ROO
	//		4  = QOOH
	//		5  = OOQOOH
	//		6  = OLE
	//		7  = CO
	//		8  = Ether
	//		9  = RO
	//		10 = KHP
	//		11 = ROOH
	//		12 = P(OOH)2
	//		13 = oleOOH
	//		14 = cycOOH
	//		15 = oleR
	//		16 = oleCO
	//		17 = etherR
	//		18 = etherCO
	//		19 = linEtheRO
	//		20 = alkenylRO

	int numberR = mol.numCrad();
	int numberCOO = mol.numCOORad();
	int numberCOOH = mol.numCOOH();
	int numberEther = mol.numberEthero();
	int numberOLE = mol.numberOle();
	int numberCO = mol.numKeto();
	int specMol = mol.isSpecialMolecule();
	int numLinEther = mol.numLinEther();
	int numCORad = mol.numCORad();

	if (specMol == 0 && mol.size() <= 4)
	{
		std::string n = checkIfNamed(&mol);
		if (n != "not found")
		{
			return n;
		}
		else
		{
			std::cerr << "WARNING: was not able to find " << mol << " in the base model list of names" << "\n";
		}
	}

	std::string prefix = "";
	if (specMol == 0 && mol.size() > 2)
	{
		Molecola baseStruct = mol.parentFuel();
		prefix = alkaneStructurePrefix(&baseStruct);
	}

	if (specMol != 0)																	 // molecule is special molecule
	{
		type = 0;
		switch (specMol)
		{
		case 1:
			name.append("N2");
			break;
		case 2:
			name.append("O2");
			break;
		case 3:
			name.append("O");
			break;
		case 4:
			name.append("OH");
			break;
		case 5:
			name.append("HO2");
			break;
		case 6:
			name.append("H2O");
			break;
		case 7:
			name.append("H2O2");
			break;
		case 8:
			name.append("H");
			break;
		case 9:
			name.append("H2");
			break;
		case 10:
			name.append("CO");
			break;

		}
	}
	else if (numberR + numberCOO + numberCOOH + numberEther + numberOLE + numberCO + numLinEther + numCORad == 0) // molecule is alkane
	{
		type = 1;
		name.append(prefix);
		name.append("C");
		name.append(std::to_string(mol.numberOfC()));
		name.append("H");
		name.append(std::to_string(mol.numberOfH()));
		if (mol.numberOfC() == 1)			// if it is methane change the name in CH4 instead of C1H4
			name = "CH4";
	}
	else if (numberR == 1 &&
		numberCOO + numberCOOH + numberEther + numberOLE + numberCO + numLinEther + numCORad == 0) // molecule is R
	{
		type = 2;
		if (prefix != "N")
			name.append(prefix);
		name.append("R");
		name.append(std::to_string(mol.numberOfC()));
		if (mol.numberOfC() == 1)			// if it is methane radical change the name in CH3 instead of R1
			name = "CH3";
		if (mol.numberOfC() == 2)			// if it is ethilene radical change the name in C2H5 instead of R2
			name = "C2H5";

	}
	else if (numberCOO == 1 &&
		numberR + numberCOOH + numberEther + numberOLE + numberCO + numLinEther + numCORad == 0)	// molecule is ROO
	{
		type = 3;
		if (prefix != "N")
			name.append(prefix);
		name.append("R");
		name.append(std::to_string(mol.numberOfC()));
		name.append("OO");
	}
	else if (numberCOOH == 1 && numberR == 1 &&
		numberCOO + numberEther + numberOLE + numberCO + numLinEther + numCORad == 0)	// molecule is QOOH
	{
		type = 4;
		if (prefix != "N")
			name.append(prefix);
		name.append("Q");
		name.append(std::to_string(mol.numberOfC()));
		name.append("OOH");
	}
	else if (numberCOO == 1 && numberCOOH == 1 &&
		numberR + numberEther + numberOLE + numberCO + numLinEther + numCORad == 0)	// molecule is OOQOOH
	{
		type = 5;
		name.append("OO");
		if (prefix != "N")
			name.append(prefix);
		name.append("Q");
		name.append(std::to_string(mol.numberOfC()));
		name.append("OOH");
	}
	else if (numberOLE == 1 &&
		numberR + numberCOO + numberCOOH + numberEther + numberCO + numLinEther + numCORad == 0)	// molecule is olefin
	{
		type = 6;
		if (prefix != "N")
			name.append(prefix);
		name.append("OLE");
		name.append(std::to_string(mol.numberOfC()));
	}
	else if (numberCO == 1 &&
		numberR + numberCOO + numberCOOH + numberEther + numberOLE + numLinEther + numCORad == 0)	// molecule is CO
	{
		type = 7;
		if (prefix != "N")
			name.append(prefix);
		name.append("C");
		name.append(std::to_string(mol.numberOfC()));
		name.append("O");
	}
	else if (numberEther == 1 &&
		numberR + numberCOO + numberCOOH + numberOLE + numberCO + numLinEther + numCORad == 0)	// molecule is an ether
	{
		type = 8;
		if (prefix != "N")
			name.append(prefix);
		name.append("CETHE");
		name.append(std::to_string(mol.numberOfC()));
	}
	else if (numberCO == 1 && numberR == 1 &&
		numberCOO + numberCOOH + numberEther + numberOLE + numLinEther + numCORad == 0)	// molecule is an RO
	{
		type = 9;
		if (prefix != "N")
			name.append(prefix);
		name.append("R");;
		name.append(std::to_string(mol.numberOfC()));
		name.append("O");
	}
	else if (numberCO == 1 && numberCOOH == 1 &&
		numberR + numberCOO + numberEther + numberOLE + numLinEther + numCORad == 0)	// molecule is KHP
	{
		type = 10;
		if (prefix != "N")
			name.append(prefix);
		name.append("KHP");
		name.append(std::to_string(mol.numberOfC()));
	}
	else if (numberCOOH == 1 &&
		numberR + numberCOO + numberEther + numberOLE + numberCO + numLinEther + numCORad == 0)	// molecule is ROOH
	{
		type = 11;
		if (prefix != "N")
			name.append(prefix);
		name.append("R");
		name.append(std::to_string(mol.numberOfC()));
		name.append("OOH");
	}
	else if (numberCOOH == 2 && numberR == 1 &&
		numberCOO + numberEther + numberOLE + numberCO + numLinEther + numCORad == 0)	// molecule is P(OOH)2
	{
		type = 12;
		if (prefix != "N")
			name.append(prefix);
		name.append("P");
		name.append(std::to_string(mol.numberOfC()));
		name.append("(OOH)2");
	}
	else if (numberCOOH == 1 && numberOLE == 1 &&
		numberR + numberCOO + numberEther + numberCO + numLinEther + numCORad == 0)	// molecule is ole-OOH
	{
		type = 13;
		if (prefix != "N")
			name.append(prefix);
		name.append("OLE");
		name.append(std::to_string(mol.numberOfC()));
		name.append("OOH");
	}
	else if (numberEther == 1 && numberCOOH == 1 &&
		numberR + numberCOO + numberOLE + numberCO + numLinEther + numCORad == 0)	// molecule is an ether-OOH
	{
		type = 14;
		if (prefix != "N")
			name.append(prefix);
		name.append("CETHE");
		name.append(std::to_string(mol.numberOfC()));
		name.append("OOH");
	}
	else if (numberR == 1 && numberOLE == 1 &&
		numberCOO + numberCOOH + numberEther + numberCO + numLinEther + numCORad == 0) // molecule is oleR
	{
		type = 15;
		if (prefix != "N")
			name.append(prefix);
		name.append("OLE");
		name.append(std::to_string(mol.numberOfC()));
		name.append("R");
	}
	else if (numberOLE == 1 && numberCO == 1 &&
		numberR + numberCOO + numberCOOH + numberEther + numLinEther + numCORad == 0)	// molecule is oleCO
	{
		type = 16;
		if (prefix != "N")
			name.append(prefix);
		name.append("OLE");
		name.append(std::to_string(mol.numberOfC()));
		name.append("O");

	}
	else if (numberEther == 1 && numberR == 1 &&
		numberCOO + numberCOOH + numberOLE + numberCO + numLinEther + numCORad == 0)	// molecule is an etherR
	{
		type = 17;
		if (prefix != "N")
			name.append(prefix);
		name.append("CETHE");
		name.append(std::to_string(mol.numberOfC()));
		name.append("R");
	}

	else if (numberCO == 1 && numberEther == 1 &&
		numberR + numberCOO + numberCOOH + numberOLE + numLinEther + numCORad == 0)	// molecule is etherCO
	{
		type = 18;
		if (prefix != "N")
			name.append(prefix);
		name.append("CETHE");
		name.append(std::to_string(mol.numberOfC()));
		name.append("O");
	}

	else if (numberCO == 1 && numLinEther == 1 && numberR == 1 &&
		numberCOO + numberCOOH + numberOLE + numberEther + numCORad == 0)	// molecule is linEtherRO
	{
		type = 19;
		if (prefix != "N")
			name.append(prefix);
		name.append("LETHE");
		name.append(std::to_string(mol.numberOfC()));
		name.append("RO");
	}
	else if (numberOLE == 1 && numCORad == 1 &&
		numLinEther + numberCO + numberR + numberCOO + numberCOOH + numberEther == 0)	// molecule is alkenyl RO
	{
		type = 20;
		name.append("ALK");
		if (prefix != "N")
			name.append(prefix);
		name.append("R");
		name.append(std::to_string(mol.numberOfC()));
		name.append("O");
	}
	if (type == -1)
	{
		std::cerr << "Error in naming the lumped molecule " << mol << ", it doesn't match any expected molecule type." << std::endl;
		exit(0);
	}
	return name;
}

void ChemkinOut::wrireaDetailed(Reaction rea)
{
	std::string firstMember;
	std::string secondMember;
	std::vector<Molecola*> reactants = rea.reactantList();
	std::vector<Molecola*> products = rea.productList();
	for (int i = 0; i < reactants.size(); i++)
	{
		firstMember.append(molToName(*(reactants[i])));
		if (i != reactants.size() - 1) firstMember.append(" + ");
	}
	for (int i = 0; i < products.size(); i++)
	{
		secondMember.append(molToName(*(products[i])));
		if (i != products.size() - 1) secondMember.append(" + ");
	}

	char cost[100];
	sprintf(cost, "   %9.3e %8.3f %8.0f ", rea.A(), rea.n(), rea.E());

	std::vector<std::string> tempStr = { firstMember, secondMember };

	std::vector<int> dupIndexes = searchDuplicates(&detReacSpecies, tempStr);	// search for duplicates
	detReacSpecies.push_back(tempStr);
	detReacConst.push_back(cost);
	detReacComments.push_back(rea.printComment());
	isDetailedReactionDuplicate.push_back(false);

	if (dupIndexes.size() > 0) isDetailedReactionDuplicate[isDetailedReactionDuplicate.size() - 1] = true;
	for (int i = 0; i < dupIndexes.size(); i++)			// flag the duplicates
		isDetailedReactionDuplicate[dupIndexes[i]] = true;
}



void ChemkinOut::writeDetailedReactions(std::vector<Reaction>* vec)
{

	addUniqueName(&genericSpecies, "O2");
	addUniqueName(&genericSpecies, "H2");
	addUniqueName(&genericSpecies, "H2O2");
	addUniqueName(&genericSpecies, "H2O");
	addUniqueName(&genericSpecies, "CH4");
	addUniqueName(&genericSpecies, "C2H6");
	addUniqueName(&genericSpecies, "N2");
	addUniqueName(&genericSpecies, "AR");
	addUniqueName(&genericSpecies, "HE");
	addUniqueName(&genericSpecies, "CO");
	addUniqueName(&genericRadicals, "HO2");
	addUniqueName(&genericRadicals, "OH");
	addUniqueName(&genericRadicals, "H");
	addUniqueName(&genericRadicals, "O");
	addUniqueName(&RList, "CH3");
	addUniqueName(&RList, "C2H5");


	std::string label = (*vec)[0].reactionLabel();
	detReacName.push_back(label);
	detReacNamePosition.push_back(detReacConst.size());
	wrireaDetailed((*vec)[0]);

	for (int i = 1; i < vec->size(); i++)
	{
		if ((*vec)[i].reactionLabel() != label)
		{
			label = (*vec)[i].reactionLabel();
			detReacName.push_back(label);
			detReacNamePosition.push_back(detReacConst.size());
		}
		wrireaDetailed((*vec)[i]);
	}
}

void ChemkinOut::wrireaDetailed(int nreaz, double A, double n, double E, Molecola m1, ...)
/////////////////////////////////////////////////////////////////////////////
//  SCRIVE le reazioni sul file stream
//
//  1  A --> B
//  2  A --> B + C
//  3  A + O2 --> B
//  4  A --> B +O2
//  5  A +O2 --> B + HO2
//  6  A --> B + HO2
//  7  A --> B + OH
//  8  A --> B + C + OH
//  9  A  + RH --> B + R*
// 10  A --> B + C + 2 OH
// 11  A --> B + C + HO2
// 12  A --> B + C + D + OH  
// 13  R* + A --> B
// 14  O2 + A --> B + HO2
// 15  OH + A --> B + H2O
// 16  H + A --> B + H2
// 17  O + A --> B + OH
// 18  HO2 + A --> B + H2O2
// 19  CH3 + A --> B + CH3
// 20  C2H5 + A --> B + C2H6
//
/////////////////////////////////////////////////////////////////////////////
{
	char cost[100];   // per contenere le costanti
	m1.output = FORMULA;
	int size = m1.size() * 9;
	if (size > 97) size = 97;

	va_list elenco;
	va_start(elenco, m1);

	std::string firstMember;
	std::string secondMember;

	switch (nreaz)
	{
	case 3: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" + O2");
		secondMember.append(molToName(m));

		addUniqueName(&genericSpecies, "O2");
		break; };
	case 1: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		secondMember.append(molToName(m));
		break; };
	case 2: { Molecola m[2];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		secondMember.append(molToName(m[0]));
		secondMember.append(" + ");
		secondMember.append(molToName(m[1]));
		break; };
	case 4: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		secondMember.append(molToName(m));
		secondMember.append(" + O2");

		addUniqueName(&genericSpecies, "O2");
		break; };
	case 5: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" + O2");
		secondMember.append(molToName(m));
		secondMember.append(" + HO2");

		addUniqueName(&genericSpecies, "O2");
		addUniqueName(&genericRadicals, "HO2");
		break; };
	case 6: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		secondMember.append(molToName(m));
		secondMember.append(" + HO2");

		addUniqueName(&genericRadicals, "HO2");
		break; };
	case 7: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		secondMember.append(molToName(m));
		secondMember.append(" + OH");

		addUniqueName(&genericRadicals, "OH");
		break; };
	case 8: { Molecola m[2];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		secondMember.append(molToName(m[0]));
		secondMember.append(" + ");
		secondMember.append(molToName(m[1]));
		secondMember.append(" + OH");

		addUniqueName(&genericRadicals, "OH");
		break; };
	case 9: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" + RH");
		secondMember.append(molToName(m));
		secondMember.append(" + R*");
		break; };
	case 10: { Molecola m[2];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		secondMember.append(molToName(m[0]));
		secondMember.append(" + ");
		secondMember.append(molToName(m[1]));
		secondMember.append(" + 2 OH");

		addUniqueName(&genericRadicals, "OH");
		break; };
	case 11: { Molecola m[2];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		secondMember.append(molToName(m[0]));
		secondMember.append(" + ");
		secondMember.append(molToName(m[1]));
		secondMember.append(" + HO2");

		addUniqueName(&genericRadicals, "HO2");
		break; };
	case 12: { Molecola m[3];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		m[2] = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		secondMember.append(molToName(m[0]));
		secondMember.append(" + ");
		secondMember.append(molToName(m[1]));
		secondMember.append(" + ");
		secondMember.append(molToName(m[2]));
		secondMember.append(" + OH");

		addUniqueName(&genericRadicals, "OH");
		break; };
	case 13: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" + R*");
		secondMember.append(molToName(m));
		break; };
	case 14: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" +   O2");
		secondMember.append(molToName(m));
		secondMember.append(" +  HO2");

		addUniqueName(&genericSpecies, "O2");
		addUniqueName(&genericRadicals, "HO2");
		break; };
	case 15: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" +   OH");
		secondMember.append(molToName(m));
		secondMember.append(" +  H2O");

		addUniqueName(&genericRadicals, "OH");
		addUniqueName(&genericRadicals, "HO2");
		break; };
	case 16: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" +    H");
		secondMember.append(molToName(m));
		secondMember.append(" +   H2");

		addUniqueName(&genericRadicals, "H");
		addUniqueName(&genericSpecies, "H2");
		break; };
	case 17: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" +    O");
		secondMember.append(molToName(m));
		secondMember.append(" +   OH");

		addUniqueName(&genericRadicals, "O");
		addUniqueName(&genericRadicals, "OH");
		break; };
	case 18: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" +  HO2");
		secondMember.append(molToName(m));
		secondMember.append(" + H2O2");

		addUniqueName(&genericRadicals, "HO2");
		addUniqueName(&genericSpecies, "H2O2");
		break; };
	case 19: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" +  CH3");
		secondMember.append(molToName(m));
		secondMember.append(" +  CH4");

		addUniqueName(&genericRadicals, "CH3");
		addUniqueName(&genericSpecies, "CH4");
		break; };
	case 20: { Molecola m = va_arg(elenco, Molecola);
		firstMember.append(molToName(m1));
		firstMember.append(" + C2H5");
		secondMember.append(molToName(m));
		secondMember.append(" + C2H6");

		addUniqueName(&genericRadicals, "C2H5");
		addUniqueName(&genericSpecies, "C2H6");
		break; };
	default:
		cerr << "Error in using wrirea! Not such case present in the options.";
		exit(200);
		break;
	};
	sprintf(cost, "   %9.3e %8.3f %8.0f ", A, n, E);

	std::vector<std::string> tempStr = { firstMember, secondMember };

	std::vector<int> dupIndexes = searchDuplicates(&detReacSpecies, tempStr);	// search for duplicates
	detReacSpecies.push_back(tempStr);
	detReacConst.push_back(cost);
	isDetailedReactionDuplicate.push_back(false);

	if (dupIndexes.size() > 0) isDetailedReactionDuplicate[isDetailedReactionDuplicate.size() - 1] = true;
	for (int i = 0; i < dupIndexes.size(); i++)			// flag the duplicates
		isDetailedReactionDuplicate[dupIndexes[i]] = true;

	//if (!addUniqueName(&reactionListDUP, firstMember + secondMember))
	//	detailedReactions << "DUP" << std::endl;
	//detailedReactions << std::left << setw(20) << firstMember << " => " << setw(30) << secondMember << "    " << cost << endl;
}

void ChemkinOut::writeDetReacName(std::string name)					// add the name of the detailed reaction in a vector and add the index of the first reaction of this family -1 to the index vector
{
	detReacName.push_back(name);
	detReacNamePosition.push_back(detReacConst.size());
}

void ChemkinOut::endlineDetailed()
{
	detailedReactions << std::endl;
}

void ChemkinOut::printAllListsDetailed()
{
	std::cout << "" << std::endl;
	std::cout << "Fuel" << std::endl;
	std::cout << fuel << std::endl;
	std::cout << "generic species list" << std::endl;
	printVector(&genericSpecies);
	std::cout << "generic radicals list" << std::endl;
	printVector(&genericRadicals);
	std::cout << "R list" << std::endl;
	printVector(&RList);
	std::cout << "ROO list" << std::endl;
	printVector(&ROOList);
	std::cout << "QOOH list" << std::endl;
	printVector(&QOOHList);
	std::cout << "OLE list" << std::endl;
	printVector(&OLEList);
	std::cout << "CO list" << std::endl;
	printVector(&COList);
	std::cout << "Ether list" << std::endl;
	printVector(&EtherList);
	std::cout << "OOQOOH list" << std::endl;
	printVector(&OOQOOHList);
	std::cout << "RO list" << std::endl;
	printVector(&ROList);
	std::cout << "KHP list" << std::endl;
	printVector(&KHPList);
	std::cout << "ROOH list" << std::endl;
	printVector(&ROOHList);
}


void ChemkinOut::completeDetailed()
{
	outFileDetailed << detailedHeading.rdbuf();

	// print elements
	outFileDetailed << "ELEMENTS" << std::endl << std::endl;
	outFileDetailed << "C \nH \nN \nO \nAR \nHE" << std::endl;
	outFileDetailed << std::endl << "END" << std::endl << std::endl << std::endl;

	// print species list
	outFileDetailed << "SPECIES" << std::endl << std::endl;
	addUniqueName(&genericSpecies, "N2");
	addUniqueName(&genericSpecies, "AR");
	addUniqueName(&genericSpecies, "HE");

	std::vector<std::string> listNames = { "Generic species", "Generic radicals", "R species", "ROO species", "QOOH species", "olefins", "aldehydes/ketones",
											"Cyclic ethers", "OOQOOH species", "RO species", "ketohydroperoxides", "P(OOH)2 species", "Olefins-OOH", "Cyc-ethers-OOH",
											"Radical olefins", "CO olefins", "Radical cyclic ethers", "Cyclyc ethers CO", "Linear ethers RO", "Alkenyl RO" };
	std::vector<std::vector<std::string>*> listVec = { &genericSpecies, &genericRadicals, &RList, &ROOList, &QOOHList, &OLEList, &COList,
															&EtherList, &OOQOOHList, &ROList, &KHPList, &POOH2List, &oleOOHList, &etherOOHList,
															&oleRList, &oleCOList, &etherRList, &etherCOList, &linEtheROList, &alkROList };

	outFileDetailed << "! Fuel" << std::endl;
	outFileDetailed << fuel << std::endl << std::endl;
	for (int i = 0; i < listVec.size(); i++)
	{
		outFileDetailed << "!" << listNames[i] << std::endl;
		for (int j = 0; j < listVec[i]->size(); j++)
		{
			outFileDetailed << std::left << std::setw(15) << (*listVec[i])[j] << " ";
			if ((j + 1) % 4 == 0 || (j + 1) == listVec[i]->size()) outFileDetailed << std::endl;
		}
		outFileDetailed << std::endl;
	}
	outFileDetailed << "END" << std::endl << std::endl << std::endl;

	// print reactions
	outFileDetailed << "REACTIONS" << std::endl;
	//outFileDetailed << detailedReactions.rdbuf();

	std::string lastName;
	for (int i = 0; i < detReacSpecies.size(); i++)
	{
		int namepos = findPosition(&detReacNamePosition, i);							// check if at this point there is a reaction name to print, if yes save the index of the name to print
		
		if (namepos != -1)
		{
			outFileDetailed << std::endl << "!" << detReacName[namepos] << std::endl;
			lastName = detReacName[namepos];
		}

#ifndef EQUILIBRIUM_MECH
		outFileDetailed << std::left << setw(20) << detReacSpecies[i][0] << " => " << setw(30) << detReacSpecies[i][1] << "    " << detReacConst[i] 
			<< "   ! " << detReacComments[i] <<	endl;	// print the reaction
		if (isDetailedReactionDuplicate[i])					// if the reaction has some duplicates add the DUP kaywork afterwards
			outFileDetailed << "DUP" << std::endl;
#else // !EQUILIBRIUM_MECH
		int stateReac = 0;   // 0: normal reaction, 1: equilibrium reaction, 2: commented reaction
		if (lastName == "O2 addition to R:" ||
			lastName == "O2 addition to QOOH:" ||
			lastName == "Isomerization ROO:" ||
			lastName == "Isomerization OOQOOH:")
		{
			stateReac = 1;
		}
		if (lastName == "O2 elimination from ROO:" ||
			lastName == "O2 elimination from OOQOOH:" ||
			lastName == "Isomerization P(OOH)2:" ||
			lastName == "Isomerization QOOH : ")
		{
			stateReac = 2;
		}
		
		if (stateReac == 2)
			outFileDetailed << "!";
		outFileDetailed << std::left << setw(20) << detReacSpecies[i][0];
		if(stateReac == 1)
			outFileDetailed << " <=> ";
		else
			outFileDetailed << " => ";
		
		outFileDetailed << setw(30) << detReacSpecies[i][1] << "    " << detReacConst[i]
			<< "   ! " << detReacComments[i] << endl;	// print the reaction
		if (isDetailedReactionDuplicate[i])					// if the reaction has some duplicates add the DUP kaywork afterwards
		{
			if (stateReac == 2)
				outFileDetailed << "!";
			outFileDetailed << "DUP" << std::endl;
		}
#endif
	}

	outFileDetailed << std::endl << "END" << std::endl;		// after printing all the reactions add the keyword END
	outFileDetailed.close();
}


std::vector<double> ChemkinOut::roundStoicCoeff(std::vector<double>* vect, int numOfDigit)		// return a vector with the element of vect rounded leaving numOfDifit of decimal digits
{
	std::vector<double> roundedVec;
	for (int i = 0; i < vect->size(); i++) roundedVec.push_back(round((*vect)[i] * pow(10.0, numOfDigit)) / pow(10.0, numOfDigit));
	return roundedVec;
}

std::string toString(double num, int numOfDigit)
{
	std::string res = std::to_string(num);
	for (int i = 0; i < 6 - numOfDigit; i++) res.pop_back();
	return res;
}

void ChemkinOut::wrireaLump(Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E)
{
	// write the lumped kinetic reaction in the form reag -> stoicCoeff[0] prod[1] + stoicCoeff[1] prod[2] + ...      A n E
	std::string reac;
	std::string species = molToName(reag);
	addUniqueName(&lumpedSpeciesList, species);
	reac += species + " => ";
	for (int i = 0; i < stoicCoeff.size(); i++)
	{
		species = "R" + std::to_string(prod[i + 1].numberOfC());
		if (species == "R1") species = "CH3";
		else if (species == "R2") species = "C2H5";
		else if (species == "R3") species = "NC3H7";
		addUniqueName(&lumpedSpeciesList, species);
		reac += std::to_string(stoicCoeff[i]) + " " + species;
		if (i != stoicCoeff.size() - 1) reac += " + ";
	}
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f", A, n, E);
	reac = reac + con + "\n";

	if (reac.size() < CHEMKIN_MAX_CHAR_PER_ROW + 2)  // if the reaction is short enough to be read by chemkin print it directly
	{
		lumpedReactions << reac;
	}
	else											// else eliminate all the spaces and round the stoic coeff until it is short enough
	{
		lumpedReactions << "!" << reac;		// write the complete reaction as a comment
		bool sizeReached = false;
		int numOfDigit = 6;
		while (!sizeReached)
		{
			reac.clear();
			std::vector<double> roundedStoich = roundStoicCoeff(&stoicCoeff, numOfDigit);
			species = molToName(reag);
			addUniqueName(&lumpedSpeciesList, species);
			reac = species + "=>";
			for (int i = 0; i < stoicCoeff.size(); i++)
			{
				species = "R" + std::to_string(prod[i + 1].numberOfC());
				if (species == "R1") species = "CH3";
				else if (species == "R2") species = "C2H5";
				else if (species == "R3") species = "NC3H7";
				reac += toString(roundedStoich[i], numOfDigit) + species;
				if (i != stoicCoeff.size() - 1) reac += "+";
			}
			char con[50];
			sprintf(con, " %.3e %.3f %.0f", A, n, E);
			reac = reac + con + "\n";
			if (reac.size() < CHEMKIN_MAX_CHAR_PER_ROW + 2)
			{
				sizeReached = true;
			}
			else numOfDigit--;
			if (numOfDigit == 1)
			{
				std::cerr << "WARNING: In order to shorten the initiation reaction under the maximum length admitted by chemkin, a number of significative digit less than 2 is required!" << std::endl;
				break;
			}
		}
		lumpedReactions << reac;
	}

	if (stoicCoeff.size() + 1 > max_species) max_species = stoicCoeff.size() + 1;
}

int ChemkinOut::addUniqueLumpedSpecies(Molecola* mol)
{
	std::string name = molToNameLump(*mol);
	for (int i = 0; i < lumpedNames.size(); i++)	// check if it is already present
		if (name == lumpedNames[i])					// if it is ...
			return i;								// ... return the index

	lumpedNames.push_back(name);					// otherwise add it
	lumpedMolecules.push_back(mol);
	return lumpedNames.size() - 1;					// and return its index
}

void ChemkinOut::wrireaLump(std::string label, std::vector<Molecola*> reactants, std::vector<Molecola*> products, double A, double n, double E)
{
	std::vector<double> stoic(products.size(), 1.0);
	wrireaLump(label, reactants, stoic, products, A, n, E);
}

void ChemkinOut::wrireaLump(std::string label, std::vector<Molecola*> reactants, std::vector<Molecola*> products, vector<double> Press, vector<double> A, vector<double> n, vector<double> E, int indexP)
{
	std::vector<double> stoic(products.size(), 1.0);
	wrireaLump(label, reactants, stoic, products, Press, A, n, E, indexP);
}

double round(double num, int numOfDecimals)
{
	return std::round(num * std::pow(10.0, numOfDecimals)) / pow(10.0, numOfDecimals);
}

void ChemkinOut::wrireaLump(std::string label, std::vector<Molecola*> reactants, std::vector<double> stoicCoeff, std::vector<Molecola*> products, double A, double n, double E)
{
	// write the reaction label
	lumpedReactions << "! " << label << std::endl;
	// approximate the stoichiometric coefficients
	int numOfDecimalDigits = 4;
	for (int i = 0; i < stoicCoeff.size(); i++)
	{
		stoicCoeff[i] = round(stoicCoeff[i], numOfDecimalDigits);
	}

	std::vector<std::string> namesLumped;		// the list of the names of product after grouping toghether all the products of the same family
	std::vector<double> stoicLumped;			// the stoichiometric coefficients of the species of namesLumped (after grouping together the species belonging to the same class)
	std::vector<Molecola*> refSpeciesLumped;	// a reference species of the groups of the vector namesLumped

	for (int i = 0; i < products.size(); i++)
	{
		std::string name = molToNameLump(*(products[i]));
		bool wasPresent = false;
		for (int j = 0; j < namesLumped.size(); j++)
		{
			if (name == namesLumped[j])
			{
				stoicLumped[j] += stoicCoeff[i];
				wasPresent = true;
				break;
			}
		}
		if (wasPresent == false)
		{
			namesLumped.push_back(name);
			stoicLumped.push_back(stoicCoeff[i]);
			refSpeciesLumped.push_back(products[i]);
		}
	}


	for (int i = 0; i < reactants.size(); i++)
	{
		std::string name = molToNameLump(*(reactants[i]));
		addUniqueLumpedSpecies(reactants[i]);
		lumpedReactions << name;
		if (i != reactants.size() - 1)
			lumpedReactions << "+";
	}
	lumpedReactions << "=>";
	int numSpecies = 0;
	bool fistProdPrinted = false;
	for (int i = 0; i < namesLumped.size(); i++)
	{
		if (stoicLumped[i] == 0.0)
			continue;
		if (fistProdPrinted == true)
			lumpedReactions << "+";
		else
			fistProdPrinted = true;
		std::string name = namesLumped[i];
		addUniqueLumpedSpecies(refSpeciesLumped[i]);
		if (stoicLumped[i] != 1.0)
			lumpedReactions << stoicLumped[i];
		lumpedReactions << name;
		numSpecies++;
	}
	char con[50];
	sprintf(con, "  %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	lumpedReactions << std::endl;


	if (numSpecies + 1 > max_species) max_species = numSpecies + 1;		// update maxSpecies
}

void ChemkinOut::wrireaLump(std::string label, std::vector<Molecola*> reactants, std::vector<double> stoicCoeff, std::vector<Molecola*> products, vector<double> Press, vector<double> A, vector<double> n, vector<double> E, int indexP)
{
	// write the reaction label
	lumpedReactions << "! " << label << std::endl;
	// approximate the stoichiometric coefficients
	int numOfDecimalDigits = 4;
	for (int i = 0; i < stoicCoeff.size(); i++)
	{
		stoicCoeff[i] = round(stoicCoeff[i], numOfDecimalDigits);
	}

	std::vector<std::string> namesLumped;		// the list of the names of product after grouping toghether all the products of the same family
	std::vector<double> stoicLumped;			// the stoichiometric coefficients of the species of namesLumped (after grouping together the species belonging to the same class)
	std::vector<Molecola*> refSpeciesLumped;	// a reference species of the groups of the vector namesLumped

	for (int i = 0; i < products.size(); i++)
	{
		std::string name = molToNameLump(*(products[i]));
		bool wasPresent = false;
		for (int j = 0; j < namesLumped.size(); j++)
		{
			if (name == namesLumped[j])
			{
				stoicLumped[j] += stoicCoeff[i];
				wasPresent = true;
				break;
			}
		}
		if (wasPresent == false)
		{
			namesLumped.push_back(name);
			stoicLumped.push_back(stoicCoeff[i]);
			refSpeciesLumped.push_back(products[i]);
		}
	}

	std::cout << label << std::endl;
	for (int i = 0; i < reactants.size(); i++)
	{
		std::cout << *(reactants[i]) << std::endl;
		std::string name = molToNameLump(*(reactants[i]));
		addUniqueLumpedSpecies(reactants[i]);
		lumpedReactions << name;
		if (i != reactants.size() - 1)
			lumpedReactions << "+";
	}
	lumpedReactions << "=>";
	int numSpecies = 0;
	bool fistProdPrinted = false;
	for (int i = 0; i < namesLumped.size(); i++)
	{
		if (stoicLumped[i] == 0.0)
			continue;
		if (fistProdPrinted == true)
			lumpedReactions << "+";
		else
			fistProdPrinted = true;
		std::string name = namesLumped[i];
		addUniqueLumpedSpecies(refSpeciesLumped[i]);
		if (stoicLumped[i] != 1.0)
			lumpedReactions << stoicLumped[i];
		lumpedReactions << name;
		numSpecies++;
	}
	char con[50];
	sprintf(con, "  %9.3e %8.3f %8.0f ", A[indexP], n[indexP], E[indexP]);
	
	lumpedReactions << con << std::endl;
	for (int i = 0; i < Press.size(); i++)
	{
		char plo[60];
		sprintf(plo, "PLOG / %7.3f       %9.3e %8.3f %8.0f  /",Press[i], A[i], n[i], E[i]);
		lumpedReactions << plo << std::endl;
	}
	lumpedReactions << std::endl;


	if (numSpecies + 1 > max_species) max_species = numSpecies + 1;		// update maxSpecies
}

void ChemkinOut::wrireaLump1(Molecola reag, HAbsRad rad, int numC, double A, double n, double E)
{
	// write the lumped kinetic reaction in the form reag + rad -> R_numC      A n E
	std::string species = molToName(reag);
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " + ";
	species = nameHAbsRad(rad);
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " => ";
	species = "R" + std::to_string(numC);
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " + ";
	species = nameHAbsRadPlusH(rad);
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species;
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;

}

void ChemkinOut::wrireaLump3(Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E)
{
	// write the lumped kinetic reaction in the form reag -> stoicCoeff[0] prod[1] + stoicCoeff[1] prod[2] + ...      A n E
	std::string reac;
	std::string species = "R" + std::to_string(reag.numberOfC());
	addUniqueName(&lumpedSpeciesList, species);
	reac += species + " => ";
	for (int i = 0; i < stoicCoeff.size(); i++)
	{
		if (prod[i + 1].isole()) species = "OLE" + std::to_string(prod[i + 1].numberOfC());
		else if (prod[i + 1].isCRad()) species = "R" + std::to_string(prod[i + 1].numberOfC());
		else
		{
			std::cerr << "Error: product of betadecomposition shold be either olefins or C * radical!" << std::endl;
			species = "ERR" + std::to_string(prod[i + 1].numberOfC());
		}
		addUniqueName(&lumpedSpeciesList, species);
		reac += std::to_string(stoicCoeff[i]) + " " + species;
		if (i != stoicCoeff.size() - 1) reac += " + ";
	}
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	reac = reac + con + "\n";

	if (reac.size() < CHEMKIN_MAX_CHAR_PER_ROW + 2)  // if the reaction is short enough to be read by chemkin print it directly
	{
		lumpedReactions << reac;
	}
	else											// else eliminate all the spaces and round the stoic coeff until it is short enough
	{
		lumpedReactions << "!" << reac;		// write the complete reaction as a comment
		bool sizeReached = false;
		int numOfDigit = 6;
		while (!sizeReached)
		{
			reac.clear();
			std::vector<double> roundedStoich = roundStoicCoeff(&stoicCoeff, numOfDigit);
			species = "R" + std::to_string(reag.numberOfC());
			addUniqueName(&lumpedSpeciesList, species);
			reac += species + "=>";
			for (int i = 0; i < stoicCoeff.size(); i++)
			{
				if (prod[i + 1].isole()) species = "OLE" + std::to_string(prod[i + 1].numberOfC());
				else if (prod[i + 1].isCRad()) species = "R" + std::to_string(prod[i + 1].numberOfC());
				else
				{
					std::cerr << "Error: product of betadecomposition shold be either olefins or C * radical!" << std::endl;
					species = "ERR" + std::to_string(prod[i + 1].numberOfC());
				}
				addUniqueName(&lumpedSpeciesList, species);
				reac += toString(roundedStoich[i], numOfDigit) + species;
				if (i != stoicCoeff.size() - 1) reac += "+";
			}
			char con[50];
			sprintf(con, " %.3e %.3f %.0f ", A, n, E);
			reac = reac + con + "\n";

			if (reac.size() < CHEMKIN_MAX_CHAR_PER_ROW + 2)
			{
				sizeReached = true;
			}
			else numOfDigit--;
			if (numOfDigit == 1)
			{
				std::cerr << "WARNING: In order to shorten the radical betadecomposition reaction under the maximum length admitted by chemkin, a number of significative digit less than 2 is required!" << std::endl;
				break;
			}
		}
		lumpedReactions << reac;

	}

	if (stoicCoeff.size() + 1 > max_species) max_species = stoicCoeff.size() + 1;
}

void ChemkinOut::wrireaLump4(int numC, double A, double n, double E)     // RnumC + O2 -> RnumCOO
{
	std::string species = "R" + std::to_string(numC);
	addUniqueName(&lumpedSpeciesList, species);
	addUniqueName(&lumpedSpeciesList, "O2");
	lumpedReactions << species << " + " << "O2" << " => ";
	species = "R" + std::to_string(numC) + "OO";
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	//lumpedReactions<< "R" << numC << " + " << "O2" << " -> " << " R" << numC << "OO";
	//lumpedReactions<< "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	//lumpedReactions<< std::endl;
}

void ChemkinOut::wrireaLump5(int numC, double A, double n, double E)     // RnumCOO -> RnumC + O2
{
	std::string species = "R" + std::to_string(numC) + "OO";
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " => ";
	species = "R" + std::to_string(numC);
	addUniqueName(&lumpedSpeciesList, species);
	addUniqueName(&lumpedSpeciesList, "O2");
	lumpedReactions << species << " + O2";
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	//lumpedReactions<< "R" << numC << "OO" << " -> " << "R" << numC << " + " << " O2";
	//lumpedReactions<< "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	//lumpedReactions<< std::endl;
}

void ChemkinOut::wrireaLump6(int numC, double A, double n, double E)     // RnumC + O2 -> OLEnumC + HO2
{
	std::string species = "R" + std::to_string(numC);
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " + O2 => ";
	species = "OLE" + std::to_string(numC);
	addUniqueName(&lumpedSpeciesList, species);
	addUniqueName(&lumpedSpeciesList, "O2");
	addUniqueName(&lumpedSpeciesList, "HO2");
	lumpedReactions << species << " + HO2";
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	//lumpedReactions<< "R" << numC << " + O2 -> " << "OLE" << numC << " + " << " HO2";
	//lumpedReactions<< "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	//lumpedReactions<< std::endl;
}

void ChemkinOut::wrireaLump7(int numC, double A, double n, double E)     // RnumCOO -> QnumCOOH
{
	std::string species = "R" + std::to_string(numC) + "OO";
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " => ";
	species = "Q" + std::to_string(numC) + "OOH";
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species;
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	//lumpedReactions<< "R" << numC << "OO -> " << "Q" << numC << "OOH";
	//lumpedReactions<< "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	//lumpedReactions<< std::endl;
}

void ChemkinOut::wrireaLump8(int numC, double A, double n, double E)     // QnumCOOH -> RnumCOO
{
	std::string species = "Q" + std::to_string(numC) + "OOH";
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " => ";
	species = "R" + std::to_string(numC) + "OO";
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species;
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	//lumpedReactions<< "Q" << numC << "OOH -> " << "R" << numC << "OO";
	//lumpedReactions<< "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	//lumpedReactions<< std::endl;
}

void ChemkinOut::wrireaLump9(int numC, double A, double n, double E)     // Q'numC'OOH -> ETER'numC' + OH
{
	std::string species = "Q" + std::to_string(numC) + "OOH";
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " => ";
	species = "ETER" + std::to_string(numC);
	addUniqueName(&lumpedSpeciesList, species);
	addUniqueName(&lumpedSpeciesList, "OH");
	lumpedReactions << species << " + OH";
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	//lumpedReactions<< "Q" << numC << "OOH -> " << "ETER" << numC << " + OH";
	//lumpedReactions<< "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	//lumpedReactions<< std::endl;
}

void ChemkinOut::wrireaLump10(int numC, double A, double n, double E)     // Q'numC'OOH -> OLE'numC' + HO2
{
	std::string species = "Q" + std::to_string(numC) + "OOH";
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " => ";
	species = "OLE" + std::to_string(numC);
	addUniqueName(&lumpedSpeciesList, species);
	addUniqueName(&lumpedSpeciesList, "HO2");
	lumpedReactions << species << " + HO2";
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	//lumpedReactions<< "Q" << numC << "OOH -> " << "OLE" << numC << " + HO2";
	//lumpedReactions<< "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	//lumpedReactions<< std::endl;
}

void ChemkinOut::wrireaLump11(Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E)
{
	std::string reac;
	std::string species = "Q" + std::to_string(reag.numberOfC()) + "OOH";
	addUniqueName(&lumpedSpeciesList, species);
	reac += species + " => 1.0 OH + ";
	addUniqueName(&lumpedSpeciesList, "OH");
	for (int i = 0; i < stoicCoeff.size(); i++)
	{
		if (prod[i + 1].isole()) species = "OLE" + std::to_string(prod[i + 1].numberOfC());
		else if (prod[i + 1].ischeto()) species = "C" + std::to_string(prod[i + 1].numberOfC()) + "HO";
		else
		{
			std::cerr << "Error: product of QOOH beta decomposition shold be either olefins or C* radical!" << std::endl;
			species = "ERR" + std::to_string(prod[i + 1].numberOfC());
		}
		addUniqueName(&lumpedSpeciesList, species);
		reac += std::to_string(stoicCoeff[i]) + " " + species;
		if (i != stoicCoeff.size() - 1) reac += " + ";
	}
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	reac = reac + con + "\n";


	if (reac.size() < CHEMKIN_MAX_CHAR_PER_ROW + 2)  // if the reaction is short enough to be read by chemkin print it directly
	{
		lumpedReactions << reac;
	}
	else											// else eliminate all the spaces and round the stoic coeff until it is short enough
	{
		lumpedReactions << "!" << reac;		// write the complete reaction as a comment
		bool sizeReached = false;
		int numOfDigit = 6;
		while (!sizeReached)
		{
			reac.clear();
			std::vector<double> roundedStoich = roundStoicCoeff(&stoicCoeff, numOfDigit);
			species = "Q" + std::to_string(reag.numberOfC()) + "OOH";
			addUniqueName(&lumpedSpeciesList, species);
			reac += species + "=>1.OH+";
			addUniqueName(&lumpedSpeciesList, "OH");
			for (int i = 0; i < stoicCoeff.size(); i++)
			{
				if (prod[i + 1].isole()) species = "OLE" + std::to_string(prod[i + 1].numberOfC());
				else if (prod[i + 1].ischeto()) species = "C" + std::to_string(prod[i + 1].numberOfC()) + "HO";
				else
				{
					std::cerr << "Error: product of QOOH beta decomposition shold be either olefins or C* radical!" << std::endl;
					species = "ERR" + std::to_string(prod[i + 1].numberOfC());
				}
				addUniqueName(&lumpedSpeciesList, species);
				reac += toString(roundedStoich[i], numOfDigit) + species;
				if (i != stoicCoeff.size() - 1) reac += "+";
			}
			char con[50];
			sprintf(con, " %.3e %.3f %.0f ", A, n, E);
			reac = reac + con + "\n";

			if (reac.size() < CHEMKIN_MAX_CHAR_PER_ROW + 2)
			{
				sizeReached = true;
			}
			else numOfDigit--;
			if (numOfDigit == 1)
			{
				std::cerr << "WARNING: In order to shorten the QOOH beta decomposition reaction under the maximum length admitted by chemkin, a number of significative digit less than 2 is required!" << std::endl;
				break;
			}
		}
		lumpedReactions << reac;

	}
	if (stoicCoeff.size() + 2 > max_species) max_species = stoicCoeff.size() + 2;
}

void ChemkinOut::wrireaLump11b(Molecola reag, Molecola prod[], std::vector<double> stoicCoeff1, double A1, double n1, double E1,
	std::vector<double> stoicCoeff2, double A2, double n2, double E2)
{
	wrireaLump11(reag, prod, stoicCoeff1, A1, n1, E1);
	lumpedReactions << "DUP" << std::endl;
	wrireaLump11(reag, prod, stoicCoeff2, A2, n2, E2);
	lumpedReactions << "DUP" << std::endl;
}

void ChemkinOut::wrireaLump12(int numC, double A, double n, double E)     // Q'numC'OOH + O2 -> OOQ'numC'OOH
{
	std::string species = "Q" + std::to_string(numC) + "OOH";
	addUniqueName(&lumpedSpeciesList, species);
	addUniqueName(&lumpedSpeciesList, "O2");
	lumpedReactions << species << " + O2 => ";
	species = "OOQ" + std::to_string(numC) + "OOH";
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species;
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	//lumpedReactions<< "Q" << numC << "OOH + O2-> " << "OOQ" << numC << "OOH";
	//lumpedReactions<< "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	//lumpedReactions<< std::endl;
}
void ChemkinOut::wrireaLump13(int numC, double A, double n, double E)     // OOQ'numC'OOH -> Q'numC'OOH + O2
{
	std::string species = "OOQ" + std::to_string(numC) + "OOH";
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " => ";
	species = "Q" + std::to_string(numC) + "OOH";
	addUniqueName(&lumpedSpeciesList, "O2");
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " + O2";
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	//lumpedReactions<< "OOQ" << numC << "OOH -> " << "Q" << numC << "OOH + O2";
	//lumpedReactions<< "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	//lumpedReactions<< std::endl;
}

void ChemkinOut::wrireaLump14(int numC, double A, double n, double E)     // OOQ'numC'OOH -> OQ'numC'OOH + OH
{
	std::string species = "OOQ" + std::to_string(numC) + "OOH";
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " => ";
	species = "OQ" + std::to_string(numC) + "OOH";
	addUniqueName(&lumpedSpeciesList, "OH");
	addUniqueName(&lumpedSpeciesList, species);
	lumpedReactions << species << " + OH";
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	lumpedReactions << con << std::endl;
	//lumpedReactions<< "OOQ" << numC << "OOH -> " << "OQ" << numC << "OOH + OH";
	//lumpedReactions<< "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	//lumpedReactions<< std::endl;
}

void ChemkinOut::wrireaLump15(Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E)
{
	std::string reac;
	std::string species = "OQ" + std::to_string(reag.numberOfC()) + "OOH";
	addUniqueName(&lumpedSpeciesList, species);
	reac += species + " => 1.0 OH + ";
	addUniqueName(&lumpedSpeciesList, "OH");
	for (int i = 0; i < stoicCoeff.size(); i++)
	{
		if (prod[i + 1].isRO()) species = "R" + std::to_string(prod[i + 1].numberOfC()) + "O";
		else if (prod[i + 1].numCrad() == 0 && prod[i + 1].numKeto() == 1) species = "C" + std::to_string(prod[i + 1].numberOfC()) + "HO";
		else
		{
			std::cerr << "Error: product of ketohydroperoxides decomposition shold be either CHO or RO radical!" << std::endl;
			species = "ERR" + std::to_string(prod[i + 1].numberOfC());
		}
		addUniqueName(&lumpedSpeciesList, species);
		reac += std::to_string(stoicCoeff[i]) + " " + species;
		if (i != stoicCoeff.size() - 1) reac += " + ";
	}
	char con[50];
	sprintf(con, "   %9.3e %8.3f %8.0f ", A, n, E);
	reac = reac + con + "\n";

	if (reac.size() < CHEMKIN_MAX_CHAR_PER_ROW + 2)  // if the reaction is short enough to be read by chemkin print it directly
	{
		lumpedReactions << reac;
	}
	else											// else eliminate all the spaces and round the stoic coeff until it is short enough
	{
		lumpedReactions << "!" << reac;		// write the complete reaction as a comment
		bool sizeReached = false;
		int numOfDigit = 6;
		while (!sizeReached)
		{
			reac.clear();
			std::vector<double> roundedStoich = roundStoicCoeff(&stoicCoeff, numOfDigit);

			species = "OQ" + std::to_string(reag.numberOfC()) + "OOH";
			addUniqueName(&lumpedSpeciesList, species);
			reac += species + "=>1.OH+";
			addUniqueName(&lumpedSpeciesList, "OH");
			for (int i = 0; i < stoicCoeff.size(); i++)
			{
				if (prod[i + 1].isRO()) species = "R" + std::to_string(prod[i + 1].numberOfC()) + "O";
				else if (prod[i + 1].numCrad() == 0 && prod[i + 1].numKeto() == 1) species = "C" + std::to_string(prod[i + 1].numberOfC()) + "HO";
				else
				{
					std::cerr << "Error: product of ketohydroperoxides decomposition shold be either CHO or RO radical!" << std::endl;
					species = "ERR" + std::to_string(prod[i + 1].numberOfC());
				}
				addUniqueName(&lumpedSpeciesList, species);
				reac += toString(roundedStoich[i], numOfDigit) + species;
				if (i != stoicCoeff.size() - 1) reac += "+";
			}
			char con[50];
			sprintf(con, " %.3e %.3f %.0f ", A, n, E);
			reac = reac + con + "\n";

			if (reac.size() < CHEMKIN_MAX_CHAR_PER_ROW + 2)
			{
				sizeReached = true;
			}
			else numOfDigit--;
			if (numOfDigit == 1)
			{
				std::cerr << "WARNING: In order to shorten the ketohydroperoxides decomposition  reaction under the maximum length admitted by chemkin, a number of significative digit less than 2 is required!" << std::endl;
				break;
			}
		}
		lumpedReactions << reac;

	}
	if (stoicCoeff.size() + 2 > max_species) max_species = stoicCoeff.size() + 2;
}


std::string ChemkinOut::nameHAbsRad(HAbsRad rad)
{
	switch (rad)
	{
	case o2:
		return "O2";
		break;
	case oh:
		return "OH";
		break;
	case h:
		return "H";
		break;
	case o:
		return "O";
		break;
	case ho2:
		return "HO2";
		break;
	case ch3:
		return "CH3";
		break;
	case c2h5:
		return "C2H5";
		break;
	case ch3oo:
		return "CH3O2";
		break;
	default:
		std::cerr << "ERROR: " << rad << " is not a valid value for HAbsRad radical!" << std::endl;
		return "ERR.";
		break;
	}
}

std::string ChemkinOut::nameHAbsRadPlusH(HAbsRad rad)
{
	switch (rad)
	{
	case o2:
		return "HO2";
		break;
	case oh:
		return "H2O";
		break;
	case h:
		return "H2";
		break;
	case o:
		return "OH";
		break;
	case ho2:
		return "H2O2";
		break;
	case ch3:
		return "CH4";
		break;
	case c2h5:
		return "C2H6";
		break;
	case ch3oo:
		return "CH3O2H";
		break;
	default:
		std::cerr << "ERROR: " << rad << " is not a valid value for HAbsRad radical!" << std::endl;
		return "ERR.";
		break;
	}
}

void ChemkinOut::printLumpName(std::string name)				// print the name of the reaction to lumpedReactions
{
	lumpedReactions << std::endl;
	lumpedReactions << name << std::endl;
}

void ChemkinOut::makeThermoOut(std::string path, std::vector<Molecola*>* specList)
{
	ofstream thermoOut;   // cumulated selectivity output file
	thermoOut.open(path);
	if (!thermoOut)
	{
		cerr << "Error in opening the thermo out file!" << endl << endl;
		exit(1);
	};

	std::cout << "Printing thermo file ..." << std::endl;

	// print special molecules
	thermoOut << "! WARNING: THE PROPERTIES IN THIS FILE ARE NOT CORRECT                          " << std::endl;
	thermoOut << "!          ONLY SUITABLE FOR CSTR SIMULATIONS OR OTHER                          " << std::endl;
	thermoOut << "!          SIMULATIONS WHERE THE THERMO DATA DOES NOT                           " << std::endl;
	thermoOut << "!          AFFECT THE RESULT!!!                                                 " << std::endl;
	thermoOut << "THERMO ALL                                                                      " << std::endl;
	thermoOut << "300.   1000.   3500.                                                            " << std::endl;
	thermoOut << "O                       O   1               G    300.00   3500.00  950.00      1" << std::endl;
	thermoOut << " 2.57318360e+00-8.95609984e-05 4.05096303e-08-8.39812674e-12 9.43621991e-16    2" << std::endl;
	thermoOut << " 2.92191409e+04 4.74952023e+00 2.95200330e+00-1.68459131e-03 2.55897854e-06    3" << std::endl;
	thermoOut << "-1.77574473e-09 4.66034833e-13 2.91471652e+04 2.94136507e+00                   4" << std::endl;
	thermoOut << "O2                      O   2               G    300.00   3500.00  760.00      1" << std::endl;
	thermoOut << " 2.81750648e+00 2.49838007e-03-1.52493521e-06 4.50547608e-10-4.87702792e-14    2" << std::endl;
	thermoOut << "-9.31713392e+02 7.94729337e+00 3.46035080e+00-8.85011121e-04 5.15281056e-06    3" << std::endl;
	thermoOut << "-5.40712413e-09 1.87809542e-12-1.02942573e+03 5.02236126e+00                   4" << std::endl;
	thermoOut << "N2                      N   2               G    300.00   3500.00 1050.00      1" << std::endl;
	thermoOut << " 2.71287897e+00 1.90359754e-03-8.54297556e-07 1.84170938e-10-1.54715988e-14    2" << std::endl;
	thermoOut << "-8.40225273e+02 7.15926558e+00 3.85321336e+00-2.44053349e-03 5.35160392e-06    3" << std::endl;
	thermoOut << "-3.75608397e-09 9.22684330e-13-1.07969550e+03 1.60217419e+00                   4" << std::endl;
	thermoOut << "H                       H   1               G    300.00   3500.00 1490.00      1" << std::endl;
	thermoOut << " 2.50000000e+00 7.40336223e-15-5.56967416e-18 1.73924876e-21-1.92673709e-25    2" << std::endl;
	thermoOut << " 2.54716200e+04-4.60117600e-01 2.50000000e+00-4.07455160e-15 5.98527266e-18    3" << std::endl;
	thermoOut << "-3.43074982e-21 6.74775716e-25 2.54716200e+04-4.60117600e-01                   4" << std::endl;
	thermoOut << "H2                      H   2               G    300.00   3500.00  750.00      1" << std::endl;
	thermoOut << " 3.73110902e+00-8.86706214e-04 1.12286897e-06-3.74349782e-10 4.17963674e-14    2" << std::endl;
	thermoOut << "-1.08851547e+03-5.35285855e+00 3.08866003e+00 2.53968841e-03-5.72992027e-06    3" << std::endl;
	thermoOut << " 5.71701843e-09-1.98865970e-12-9.92148124e+02-2.43823459e+00                   4" << std::endl;
	thermoOut << "H2O                     H   2O   1          G    300.00   3500.00 1590.00      1" << std::endl;
	thermoOut << " 2.30940463e+00 3.65433887e-03-1.22983871e-06 2.11931683e-10-1.50333493e-14    2" << std::endl;
	thermoOut << "-2.97294901e+04 8.92765177e+00 4.03530937e+00-6.87559833e-04 2.86629214e-06    3" << std::endl;
	thermoOut << "-1.50552360e-09 2.55006790e-13-3.02783278e+04-1.99201641e-01                   4" << std::endl;
	thermoOut << "H2O2                    H   2O   2          G    300.00   3500.00 1180.00      1" << std::endl;
	thermoOut << " 4.56163072e+00 4.35560969e-03-1.48694629e-06 2.38275424e-10-1.46610352e-14    2" << std::endl;
	thermoOut << "-1.80016693e+04 5.66597119e-01 2.91896355e+00 9.92397296e-03-8.56537418e-06    3" << std::endl;
	thermoOut << " 4.23738723e-09-8.61930485e-13-1.76139998e+04 8.76340177e+00                   4" << std::endl;
	thermoOut << "CH3                     C   1H   3          G    300.00   3500.00 1270.00      1" << std::endl;
	thermoOut << " 2.57723974e+00 6.62601164e-03-2.54906392e-06 4.67320141e-10-3.34867663e-14    2" << std::endl;
	thermoOut << " 1.65488693e+04 6.94195966e+00 3.53327401e+00 3.61488008e-03 1.00739068e-06    3" << std::endl;
	thermoOut << "-1.39958516e-09 3.34014277e-13 1.63060366e+04 2.10113860e+00                   4" << std::endl;
	thermoOut << "HO2                     H   1O   2          G    300.00   3500.00 1540.00      1" << std::endl;
	thermoOut << " 4.16318067e+00 1.99798265e-03-4.89192086e-07 7.71153172e-11-7.30772104e-15    2" << std::endl;
	thermoOut << " 4.41348948e+01 2.95517985e+00 2.85241381e+00 5.40257188e-03-3.80535043e-06    3" << std::endl;
	thermoOut << " 1.51268170e-09-2.40354212e-13 4.47851086e+02 9.84483831e+00                   4" << std::endl;
	thermoOut << "OH                      H   1O   1          G    300.00   3500.00  880.00      1" << std::endl;
	thermoOut << " 3.62538436e+00-5.02165281e-04 8.36958463e-07-2.95714531e-10 3.30350486e-14    2" << std::endl;
	thermoOut << " 3.41380110e+03 1.55419440e+00 3.37995109e+00 6.13440526e-04-1.06464235e-06    3" << std::endl;
	thermoOut << " 1.14489214e-09-3.76228211e-13 3.45699735e+03 2.70689352e+00                   4" << std::endl;
	thermoOut << "AR                      AR  1               G    300.00   3500.00 1490.00      1" << std::endl;
	thermoOut << " 2.50000000e+00 7.40336223e-15-5.56967416e-18 1.73924876e-21-1.92673709e-25    2" << std::endl;
	thermoOut << "-7.45375000e+02 4.36600000e+00 2.50000000e+00-4.07455160e-15 5.98527266e-18    3" << std::endl;
	thermoOut << "-3.43074982e-21 6.74775716e-25-7.45375000e+02 4.36600000e+00                   4" << std::endl;
	thermoOut << "HE                      HE  1               G    300.00   3500.00 1490.00      1" << std::endl;
	thermoOut << " 2.50000000e+00 7.40336223e-15-5.56967416e-18 1.73924876e-21-1.92673709e-25    2" << std::endl;
	thermoOut << "-7.45375000e+02 9.28723974e-01 2.50000000e+00-4.07455160e-15 5.98527266e-18    3" << std::endl;
	thermoOut << "-3.43074982e-21 6.74775716e-25-7.45375000e+02 9.28723974e-01                   4" << std::endl;

	// print all the regular molecules
	for (int i = 0; i < specList->size(); i++)
	{
		Molecola* mol = (*specList)[i];
		if (mol->isSpecialMolecule() != 0)
			continue;

		int numberOfC = mol->numberOfC();
		int numberOfH = mol->numberOfH();
		int numberOfO = mol->numberOfO();
		std::string name = molToName(*mol);

		thermoOut << std::left << std::setw(24) << name;

		if (numberOfC > 0)
			thermoOut << "C" << std::right << std::setw(4) << numberOfC;
		else
			thermoOut << "     ";

		if (numberOfH > 0)
			thermoOut << "H" << std::right << std::setw(4) << numberOfH;
		else
			thermoOut << "     ";

		if (numberOfO > 0)
			thermoOut << "O" << std::right << std::setw(4) << numberOfO;
		else
			thermoOut << "     ";

		thermoOut << "     G    300.00   3500.00 1790.00      1" << std::endl;
		thermoOut << " 2.73219585e+01 3.55870996e-02-1.24650774e-05 2.00089265e-09-1.21997682e-13    2" << std::endl;
		thermoOut << "-4.84634178e+04-1.12262056e+02 1.16215154e+00 9.40447688e-02-6.14519510e-05    3" << std::endl;
		thermoOut << " 2.02455383e-08-2.67013255e-12-3.90982069e+04 2.91745387e+01                   4" << std::endl;
	}

	std::cout << "                        ... COMPLETED" << std::endl;

}

void ChemkinOut::makeThermoOutLumped(std::string path)
{
	std::vector<std::string> names = lumpedNames;
	std::vector<Molecola*>* specList = &lumpedMolecules;


	ofstream thermoOut;   // cumulated selectivity output file
	thermoOut.open(path);
	if (!thermoOut)
	{
		cerr << "Error in opening the thermo out file!" << endl << endl;
		exit(1);
	};

	std::cout << "Printing thermo file ..." << std::endl;

	// print special molecules
	thermoOut << "! WARNING: THE PROPERTIES IN THIS FILE ARE NOT CORRECT                          " << std::endl;
	thermoOut << "!          ONLY SUITABLE FOR CSTR SIMULATIONS OR OTHER                          " << std::endl;
	thermoOut << "!          SIMULATIONS WHERE THE THERMO DATA DOES NOT                           " << std::endl;
	thermoOut << "!          AFFECT THE RESULT!!!                                                 " << std::endl;
	thermoOut << "THERMO ALL                                                                      " << std::endl;
	thermoOut << "300.   1000.   3500.                                                            " << std::endl;
	thermoOut << "O                       O   1               G    300.00   3500.00  950.00      1" << std::endl;
	thermoOut << " 2.57318360e+00-8.95609984e-05 4.05096303e-08-8.39812674e-12 9.43621991e-16    2" << std::endl;
	thermoOut << " 2.92191409e+04 4.74952023e+00 2.95200330e+00-1.68459131e-03 2.55897854e-06    3" << std::endl;
	thermoOut << "-1.77574473e-09 4.66034833e-13 2.91471652e+04 2.94136507e+00                   4" << std::endl;
	thermoOut << "O2                      O   2               G    300.00   3500.00  760.00      1" << std::endl;
	thermoOut << " 2.81750648e+00 2.49838007e-03-1.52493521e-06 4.50547608e-10-4.87702792e-14    2" << std::endl;
	thermoOut << "-9.31713392e+02 7.94729337e+00 3.46035080e+00-8.85011121e-04 5.15281056e-06    3" << std::endl;
	thermoOut << "-5.40712413e-09 1.87809542e-12-1.02942573e+03 5.02236126e+00                   4" << std::endl;
	thermoOut << "N2                      N   2               G    300.00   3500.00 1050.00      1" << std::endl;
	thermoOut << " 2.71287897e+00 1.90359754e-03-8.54297556e-07 1.84170938e-10-1.54715988e-14    2" << std::endl;
	thermoOut << "-8.40225273e+02 7.15926558e+00 3.85321336e+00-2.44053349e-03 5.35160392e-06    3" << std::endl;
	thermoOut << "-3.75608397e-09 9.22684330e-13-1.07969550e+03 1.60217419e+00                   4" << std::endl;
	thermoOut << "H                       H   1               G    300.00   3500.00 1490.00      1" << std::endl;
	thermoOut << " 2.50000000e+00 7.40336223e-15-5.56967416e-18 1.73924876e-21-1.92673709e-25    2" << std::endl;
	thermoOut << " 2.54716200e+04-4.60117600e-01 2.50000000e+00-4.07455160e-15 5.98527266e-18    3" << std::endl;
	thermoOut << "-3.43074982e-21 6.74775716e-25 2.54716200e+04-4.60117600e-01                   4" << std::endl;
	thermoOut << "H2                      H   2               G    300.00   3500.00  750.00      1" << std::endl;
	thermoOut << " 3.73110902e+00-8.86706214e-04 1.12286897e-06-3.74349782e-10 4.17963674e-14    2" << std::endl;
	thermoOut << "-1.08851547e+03-5.35285855e+00 3.08866003e+00 2.53968841e-03-5.72992027e-06    3" << std::endl;
	thermoOut << " 5.71701843e-09-1.98865970e-12-9.92148124e+02-2.43823459e+00                   4" << std::endl;
	thermoOut << "H2O                     H   2O   1          G    300.00   3500.00 1590.00      1" << std::endl;
	thermoOut << " 2.30940463e+00 3.65433887e-03-1.22983871e-06 2.11931683e-10-1.50333493e-14    2" << std::endl;
	thermoOut << "-2.97294901e+04 8.92765177e+00 4.03530937e+00-6.87559833e-04 2.86629214e-06    3" << std::endl;
	thermoOut << "-1.50552360e-09 2.55006790e-13-3.02783278e+04-1.99201641e-01                   4" << std::endl;
	thermoOut << "H2O2                    H   2O   2          G    300.00   3500.00 1180.00      1" << std::endl;
	thermoOut << " 4.56163072e+00 4.35560969e-03-1.48694629e-06 2.38275424e-10-1.46610352e-14    2" << std::endl;
	thermoOut << "-1.80016693e+04 5.66597119e-01 2.91896355e+00 9.92397296e-03-8.56537418e-06    3" << std::endl;
	thermoOut << " 4.23738723e-09-8.61930485e-13-1.76139998e+04 8.76340177e+00                   4" << std::endl;
	thermoOut << "CH3                     C   1H   3          G    300.00   3500.00 1270.00      1" << std::endl;
	thermoOut << " 2.57723974e+00 6.62601164e-03-2.54906392e-06 4.67320141e-10-3.34867663e-14    2" << std::endl;
	thermoOut << " 1.65488693e+04 6.94195966e+00 3.53327401e+00 3.61488008e-03 1.00739068e-06    3" << std::endl;
	thermoOut << "-1.39958516e-09 3.34014277e-13 1.63060366e+04 2.10113860e+00                   4" << std::endl;
	thermoOut << "HO2                     H   1O   2          G    300.00   3500.00 1540.00      1" << std::endl;
	thermoOut << " 4.16318067e+00 1.99798265e-03-4.89192086e-07 7.71153172e-11-7.30772104e-15    2" << std::endl;
	thermoOut << " 4.41348948e+01 2.95517985e+00 2.85241381e+00 5.40257188e-03-3.80535043e-06    3" << std::endl;
	thermoOut << " 1.51268170e-09-2.40354212e-13 4.47851086e+02 9.84483831e+00                   4" << std::endl;
	thermoOut << "OH                      H   1O   1          G    300.00   3500.00  880.00      1" << std::endl;
	thermoOut << " 3.62538436e+00-5.02165281e-04 8.36958463e-07-2.95714531e-10 3.30350486e-14    2" << std::endl;
	thermoOut << " 3.41380110e+03 1.55419440e+00 3.37995109e+00 6.13440526e-04-1.06464235e-06    3" << std::endl;
	thermoOut << " 1.14489214e-09-3.76228211e-13 3.45699735e+03 2.70689352e+00                   4" << std::endl;
	thermoOut << "AR                      AR  1               G    300.00   3500.00 1490.00      1" << std::endl;
	thermoOut << " 2.50000000e+00 7.40336223e-15-5.56967416e-18 1.73924876e-21-1.92673709e-25    2" << std::endl;
	thermoOut << "-7.45375000e+02 4.36600000e+00 2.50000000e+00-4.07455160e-15 5.98527266e-18    3" << std::endl;
	thermoOut << "-3.43074982e-21 6.74775716e-25-7.45375000e+02 4.36600000e+00                   4" << std::endl;
	thermoOut << "HE                      HE  1               G    300.00   3500.00 1490.00      1" << std::endl;
	thermoOut << " 2.50000000e+00 7.40336223e-15-5.56967416e-18 1.73924876e-21-1.92673709e-25    2" << std::endl;
	thermoOut << "-7.45375000e+02 9.28723974e-01 2.50000000e+00-4.07455160e-15 5.98527266e-18    3" << std::endl;
	thermoOut << "-3.43074982e-21 6.74775716e-25-7.45375000e+02 9.28723974e-01                   4" << std::endl;

	// print all the regular molecules
	for (int i = 0; i < specList->size(); i++)
	{
		Molecola* mol = (*specList)[i];
		if (mol->isSpecialMolecule() != 0)
			continue;

		int numberOfC = mol->numberOfC();
		int numberOfH = mol->numberOfH();
		int numberOfO = mol->numberOfO();
		std::string name = names[i];

		thermoOut << std::left << std::setw(24) << name;

		if (numberOfC > 0)
			thermoOut << "C" << std::right << std::setw(4) << numberOfC;
		else
			thermoOut << "     ";

		if (numberOfH > 0)
			thermoOut << "H" << std::right << std::setw(4) << numberOfH;
		else
			thermoOut << "     ";

		if (numberOfO > 0)
			thermoOut << "O" << std::right << std::setw(4) << numberOfO;
		else
			thermoOut << "     ";

		thermoOut << "     G    300.00   3500.00 1790.00      1" << std::endl;
		thermoOut << " 2.73219585e+01 3.55870996e-02-1.24650774e-05 2.00089265e-09-1.21997682e-13    2" << std::endl;
		thermoOut << "-4.84634178e+04-1.12262056e+02 1.16215154e+00 9.40447688e-02-6.14519510e-05    3" << std::endl;
		thermoOut << " 2.02455383e-08-2.67013255e-12-3.90982069e+04 2.91745387e+01                   4" << std::endl;
	}

	std::cout << "                        ... COMPLETED" << std::endl;

}

std::string ChemkinOut::checkIfNamed(Molecola* mol)
{
	std::string name = "not found";
	for (int i = 0; i < namedSpecies.size(); i++)
	{
		if (*mol == namedSpecies[i])
		{
			name = namesList[i];
			break;
		}
	}

	return name;
}

