
#include "chemkinOut.h"

ChemkinOut::ChemkinOut(CKMechReader* baseMechanism)
{
	baseMech = baseMechanism;
}

void ChemkinOut::removeH(std::string* name) // if it is present in the name, removes one hydrogen to a alkane name in order to make it a substitute group
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
				(*name) = name->substr(0, Hpos + 1);
				name->append(std::to_string(val));
			}
		}
	}
}

std::string ChemkinOut::nameCOBranch(Molecola branch)
{
	std::string name;
	std::vector<int> prodMainChain(0);
	std::vector<int> prodMePos(0);
	std::vector<int> prodEtPos(0);
	std::vector<int> prodMap(0);
	branch.orderedListOfMeAndEt(&prodMainChain, &prodMePos,
		&prodEtPos, &prodMap);
	int posR = prodMap[branch.posCrad() - 1];
	Molecola parProd = branch.parentFuel();
	name.append(alkaneStructurePrefix(&parProd));
	name.append("C");
	name.append(std::to_string(parProd.numberOfC()));
	if (posR != 1)
	{
		name.append("-");
		name.append(std::to_string(posR));
	}
	return name;
}

std::string ChemkinOut::molToName(Molecola mol)
{

	std::string name;
	species type = mol.kindOfSPecies();
	//std::cout << "1" << std::endl;
	if (type == unidentified_)
	{
		std::cerr << "Error in naming the molecule " << mol << ", it doesn't match any expected molecule type." << std::endl;
		std::cerr << "numberR = " << mol.numCrad()<< std::endl;
		std::cerr << "numberCOO = " << mol.numCOORad()<< std::endl;
		std::cerr << "numberCOOH = " << mol.numCOOH()<< std::endl;
		std::cerr << "numberEther = " << mol.numberEthero()<< std::endl;
		std::cerr << "numberOLE = " << mol.numberOle()<< std::endl;
		std::cerr << "numberCO = " << mol.numKeto()<< std::endl;
		std::cerr << "specMol = " << mol.isSpecialMolecule()<< std::endl;
		std::cerr << "numLinEth = " << mol.numLinEther()<< std::endl;
		std::cerr << "numCOR = " << mol.numCORad()<< std::endl;
		mol.printStruttura();
		mol.printDoppioO();
		exit(0);
	}
	
	std::string n = "not found";
	if (type != special_)
	{
		n = checkIfNamed(&mol);
	}
	
	if (n != "not found")
	{
		name = n;
		//std::cout << "4" << std::endl;
	}
	else
	{
		//std::cout << "5" << std::endl;
		


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
			case 11:
				name.append("HCCO");
				break;
			case 12:
				name.append("CH3O");
				break;
			case 13:
				name.append("CH3OH");
				break;
			case 14:
				name.append("CH3O2");
				break;
			case 15:
				name.append("CH3O2H");
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


			name.append(alkaneStructurePrefix(&mol));

			name.append("C");
			name.append(std::to_string(mol.numberOfC()));


			switch (type)
			{
			case R_: // molecule is R
			{
				name.append("-");
				name.append(std::to_string(map[mol.posCrad() - 1]));
				break;
			}

			case ROO_:	// molecule is ROO
			{
				name.append("-");
				name.append(std::to_string(map[mol.posCOOrad() - 1]));
				name.append("O2");
				break;
			}

			case oleOH_:
			{
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

				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 1));
				name.append("OH");
				name.append(std::to_string(map[mol.posCOH() - 1]));
				name.append("-");
				name.append(std::to_string(pos));
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}
			case QOOH_:	// molecule is QOOH
			{
				int posCOOH = map[mol.posCOOH() - 1];
				int posCrad = map[mol.posCrad() - 1];
				if (posCrad > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.tipo(mol.posCOOH(), 2);
					ref.tipo(mol.posCrad(), 4);
					if (ref == mol)
					{
						int temp = posCOOH;
						posCOOH = posCrad;
						posCrad = temp;
					}
				}
				name.append("OOH");
				name.append(std::to_string(posCOOH));
				name.append("-");
				name.append(std::to_string(posCrad));
				break;
			}

			case OOQOOH_:	// molecule is OOQOOH
			{
				int posCOOH = map[mol.posCOOH() - 1];
				int posCOOrad = map[mol.posCOOrad() - 1];
				if (posCOOH > posCOOrad)
				{
					Molecola ref = mol.parentFuel();
					ref.tipo(mol.posCOOH(),3);
					ref.tipo(mol.posCOOrad(),4);
					if (ref == mol)
					{
						int temp = posCOOH;
						posCOOH = posCOOrad;
						posCOOrad = temp;
					}
				}
				name.append("OOH");
				name.append(std::to_string(posCOOH));
				name.append("-");
				name.append(std::to_string(posCOOrad));
				name.append("O2");
				break;
			}

			case OLE_:	// molecule is olefin
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
				break;
			}

			case CO_:	// molecule is CO
			{
				if (mol.isAldehyde())
				{
					if (mol.isLinear())
					{
						name.append("C");
						name.append(std::to_string(mol.numberOfC()-1));
						name.append("H");
						name.append(std::to_string(mol.numberOfH() - 1));
						name.append("CHO");
					}
					else
					{
						Molecola par = mol.parentFuel();
						Molecola prod1;
						Molecola prod2;
						int posKeto = mol.posKeto();
						int alpha[SIZEMAX + 1];
						mol.scorri(posKeto, 1, alpha);
						int posAlpha = 2;
						for (int i = 0; i < SIZEMAX + 1; i++)
						{
							if (alpha[i] == 1)
							{
								posAlpha = i;
								break;
							}
						}
						par.spezza(posKeto, posAlpha, &prod1, &prod2);
						Molecola Prod;
						if (prod1.size() > prod2.size())
							Prod = prod1;
						else
							Prod = prod2;
						name.append(nameCOBranch(Prod));
						name.append("-CHO");
					}
				}
				else
				{
					int posKeto = mol.posKeto();
					int alpha[SIZEMAX + 1];
					mol.scorri(posKeto, 1, alpha);
					int posAlpha1 = -1;
					int posAlpha2 = -1;
					for (int i = 0; i < SIZEMAX + 1; i++)
					{
						if (alpha[i] == 1 && posAlpha1 == -1)
						{
							posAlpha1 = i;
							continue;
						}
						if (alpha[i] == 1 && posAlpha1 != -1)
						{
							posAlpha2 = i;
							break;
						}
					}
					Molecola prod1, prod2, buff;
					Molecola par = mol.parentFuel();
					par.spezza(posKeto, posAlpha1, &buff, &prod1);
					par.spezza(posKeto, posAlpha2, &buff, &prod2);
					std::string nameProd1 = nameCOBranch(prod1);
					std::string nameProd2 = nameCOBranch(prod2);
					if (map[posAlpha1 - 1] < map[posAlpha2 - 1])
					{
						name.append(nameProd1);
						name.append("CO");
						name.append(nameProd2);
					}
					else
					{
						name.append(nameProd2);
						name.append("CO");
						name.append(nameProd1);
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
			case OHcEth_:	// molecule is an ether
			{
				name.append("H");
				name.append(std::to_string(mol.numberOfH()-1));
				name.append("OH");
				name.append(std::to_string(map[mol.posCOH() - 1]));

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
				int posCrad = map[mol.posCrad() - 1];
				int posKeto = map[mol.posKeto() - 1];
				if (posCrad > posKeto)
				{
					Molecola ref = mol.parentFuel();
					ref.addcheto(mol.posCrad());
					ref.tipo(mol.posKeto(), 2);
					if (ref == mol)
					{
						int temp = posCrad;
						posCrad = posKeto;
						posKeto = temp;
					}
				}
				name.append("H");
				name.append(std::to_string(mol.numberOfH()));
				name.append("-");
				name.append(std::to_string(posCrad));
				name.append("-");
				name.append(std::to_string(posKeto));
				name.append("O");
				break;
			}

			case KHP_:	// molecule is KHP
			{
				int posCOOH = map[mol.posCOOH() - 1];
				int posKeto = map[mol.posKeto() - 1];
				if (posKeto > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.addcheto(mol.posCOOH());
					ref.tipo(mol.posKeto(), 4);
					if (ref == mol)
					{
						int temp = posCOOH;
						posCOOH = posKeto;
						posKeto = temp;
					}
				}
				name.append("KET");
				name.append(std::to_string(posKeto));
				name.append("-");
				name.append(std::to_string(posCOOH));
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
				std::vector<int> posAtsOOH = mol.posOOHinPOOH2();
				int posAtOOH1 = posAtsOOH[0];
				int posAtOOH2 = posAtsOOH[1];
				int posOOH1 = map[posAtOOH1 - 1];
				int posOOH2 = map[posAtOOH2 - 1];
				int posC = map[mol.posCrad() - 1];
				if (posC > posOOH1)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(posAtOOH2, mol.posCrad());
					ref.tipo(posAtOOH1, 2);
					if (ref == mol)
					{
						int temp = posOOH1;
						posOOH1 = posC;
						posC = temp;
					}
				}
				if (posC > posOOH2)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(posAtOOH1, mol.posCrad());
					ref.tipo(posAtOOH2, 2);
					if (ref == mol)
					{
						int temp = posOOH2;
						posOOH2 = posC;
						posC = temp;
					}
				}
				name.append("Q");
				name.append(std::to_string(std::min(posOOH1, posOOH2)));
				name.append(std::to_string(std::max(posOOH1, posOOH2)));
				name.append("-");
				name.append(std::to_string(posC));
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
				int atPos;
				int atOther;
				if (pos1 <= mainChain.size() && pos2 <= mainChain.size())
				{
					if (pos1 > pos2)
					{
						pos = pos2;
						atPos = posOle[1];
						atOther = posOle[0];
					}
					else
					{
						pos = pos1;
						atPos = posOle[0];
						atOther = posOle[1];
					}

				}
				else if (pos1 <= mainChain.size() && pos2 > mainChain.size())
				{
					pos = pos2;
					atPos = posOle[1];
					atOther = posOle[0];
				}
				else if (pos1 > mainChain.size() && pos2 <= mainChain.size())
				{
					pos = pos1;
					atPos = posOle[0];
					atOther = posOle[1];
				}
				else if (pos1 > mainChain.size() && pos2 > mainChain.size())
				{
					if (pos1 < pos2)
					{
						pos = pos2;
						atPos = posOle[1];
						atOther = posOle[0];
					}
					else
					{
						pos = pos1;
						atPos = posOle[0];
						atOther = posOle[1];
					}
				}
				
				int posCOOH = map[mol.posCOOH() - 1];
				if (pos > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.addole(atOther, mol.posCOOH());
					ref.tipo(atPos, 4);
					if (ref == mol)
					{
						int temp = posCOOH;
						posCOOH = pos;
						pos = std::min(map[atOther - 1], temp);
					}
				}
				name.append("D");
				name.append(std::to_string(pos));
				name.append("-");
				name.append(std::to_string(posCOOH));
				name.append("OOH");
				break;
			}

			case cEthOOH_:	// molecule is an ether-OOH
			{

				name.append("O");
				std::vector<int> pos = mol.posEthero();
				//int pos1 = std::min(map[pos[0] - 1], map[pos[1] - 1]);
				//int pos2 = std::max(map[pos[0] - 1], map[pos[1] - 1]);
				int posO1;
				int posO2;
				int numAtO1;
				int numAtO2;
				if (map[pos[0] - 1] > map[pos[1] - 1])
				{
					numAtO1 = pos[1];
					posO1 = map[pos[1] - 1];
					numAtO2 = pos[0];
					posO2 = map[pos[0] - 1];
				}
				else
				{
					numAtO1 = pos[0];
					posO1 = map[pos[0] - 1];
					numAtO2 = pos[1];
					posO2 = map[pos[1] - 1];
				}
				int posCOOH = map[mol.posCOOH() - 1];
				if (posO1 > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(numAtO2, mol.posCOOH());
					ref.tipo(numAtO1, 4);
					if (ref == mol)
					{
						int temp = posO1;
						posO1 = posCOOH;
						posCOOH = temp;
					}
				}
				if (posO2 > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(numAtO1, mol.posCOOH());
					ref.tipo(numAtO2, 4);
					if (ref == mol)
					{
						int temp = posO2;
						posO2 = posCOOH;
						posCOOH = temp;
					}
				}
				name.append(std::to_string(std::min(posO1,posO2)));
				name.append(std::to_string(std::max(posO1, posO2)));
				name.append("-");
				name.append(std::to_string(posCOOH));
				name.append("OOH");
				break;
			}

			case olecEthOOH_:	// molecule is an ether-OOH
			{
				name.append("H");
				name.append(std::to_string(mol.numberOfH()-1));


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

				std::vector<int> pose = mol.posEthero();
				int pose1 = std::min(map[pose[0] - 1], map[pose[1] - 1]);
				int pose2 = std::max(map[pose[0] - 1], map[pose[1] - 1]);


				name.append(std::to_string(pos));

				name.append("O");
				std::vector<int> posl = mol.posEthero();
				//int pos1 = std::min(map[pos[0] - 1], map[pos[1] - 1]);
				//int pos2 = std::max(map[pos[0] - 1], map[pos[1] - 1]);
				int posO1;
				int posO2;
				int numAtO1;
				int numAtO2;
				if (map[posl[0] - 1] > map[posl[1] - 1])
				{
					numAtO1 = posl[1];
					posO1 = map[posl[1] - 1];
					numAtO2 = posl[0];
					posO2 = map[posl[0] - 1];
				}
				else
				{
					numAtO1 = posl[0];
					posO1 = map[posl[0] - 1];
					numAtO2 = posl[1];
					posO2 = map[posl[1] - 1];
				}
				int posCOOH = map[mol.posCOOH() - 1];
				if (posO1 > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(numAtO2, mol.posCOOH());
					ref.tipo(numAtO1, 4);
					if (ref == mol)
					{
						int temp = posO1;
						posO1 = posCOOH;
						posCOOH = temp;
					}
				}
				if (posO2 > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(numAtO1, mol.posCOOH());
					ref.tipo(numAtO2, 4);
					if (ref == mol)
					{
						int temp = posO2;
						posO2 = posCOOH;
						posCOOH = temp;
					}
				}
				name.append(std::to_string(std::min(posO1, posO2)));
				name.append(std::to_string(std::max(posO1, posO2)));
				name.append("-");
				name.append(std::to_string(posCOOH));
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
				int posR = map[mol.posCrad() - 1];
				if (pos1 <= mainChain.size() && pos2 <= mainChain.size())
					pos = std::min(pos1, pos2);
				else if (pos1 <= mainChain.size() && pos2 > mainChain.size())
					pos = pos2;
				else if (pos1 > mainChain.size() && pos2 <= mainChain.size())
					pos = pos1;
				else if (pos1 > mainChain.size() && pos2 > mainChain.size())
					pos = std::max(pos1, pos2);

				// try isomer
				Molecola copy = mol;
				std::vector<int> posOleCopy = copy.posOle(1);
				int posOle1Copy = posOleCopy[0];
				int posOle2Copy = posOleCopy[1];
				int posAtRCopy = copy.trova(2);
				if (copy.dist(posOle1Copy, posAtRCopy) == 1)
				{
					copy.removeOle(posOle1Copy, posOle2Copy);
					copy.Crad_to_C(posAtRCopy);
					copy.addole(posOle1Copy, posAtRCopy);
					copy.C_to_Crad(posOle2Copy);
				}
				else if (copy.dist(posOle2Copy, posAtRCopy) == 1)
				{
					copy.removeOle(posOle1Copy, posOle2Copy);
					copy.Crad_to_C(posAtRCopy);
					copy.addole(posOle2Copy, posAtRCopy);
					copy.C_to_Crad(posOle1Copy);
				}
				posOleCopy = copy.posOle(1);
				int pos1Copy = map[posOleCopy[0] - 1];
				int pos2Copy = map[posOleCopy[1] - 1];
				int posCopy = 0;
				int posRCopy = map[copy.posCrad() - 1];
				if (pos1Copy <= mainChain.size() && pos2Copy <= mainChain.size())
					posCopy = std::min(pos1Copy, pos2Copy);
				else if (pos1Copy <= mainChain.size() && pos2Copy > mainChain.size())
					posCopy = pos2Copy;
				else if (pos1Copy > mainChain.size() && pos2Copy <= mainChain.size())
					posCopy = pos1Copy;
				else if (pos1Copy > mainChain.size() && pos2Copy > mainChain.size())
					posCopy = std::max(pos1Copy, pos2Copy);

				if (posCopy < pos)
				{
					pos = posCopy;
					posR = posRCopy;
				}

				//name.append("D"); JIAXIN
				name.append("H");
				name.append(std::to_string(mol.numberOfH()));
				name.append(std::to_string(pos));
				name.append("-");
				name.append(std::to_string(posR));
				// JIAXIN
				//if ((mol.posCrad() == posOle[0] || mol.posCrad() == posOle[1]))
				//	name.append("a");
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
				int atPos;
				int atOther;
				int posKeto = map[mol.posKeto() - 1];
				if (pos1 <= mainChain.size() && pos2 <= mainChain.size())
				{
					//pos = std::min(pos1, pos2);
					if (pos1 < pos2)
					{
						pos = pos1;
						atPos = posOle[0];
						atOther = posOle[1];
					}
					else
					{
						pos = pos2;
						atPos = posOle[1];
						atOther = posOle[0];
					}

				}
				else if (pos1 <= mainChain.size() && pos2 > mainChain.size())
				{
					pos = pos2;
					atPos = posOle[1];
					atOther = posOle[0];
				}
				else if (pos1 > mainChain.size() && pos2 <= mainChain.size())
				{
					pos = pos1;
					atPos = posOle[0];
					atOther = posOle[1];
				}
				else if (pos1 > mainChain.size() && pos2 > mainChain.size())
				{
					//pos = std::max(pos1, pos2);
					if (pos1 > pos2)
					{
						pos = pos1;
						atPos = posOle[0];
						atOther = posOle[1];
					}
					else
					{
						pos = pos2;
						atPos = posOle[1];
						atOther = posOle[0];
					}
				}
				
				if (posKeto > pos)
				{
					Molecola ref = mol.parentFuel();
					ref.addole(atOther, mol.posKeto());
					ref.addcheto(atPos);
					if (ref == mol)
					{
						int temp = posKeto;
						posKeto = pos;
						pos = temp;
					}
				}
				name.append("D");
				name.append(std::to_string(pos));
				name.append("-");
				name.append(std::to_string(posKeto));
				name.append("O");
				if ((mol.posKeto() == posOle[0] || mol.posKeto() == posOle[1]))
					name.append("a");

				//if (isI3C6)
				//{
				//	name.append("-CO");
				//}

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
				//name.append(std::to_string(Add(mol, &linEtheROIsomersList) + 1));
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
				name.append(std::to_string(map[mol.posCOrad() - 1]));
				name.append("O");
				if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
					name.append("a");

				break;
			}
			case oleROO_:
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

				name.append("H");
				name.append(std::to_string(mol.numberOfH()));
				name.append(std::to_string(pos));
				name.append("-");
				name.append(std::to_string(map[mol.posCOOrad() - 1]));
				name.append("O2");
				if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
					name.append("a");

				break;
			}

			case oleQOOH_:
			{
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

				int posCOOH = map[mol.posCOOH() - 1];
				int posCrad = map[mol.posCrad() - 1];
				if (posCrad > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.tipo(mol.posCOOH(), 2);
					ref.tipo(mol.posCrad(), 4);
					if (ref == mol)
					{
						int temp = posCOOH;
						posCOOH = posCrad;
						posCrad = temp;
					}
				}
				name.append("H");
				name.append(std::to_string(mol.numberOfH()-1));
				name.append(std::to_string(pos));
				name.append("OOH");
				name.append(std::to_string(posCOOH));
				name.append("-");
				name.append(std::to_string(posCrad));
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");

				break;
			}
			case QOHOOH_:
			{
				int posCOOH = map[mol.posCOOH() - 1];
				int posCrad = map[mol.posCrad() - 1];

				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 2));
				name.append("OH");
				name.append(std::to_string(map[mol.posCOH() - 1]));

				name.append("Q");
				name.append(std::to_string(posCOOH));
				name.append("-");
				name.append(std::to_string(posCrad));
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}
			case O2QOHOOH_:
			{
				int posCOOH = map[mol.posCOOH() - 1];
				int posCOO = map[mol.posCOOrad() - 1];
				//int posCrad = map[mol.posCrad() - 1];

				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 2));
				name.append("OH");
				name.append(std::to_string(map[mol.posCOH() - 1]));
				name.append("Q");
				name.append(std::to_string(posCOOH));
				name.append("-");
				name.append(std::to_string(posCOO));
				name.append("O2");
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}

			case OHKHP_:
			{
				int posCOOH = map[mol.posCOOH() - 1];
				int posKeto = map[mol.posKeto() - 1];
				if (posKeto > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.addcheto(mol.posCOOH());
					ref.tipo(mol.posKeto(), 4);
					if (ref == mol)
					{
						int temp = posCOOH;
						posCOOH = posKeto;
						posKeto = temp;
					}
				}
				//int posCrad = map[mol.posCrad() - 1];
				name.append("H");
				name.append(std::to_string(mol.numberOfH()-2));
				name.append("OH");
				name.append(std::to_string(map[mol.posCOH() - 1]));
				name.append("KET");
				name.append(std::to_string(posKeto));
				name.append("-");
				name.append(std::to_string(posCOOH));

				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}

			case COHCHO_:
			{
				int posKeto = map[mol.posKeto() - 1];

				//int posCrad = map[mol.posCrad() - 1];
				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 1));
				name.append("OH");
				name.append(std::to_string(map[mol.posCOH() - 1]));
				name.append("KET");
				name.append(std::to_string(posKeto));

				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}


			case POHOOH_:
			{
				std::vector<int> posAtsOOH = mol.posOOHinPOOH2();
				int posAtOOH1 = posAtsOOH[0];
				int posAtOOH2 = posAtsOOH[1];
				int posOOH1 = map[posAtOOH1 - 1];
				int posOOH2 = map[posAtOOH2 - 1];
				int posC = map[mol.posCrad() - 1];
				if (posC > posOOH1)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(posAtOOH2, mol.posCrad());
					ref.tipo(posAtOOH1, 2);
					if (ref == mol)
					{
						int temp = posOOH1;
						posOOH1 = posC;
						posC = temp;
					}
				}
				if (posC > posOOH2)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(posAtOOH1, mol.posCrad());
					ref.tipo(posAtOOH2, 2);
					if (ref == mol)
					{
						int temp = posOOH2;
						posOOH2 = posC;
						posC = temp;
					}
				}

				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 3));
				name.append("OH");
				name.append(std::to_string(map[mol.posCOH() - 1]));
				name.append("Q");
				name.append(std::to_string(std::min(posOOH1, posOOH2)));
				name.append(std::to_string(std::max(posOOH1, posOOH2)));
				name.append("-");
				name.append(std::to_string(posC));
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				
				
				break;




			}
			case OHcEthOOH_:	// molecule is an ether-OOH
			{
				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 2));
				name.append("OH");
				name.append(std::to_string(map[mol.posCOH() - 1]));
				name.append("O");
				std::vector<int> pos = mol.posEthero();
				//int pos1 = std::min(map[pos[0] - 1], map[pos[1] - 1]);
				//int pos2 = std::max(map[pos[0] - 1], map[pos[1] - 1]);
				int posO1;
				int posO2;
				int numAtO1;
				int numAtO2;
				if (map[pos[0] - 1] > map[pos[1] - 1])
				{
					numAtO1 = pos[1];
					posO1 = map[pos[1] - 1];
					numAtO2 = pos[0];
					posO2 = map[pos[0] - 1];
				}
				else
				{
					numAtO1 = pos[0];
					posO1 = map[pos[0] - 1];
					numAtO2 = pos[1];
					posO2 = map[pos[1] - 1];
				}
				int posCOOH = map[mol.posCOOH() - 1];
				if (posO1 > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(numAtO2, mol.posCOOH());
					ref.tipo(numAtO1, 4);
					if (ref == mol)
					{
						int temp = posO1;
						posO1 = posCOOH;
						posCOOH = temp;
					}
				}
				if (posO2 > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(numAtO1, mol.posCOOH());
					ref.tipo(numAtO2, 4);
					if (ref == mol)
					{
						int temp = posO2;
						posO2 = posCOOH;
						posCOOH = temp;
					}
				}
				name.append(std::to_string(std::min(posO1, posO2)));
				//name.append("-");
				name.append(std::to_string(std::max(posO1, posO2)));
				name.append("Q");
				name.append(std::to_string(posCOOH));
				break;
			}

			case OHoleOOH_:
			{
				std::vector<int> posOle = mol.posOle(1);
				int pos1 = map[posOle[0] - 1];
				int pos2 = map[posOle[1] - 1];
				int posOH = mol.posCOH();
				int posCOOH = mol.posCOOH();
				float halfSize = (mainChain.size()+1) / 2.0;
				int posOLE = 0;
				if (pos1 > pos2)	// just to be sure thatn pos1 is smaller than pos2
				{
					int temp = pos1;
					pos1 = pos2;
					pos2 = temp;
				}
				// check if there is priority on olefins numbering
				if ((pos1 < halfSize && pos2 < halfSize) || pos2 == halfSize)
				{
					// keep same numbering
					posOLE = pos1;
				}
				else if ((pos1 > halfSize && pos2 > halfSize) || pos1 == halfSize)
				{
					//flip numbering
					posOLE = mainChain.size() - pos2 + 1;
					posOH = mainChain.size() - posOH + 1;
					posCOOH = mainChain.size() - posCOOH + 1;

				}
				else if (pos1 < halfSize && pos2 > halfSize)
				{
					//give priority to OH
					if (posOH > halfSize)
					{
						// FLIP
						posOLE = mainChain.size() - pos2 + 1;
						posOH = mainChain.size() - posOH + 1;
						posCOOH = mainChain.size() - posCOOH + 1;
					}
					else if (posOH < halfSize)
					{
						// keep same numbering
						posOLE = pos1;
					}
				}

				// std::vector<int> posOle = mol.posOle(1);
				// int pos1 = map[posOle[0] - 1];
				// int pos2 = map[posOle[1] - 1];
				// int pos = 0;
				// if (pos1 <= mainChain.size() && pos2 <= mainChain.size())
				// 	pos = std::min(pos1, pos2);
				// else if (pos1 <= mainChain.size() && pos2 > mainChain.size())
				// 	pos = pos2;
				// else if (pos1 > mainChain.size() && pos2 <= mainChain.size())
				// 	pos = pos1;
				// else if (pos1 > mainChain.size() && pos2 > mainChain.size())
				// 	pos = std::max(pos1, pos2);
				// 
				// int posCOOH = map[mol.posCOOH() - 1];
				// int posCOO = map[mol.posCOOrad() - 1];
				// //int posCrad = map[mol.posCrad() - 1];
				// 
				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 2));
				//name.append("D");
				name.append(std::to_string(posOLE));
				name.append("OH");
				//name.append(std::to_string(map[mol.posCOH() - 1]));
				name.append(std::to_string(posOH));
				name.append("Q");
				name.append(std::to_string(posCOOH));
				// //if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				// //	name.append("a");
				break;
			}
			
			case OHoleketo_:
			{
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
				
				int posOH = map[mol.posCOH() - 1];
				int posKeto = map[mol.posKeto() - 1];

				int len = mol.size();
				int posrev = len - pos;
				int posOHrev = len + 1 - posOH;
				int posKetorev = len + 1 - posKeto; // p4 < p3;


				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 1));
				//name.append("D");
				if (pos <= posrev)
				{
					name.append(std::to_string(pos));
					name.append("OH");
					name.append(std::to_string(posOH));
					name.append("KET");
					name.append(std::to_string(posKeto));
				}
				else
				{
					name.append(std::to_string(posrev));
					name.append("OH");
					name.append(std::to_string(posOHrev));
					name.append("KET");
					name.append(std::to_string(posKetorev));
				}
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}

			case diketo_:
			{
				int posKeto = map[mol.posKeto() - 1];
				int posKeto1= map[mol.posKeto1() - 1];
				name.append("H");
				name.append(std::to_string(mol.numberOfH()));
				//name.append("D");
				name.append("KETO");
				name.append(std::to_string(posKeto));
				name.append(std::to_string(posKeto1));

				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}
			case OHdiketo_:
			{
				int posOH = map[mol.posCOH() - 1];
				int posKeto = map[mol.posKeto() - 1];
				int posKeto1 = map[mol.posKeto1() - 1];
				name.append("H");
				name.append(std::to_string(mol.numberOfH()-1));
				name.append("OH");
				name.append(std::to_string(posOH));
				//name.append("D");
				name.append("KETO");
				name.append(std::to_string(posKeto));
				name.append(std::to_string(posKeto1));

				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}

			case OHoleRketo_:
			{
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

				int posOH = map[mol.posCOH() - 1];
				int posKeto = map[mol.posKeto() - 1];
				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 1));
				name.append("D");
				name.append(std::to_string(pos));
				name.append("OH");
				name.append(std::to_string(posOH));
				//name.append("D");
				name.append("KETO");
				name.append(std::to_string(posKeto));
				name.append("-");
				name.append(std::to_string(map[mol.posCrad() - 1]));

				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}

			case OHRO_:
			{
				int posR = map[mol.posCrad() - 1];
				int posKeto = map[mol.posKeto() - 1];

				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 1));
				name.append("OH");
				name.append(std::to_string(map[mol.posCOH() - 1]));
				name.append("KET");
				name.append(std::to_string(posKeto));
				name.append("R");
				name.append(std::to_string(posR));
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}

			case oleROH_:
			{
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

				int posOH = map[mol.posCOH() - 1];
				int posR = map[mol.posCrad() - 1];

				int len = mol.size();
				int posrev = len - pos;
				int posOHrev = len + 1 - posOH; // p4 < p3;
				int posRrev = len + 1 - posR;

				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 1));
				//name.append("D");
				if (pos <= posrev)
				{
					name.append(std::to_string(pos));
					name.append("OH");
					name.append(std::to_string(posOH));
					name.append("-");
					name.append(std::to_string(posR));
				}
				else
				{
					name.append(std::to_string(posrev));
					name.append("OH");
					name.append(std::to_string(posOHrev));
					name.append("-");
					name.append(std::to_string(posRrev));
				}
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}

			case oleRO_:
			{
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

				int posCrad = map[mol.posCrad() - 1];

				int posKeto = map[mol.posKeto() - 1];

				name.append("H");
				name.append(std::to_string(mol.numberOfH()));
				//name.append("D");
				name.append(std::to_string(pos));
				name.append("KET");
				name.append(std::to_string(posKeto));
				name.append("R");
				name.append(std::to_string(posCrad));
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}

			case ROHOO_:
			{
				int posCOO = map[mol.posCOOrad() - 1];
				//int posCrad = map[mol.posCrad() - 1];

				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 1));
				name.append("OH-");
				name.append(std::to_string(map[mol.posCOH() - 1]));
				name.append("O2-");
				name.append(std::to_string(posCOO));
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}
			case olecEth_:
			{
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

				std::vector<int> pose = mol.posEthero();
				int pose1 = std::min(map[pose[0] - 1], map[pose[1] - 1]);
				int pose2 = std::max(map[pose[0] - 1], map[pose[1] - 1]);


				name.append("H");
				name.append(std::to_string(mol.numberOfH()));
				name.append(std::to_string(pos));
				name.append("O");
				name.append(std::to_string(pose1));
				name.append("-");
				name.append(std::to_string(pose2));
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");

				break;
			}
			case oleO2QOOH_:
			{
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

				int posCOOH = map[mol.posCOOH() - 1];
				int posCOOrad = map[mol.posCOOrad() - 1];
				if (posCOOH > posCOOrad)
				{
					Molecola ref = mol.parentFuel();
					ref.tipo(mol.posCOOH(), 3);
					ref.tipo(mol.posCOOrad(), 4);
					if (ref == mol)
					{
						int temp = posCOOH;
						posCOOH = posCOOrad;
						posCOOrad = temp;
					}
				}

				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 1));
				name.append(std::to_string(pos));
				name.append("OOH");
				name.append(std::to_string(posCOOH));
				name.append("-");
				name.append(std::to_string(posCOOrad));
				name.append("O2");
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");

				break;
			}
			case olePOOH2_:
			{
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


				std::vector<int> posAtsOOH = mol.posOOHinPOOH2();
				int posAtOOH1 = posAtsOOH[0];
				int posAtOOH2 = posAtsOOH[1];
				int posOOH1 = map[posAtOOH1 - 1];
				int posOOH2 = map[posAtOOH2 - 1];
				int posC = map[mol.posCrad() - 1];
				if (posC > posOOH1)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(posAtOOH2, mol.posCrad());
					ref.tipo(posAtOOH1, 2);
					if (ref == mol)
					{
						int temp = posOOH1;
						posOOH1 = posC;
						posC = temp;
					}
				}
				if (posC > posOOH2)
				{
					Molecola ref = mol.parentFuel();
					ref.addetero(posAtOOH1, mol.posCrad());
					ref.tipo(posAtOOH2, 2);
					if (ref == mol)
					{
						int temp = posOOH2;
						posOOH2 = posC;
						posC = temp;
					}
				}

				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 2));
				name.append(std::to_string(pos));
				name.append("Q");
				name.append(std::to_string(std::min(posOOH1, posOOH2)));
				name.append(std::to_string(std::max(posOOH1, posOOH2)));
				name.append("-");
				name.append(std::to_string(posC));
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");

				break;
			}

			case dienesR_:
			{
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

				std::vector<int> posOle2 = mol.posOle(2);
				int pos21 = map[posOle2[0] - 1];
				int pos22 = map[posOle2[1] - 1];

				int pos_1 = 0;
				if (pos21 <= mainChain.size() && pos22 <= mainChain.size())
					pos_1 = std::min(pos21, pos22);
				else if (pos21 <= mainChain.size() && pos22 > mainChain.size())
					pos_1 = pos22;
				else if (pos21 > mainChain.size() && pos22 <= mainChain.size())
					pos_1 = pos21;
				else if (pos21 > mainChain.size() && pos22 > mainChain.size())
					pos_1 = std::max(pos21, pos22);

				int p1 = std::min(pos, pos_1);
				int p2 = std::max(pos, pos_1);
				int posR = map[mol.posCrad() - 1];

				name.append("H");
				name.append(std::to_string(mol.numberOfH()));
				name.append("D");
				name.append(std::to_string(p1));
				name.append(std::to_string(p2));
				name.append("-");
				name.append(std::to_string(posR));

				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}
			case dienes_:
			{
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

				std::vector<int> posOle2 = mol.posOle(2);
				int pos21 = map[posOle2[0] - 1];
				int pos22 = map[posOle2[1] - 1];

				int pos_1 = 0;
				if (pos21 <= mainChain.size() && pos22 <= mainChain.size())
					pos_1 = std::min(pos21, pos22);
				else if (pos21 <= mainChain.size() && pos22 > mainChain.size())
					pos_1 = pos22;
				else if (pos21 > mainChain.size() && pos22 <= mainChain.size())
					pos_1 = pos21;
				else if (pos21 > mainChain.size() && pos22 > mainChain.size())
					pos_1 = std::max(pos21, pos22);

				int p1 = std::min(pos, pos_1);
				int p2 = std::max(pos, pos_1);
				int len = mol.size();
				int p3 = len - p1;
				int p4 = len - p2; // p4 < p3;


				name.append("H");
				name.append(std::to_string(mol.numberOfH()));
				name.append("D");
				if (p1 <= p4)
				{
					name.append(std::to_string(p1));
					name.append(std::to_string(p2));
				}
				else
				{
					name.append(std::to_string(p4));
					name.append(std::to_string(p3));
				}

				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}
			
			case dienesCO_:
			{
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

				std::vector<int> posOle2 = mol.posOle(2);
				int pos21 = map[posOle2[0] - 1];
				int pos22 = map[posOle2[1] - 1];

				int pos_1 = 0;
				if (pos21 <= mainChain.size() && pos22 <= mainChain.size())
					pos_1 = std::min(pos21, pos22);
				else if (pos21 <= mainChain.size() && pos22 > mainChain.size())
					pos_1 = pos22;
				else if (pos21 > mainChain.size() && pos22 <= mainChain.size())
					pos_1 = pos21;
				else if (pos21 > mainChain.size() && pos22 > mainChain.size())
					pos_1 = std::max(pos21, pos22);

				int p1 = std::min(pos, pos_1);
				int p2 = std::max(pos, pos_1);
				int len = mol.size();
				int p3 = len - p1;
				int p4 = len - p2; // p4 < p3;

				int posKeto = map[mol.posKeto() - 1];
				int posKeto1 = len + 1 - posKeto;
				int poskt = std::min(posKeto, posKeto1);

				if (p1 <= p4)
				{
					if (p1 == p4 && p2 == p3)
					{
						name.append("H");
						name.append(std::to_string(mol.numberOfH()));
						name.append("D");
						name.append(std::to_string(p1));
						name.append(std::to_string(p2));
						name.append("KET");
						name.append(std::to_string(poskt));
					}
					else {
						name.append("H");
						name.append(std::to_string(mol.numberOfH()));
						name.append("D");
						name.append(std::to_string(p1));
						name.append(std::to_string(p2));
						name.append("KET");
						name.append(std::to_string(posKeto));
					}
				}
				else 
				{
					name.append("H");
					name.append(std::to_string(mol.numberOfH()));
					name.append("D");
					name.append(std::to_string(p4));
					name.append(std::to_string(p3));
					name.append("KET");
					name.append(std::to_string(posKeto1));
				}
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");

				break;
			}

			case dienesOH_:
			{
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

				std::vector<int> posOle2 = mol.posOle(2);
				int pos21 = map[posOle2[0] - 1];
				int pos22 = map[posOle2[1] - 1];

				int pos_1 = 0;
				if (pos21 <= mainChain.size() && pos22 <= mainChain.size())
					pos_1 = std::min(pos21, pos22);
				else if (pos21 <= mainChain.size() && pos22 > mainChain.size())
					pos_1 = pos22;
				else if (pos21 > mainChain.size() && pos22 <= mainChain.size())
					pos_1 = pos21;
				else if (pos21 > mainChain.size() && pos22 > mainChain.size())
					pos_1 = std::max(pos21, pos22);

				int p1 = std::min(pos, pos_1);
				int p2 = std::max(pos, pos_1);


				name.append("H");
				name.append(std::to_string(mol.numberOfH()-1));
				name.append("D");
				name.append(std::to_string(p1));
				name.append(std::to_string(p2));
				name.append("OH");
				name.append(std::to_string(map[mol.posCOH() - 1]));
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}

			case dieneOOH_:
			{
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

				std::vector<int> posOle2 = mol.posOle(2);
				int pos21 = map[posOle2[0] - 1];
				int pos22 = map[posOle2[1] - 1];

				int pos_1 = 0;
				if (pos21 <= mainChain.size() && pos22 <= mainChain.size())
					pos_1 = std::min(pos21, pos22);
				else if (pos21 <= mainChain.size() && pos22 > mainChain.size())
					pos_1 = pos22;
				else if (pos21 > mainChain.size() && pos22 <= mainChain.size())
					pos_1 = pos21;
				else if (pos21 > mainChain.size() && pos22 > mainChain.size())
					pos_1 = std::max(pos21, pos22);

				int p1 = std::min(pos, pos_1);
				int p2 = std::max(pos, pos_1);
				int len = mol.size();
				int p3 = len - p1;
				int p4 = len - p2; // p4 < p3;

				int posCOOH = map[mol.posCOOH() - 1];
				int posCOOH1 = len + 1 - posCOOH;




				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 1));
				name.append("D");
				if (p1 <= p4)
				{
					name.append(std::to_string(p1));
					name.append(std::to_string(p2));
					name.append("-");
					name.append(std::to_string(posCOOH));
					name.append("OOH");
				}
				else 
				{
					name.append(std::to_string(p4));
					name.append(std::to_string(p3));
					name.append("-");
					name.append(std::to_string(posCOOH1));
					name.append("OOH");
				}

				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");
				break;
			}
			case ROH_: // example: C7H14OH-1J2
			{
				name.append("H");
				name.append(std::to_string(mol.numberOfH() - 1));
				//std::cout << "number of Hs" << mol.numberOfH() << std::endl;
				name.append("OH-");
				name.append(std::to_string(map[mol.posCOH() - 1]));
				//std::cout << "OH pos " << std::to_string(map[mol.posCOH() - 1]) << std::endl;
				name.append("J");
				name.append(std::to_string(map[mol.posCrad() - 1]));
				//std::cout << "Radical pos " << std::to_string(map[mol.posCrad() - 1]) << std::endl;

				break;
			}
			case oleKHP_:
			{
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

				int posCOOH = map[mol.posCOOH() - 1];
				int posKeto = map[mol.posKeto() - 1];
				if (posKeto > posCOOH)
				{
					Molecola ref = mol.parentFuel();
					ref.addcheto(mol.posCOOH());
					ref.tipo(mol.posKeto(), 4);
					if (ref == mol)
					{
						int temp = posCOOH;
						posCOOH = posKeto;
						posKeto = temp;
					}
				}
				int len = mol.size();
				int posrev = len - pos;
				int posKetorev = len + 1 - posKeto; // p4 < p3;
				int posCOOHrev = len + 1 - posCOOH;


				name.append("H");
				name.append(std::to_string(mol.numberOfH()));
				if (pos <= posrev)
				{
					name.append(std::to_string(pos));
					name.append("KET");
					name.append(std::to_string(posKeto));
					//name.append("-");
					name.append(std::to_string(posCOOH));
				}
				else 
				{
					name.append(std::to_string(posrev));
					name.append("KET");
					name.append(std::to_string(posKetorev));
					//name.append("-");
					name.append(std::to_string(posCOOHrev));
				}
				//if ((mol.posCOrad() == posOle[0] || mol.posCOrad() == posOle[1]))
				//	name.append("a");

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

	if (name.size() > 16)
	{
		bool found = false;
		for (int i = 0; i < longNamedSpecies.size(); i++)
		{
			if (mol == longNamedSpecies[i])
			{
				name = longNamedSpeciesNames[i];
				found = true;
			}
		}
		if (found == false)
		{
			longNamedSpeciesOriginalNames.push_back(name);
			name = "MOL" + std::to_string(longNamedSpecies.size()+1);
			longNamedSpecies.push_back(mol);
			longNamedSpeciesNames.push_back(name);
		}
	}
	
	return name;
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

	std::string detName = molToName(mol);
	if (specMol == 0 && UTL::isPresent(&(baseMech)->mechSpecies, detName))
	{
		return detName;
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
		case 11:
			name.append("HCCO");
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

std::string ChemkinOut::checkIfNamed(Molecola* mol)
{
	std::string name = "not found";
	//std::cout << *mol << std::endl;

	std::string molInchi = mol->inchiName();
	//if (molInchi == "InChI=1S/C2H4/c1-2/h1-2H2")
	//	std::cout << "here" << std::endl;
	if (mol->numKYN() == 1)
	{
		molInchi = "InChI=1S/C2H2/c1-2/h1-2H"; // C2H2
		if (mol->numCOH() == 1)	molInchi = "InChI=1S/C2H2O/c1-2-3/h1,3H"; //HCCOH
	}
		

	

	for(int i = 0; i < (*baseMech).glossaryInChIs.size(); i++)
	{

		if (molInchi == (*baseMech).glossaryInChIs[i])
		{
			return (*baseMech).glossarySpecies[i];
		}
	}
	return name;
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
	// JIAXIN
	//if (mol->isLinear())
	//	prefix.append("N");
	//else
	if (mol->isLinear() == false)
	{
		bool alreadyNamed = false;
		if (numMe == 2 && numEt == 0)
		{
			if ((mePos[0] == 2 && mePos[1] == 2))
			{
				prefix.append("NE");
				alreadyNamed = true;
			}
			else if ((mePos[0] == 2 && mePos[1] == mainChain.size()-1)
				|| (mePos[0] == mainChain.size() - 1 && mePos[1] == 2))
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
		prefix = "IS2";
	}
	if ((mol->numberOfC() == 15) && (prefix == "I224668"))
	{
		prefix = "IS3";
	}
	if ((mol->numberOfC() == 15) && (prefix == "I224688"))
	{
		prefix = "IS4";
	}
	if ((mol->numberOfC() == 12) && (prefix == "I22466"))
	{
		prefix = "ISO";
	}
	if ((mol->numberOfC() == 12) && (prefix == "I22466"))
	{
		prefix = "IS2";
	}
	if ((mol->numberOfC() == 12) && (prefix == "I22446"))
	{
		prefix = "IS3";
	}
	if ((mol->numberOfC() == 14) && (prefix == "I22468"))
	{
		prefix = "IS1";
	}
	if ((mol->numberOfC() == 14) && (prefix == "I22488"))
	{
		prefix = "IS2";
	}
	if ((mol->numberOfC() == 14) && (prefix == "I22448"))
	{
		prefix = "IS3";
	}
	if ((mol->numberOfC() == 14) && (prefix == "I22668"))
	{
		prefix = "IS4";
	}
	if ((mol->numberOfC() == 14) && (prefix == "I24468"))
	{
		prefix = "IS5";
	}
	if ((mol->numberOfC() == 11) && (prefix == "I2244"))
	{
		prefix = "IS1";
	}
	if ((mol->numberOfC() == 11) && (prefix == "I2246"))
	{
		prefix = "IS2";
	}
	if ((mol->numberOfC() == 11) && (prefix == "I2446"))
	{
		prefix = "IS3";
	}
	if ((mol->numberOfC() == 10) && (prefix == "I2244"))
	{
		prefix = "IS1";
	}

	return prefix;
}

std::string ChemkinOut::specialMolName(Molecola mol)
{
	int specMolId = mol.isSpecialMolecule();
	switch (specMolId)
	{
	case 1:
		return "N2";
		break;
	case 2:
		return "O2";
		break;
	case 3:
		return "O";
		break;
	case 4:
		return "OH";
		break;
	case 5:
		return "HO2";
		break;
	case 6:
		return "H2O";
		break;
	case 7:
		return "H2O2";
		break;
	case 8:
		return "H";
		break;
	case 9:
		return "H2";
		break;
	case 10:
		return "CO";
		break;
	default:
		UTL::error("isSpecialMolecule called on a non special molecule");
		return "ERROR";
	}

	return "ERROR";

}

void ChemkinOut::printLongNameSpeciesMessage()
{
	for (int i = 0; i < longNamedSpecies.size(); i++)
	{
		std::stringstream msg;
		msg << "Species with long name found:" << std::endl;
		msg << "   Species " << longNamedSpecies[i] << " has the following automatic name "
			<< longNamedSpeciesOriginalNames[i] << " that exceeds the 16 characters limit." << std::endl;
		msg << "   The name " << longNamedSpeciesNames[i] << " has been given to it." << std::endl;
		msg << "   The preferred name can be defined in the glossary for the following InChI identifier:" << std::endl;
		msg << "       " << longNamedSpecies[i].inchiName() << std::endl;
		UTL::warning(msg.str());
	}
}
