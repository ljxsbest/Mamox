
#include "molec.h"
#include <string.h>

void Molecola::init(void)
{
	conc_v.resize(numOfTemps, 0.0);
	conc = 0;
	velo = 0;
	velo_v.resize(numOfTemps, 0.0);
	selectivities.resize(numOfTemps, 0.0);
	//velo1 = 0;
	//velo2 = 0;
	isomeri = 1;
	output = FORMULA;
	numole = numete = numO = 0;
	
	for (int r = 0; r <= SIZEMAX; r++)
	{
		for (int c = 0; c <= 5; c++) struttura[r][c] = 0;
		ole[r][0] = ole[r][1] = 0;
		ete[r][0] = ete[r][1] = 0;
		doppioO[r] = 0;
	};
	specialMolecule = 0;
};

Molecola::Molecola(void)
{
	init();
};

Molecola::Molecola(int specialMoleculeID)
{
	init();
	specialMolecule = specialMoleculeID;
}


Molecola::Molecola(char nomefile[])
{
	init();
	ifstream myfile(nomefile);
	if (!myfile)
	{
		cerr << "Errore nell' apertura del file: " << nomefile;
		exit(1);
	};
	myfile >> *this;
	myfile.close();
};

Molecola::Molecola(int size_, std::vector<int> structure, int numOle_, std::vector<int> oleBonds, int numEthe_,
	std::vector<int> etheBonds, int numKeto_, std::vector<int> ketoPos)
{
	init();

	for (int i = 0; i < size_; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			struttura[i + 1][j + 1] = structure[i * 5 + j];
		}
	}

	//numole, numete, numO

	numole = numOle_;
	for (int i = 0; i < numOle_; i++)
	{
		ole[i + 1][0] = oleBonds[i * 2];
		ole[i + 1][1] = oleBonds[i * 2 + 1];
	}

	numete = numEthe_;
	for (int i = 0; i < numEthe_; i++)
	{
		ete[i + 1][0] = etheBonds[i * 2];
		ete[i + 1][1] = etheBonds[i * 2 + 1];
	}

	numO = numKeto_;
	for (int i = 0; i < numKeto_; i++)
	{
		doppioO[i + 1] = ketoPos[i];
	}
}

int Molecola::isole(int pos1, int pos2)
{
	for (int i = 1; i <= numole; i++)
	{
		if (ole[i][0] == pos1 && ole[i][1] == pos2) return 1;
		if (ole[i][1] == pos1 && ole[i][0] == pos2) return 1;
	};
	return 0;
};

int Molecola::isole(int pos)
{
	for (int i = 1; i <= numole; i++)
	{
		if (ole[i][0] == pos || ole[i][1] == pos) return 1;
	};
	return 0;
};

int Molecola::isole(void)
{
	if (numole > 0) return 1;
	return 0;
};

int Molecola::isCRad(void)		// true if molecule is a C* radical, false otherwise.
{
	int i;
	for (i = 1; i <= SIZEMAX; i++)
		if (tipo(i) == 2)  return 1;
	return 0;
}

int Molecola::addole(int pos1, int pos2)
{
	if (dist(pos1, pos2) != 1) return 0;
	if (isole(pos1, pos2)) return 0;
	//if (numAbstractableH(pos1) < 1 || numAbstractableH(pos2) < 1)  // non ci puo' stare il legame
	//{
	//	if (tipo(pos1) == 2 && tipo(pos2) == 2)   // a meno che non ci siano due 
	//	{                                  // radicali vicini
	//		tipo(pos1, 1);
	//		tipo(pos2, 1);
	//	}
	//	else { return 0; };
	//};
	if (numAbstractableH(pos1) < 1)
	{
		if (tipo(pos1) == 2)
			tipo(pos1, 1);
		else
		{
			std::cerr << "WARNING: addole called on a carbon without hydrogens" << std::endl;
			return 0;
		}
	}
	if (numAbstractableH(pos2) < 1)
	{
		if (tipo(pos2) == 2)
			tipo(pos2, 1);
		else
		{
			std::cerr << "WARNING: addole called on a carbon without hydrogens" << std::endl;
			return 0;
		}
	}


	// se ho passato tutti i controlli...
	numole++;
	ole[numole][0] = pos1;
	ole[numole][1] = pos2;
	return 1;
};

int Molecola::isetero(int pos1, int pos2)
{
	for (int i = 1; i <= numete; i++)
	{
		if (ete[i][0] == pos1 && ete[i][1] == pos2) return 1;
		if (ete[i][1] == pos1 && ete[i][0] == pos2) return 1;
	};
	return 0;
};

int Molecola::isetero(int pos)
{
	for (int i = 1; i <= numete; i++)
		if (ete[i][0] == pos || ete[i][1] == pos) return 1;
	return 0;
};

int Molecola::isetero(void)
{
	if (numete > 0) return 1;
	return 0;
};

int Molecola::addetero(int pos1, int pos2)
{
	if (tipoC(pos1) == Cq) return 0;
	if (tipoC(pos2) == Cq) return 0;
	numete++;
	ete[numete][0] = pos1;
	ete[numete][1] = pos2;
	return 1;
};

int Molecola::removeOO(int pos)
{
	if (struttura[pos][1] != 3)
	{
		std::cerr << "ERROR: removeOO called on a non COO carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 1;
	return 1;
}

int Molecola::removeOOH(int pos)
{
	if (struttura[pos][1] != 4)
	{
		std::cerr << "ERROR: removeOOH called on a non COOH carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 1;
	return 1;
}

int Molecola::removeOle(int num)		// remove the num-th olefin (order in which they are saved in the ole array), if possible return 1, if not possible return 0
{
	if (numole == 0)
	{
		std::cerr << "WARNING: removeOle called on a non olefins!" << std::endl;
		return 0;
	}
	if (num > numole)
	{
		std::cerr << "WARNING: removeOle called on olefin number "<< num << " but there are only " << numole << "olefins!" << std::endl;
		return 0;
	}
	if (num < 1)
	{
		std::cerr << "ERROR: removeOle called with position " << num << "! aborted." << std::endl;
		exit(0);
		return 0;
	}

	// remove the olefin
	for (int i = num; i < numole;)
	{
		ole[num][0] = ole[num + 1][0];
		ole[num][1] = ole[num + 1][1];
	}
	ole[numole][1] = 0;
	ole[numole][0] = 0;

	numole--;
	return 1;
}

int Molecola::removeOle(int pos1, int pos2)		// remove the double bond between the two carbons in pos1 and pos2, if possible return 1, if not possible return 0
{
	int olePos = 0;
	if (pos1 < 1 || pos2 < 1)
	{
		std::cerr << "ERROR: removeOle called in positions " << pos1 << " and " << pos2 << " ! aborted." << std::endl;
		exit(0);
		return 0;
	}
	if (numole == 0)
	{
		std::cerr << "WARNING: removeOle called on a non olefins!" << std::endl;
		return 0;
	}
	for (int i = 1; i <= numole; i++)
	{
		if ((pos1 == ole[i][0] && pos2 == ole[i][1]) || (pos1 == ole[i][1] && pos2 == ole[i][0]))
		{
			olePos = i;
			break;
		}
	}
	if (olePos == 0)
	{
		std::cerr << "WARNING: removeOle was not able to find any double bond between the positions " << pos1 << " and " << pos2 << std::endl;
		return 0;
	}
	removeOle(olePos);
	return 1;
}

int Molecola::removeCycEther(int num)		// remove the num-th cyclic ether (order in which they are saved in the ete array), if possible return 1, if not possible return 0
{
	if (numete == 0)
	{
		std::cerr << "WARNING: removeCycEther called on a non cyclyc ether containing molecule!" << std::endl;
		return 0;
	}
	if (num > numete)
	{
		std::cerr << "WARNING: removeCycEther called on ether number " << num << " but there are only " << numete << "cyclic ethers!" << std::endl;
		return 0;
	}
	if (num < 1)
	{
		std::cerr << "ERROR: removeCycEther called with position " << num << "! aborted." << std::endl;
		exit(0);
		return 0;
	}

	// remove the cyclic ether
	for (int i = num; i < numete;)
	{
		ete[num][0] = ete[num + 1][0];
		ete[num][1] = ete[num + 1][1];
	}
	ete[numete][0] = 0;
	ete[numete][1] = 0;

	numete--;
	return 1;
}

void Molecola::removeAllKeto()
{
	for (int i = 0; i < numO; i++)
	{
		doppioO[i + 1] = 0;
	}
	numO = 0;
}

int Molecola::removeCycEther(int pos1, int pos2)		// remove the cyclic ether between the two carbons in pos1 and pos2, if possible return 1, if not possible return 0
{
	int etePos = 0;
	if (pos1 < 1 || pos2 < 1)
	{
		std::cerr << "ERROR: removeCycEther called in positions " << pos1 << " and " << pos2 << " ! aborted." << std::endl;
		exit(0);
		return 0;
	}
	if (numete == 0)
	{
		std::cerr << "WARNING: removeCycEther called on a non cyclyc ether containing molecule!" << std::endl;
		return 0;
	}
	for (int i = 1; i <= numete; i++)
	{
		if ((pos1 == ete[i][0] && pos2 == ete[i][1]) || (pos1 == ete[i][1] && pos2 == ete[i][0]))
		{
			etePos = i;
			break;
		}
	}
	if (etePos == 0)
	{
		std::cerr << "WARNING: removeCycEther was not able to find any cyclic ether bonded at the positions " << pos1 << " and " << pos2 << std::endl;
		return 0;
	}
	removeCycEther(etePos);
	return 1;
}

int Molecola::OO_to_OOH(int pos)
{
	if (struttura[pos][1] != 3)
	{
		std::cerr << "ERROR: OO_to_OOH called on a non COO carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 4;
	return 1;
}

int Molecola::C_to_Crad(int pos)
{
	if (struttura[pos][1] != 1)
	{
		std::cerr << "ERROR: C_to_Crad called on a non C carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 2;
	return 1;
}

int Molecola::Crad_to_C(int pos)
{
	if (struttura[pos][1] != 2)
	{
		std::cerr << "ERROR: Crad_to_C called on a non C* carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 1;
	return 1;
}

int Molecola::Crad_to_COOH(int pos)
{
	if (struttura[pos][1] != 2)
	{
		std::cerr << "ERROR: Crad_to_COOH called on a non C* carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 4;
	return 1;
}

int Molecola::COOrad_to_COOH(int pos)
{
	if (struttura[pos][1] != 3)
	{
		std::cerr << "ERROR: COOrad_to_COOH called on a non COOrad carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 4;
	return 1;
}

int Molecola::COrad_to_CO(int pos)
{
	if (struttura[pos][1] != 9)
	{
		std::cerr << "ERROR: COrad_to_CO called on a non COrad carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 1;
	addcheto(pos);
	return 1;
}

int Molecola::Crad_to_keto(int pos)
{
	if (struttura[pos][1] != 2)
	{
		std::cerr << "ERROR: Crad_to_keto called on a non C* carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 1;
	addcheto(pos);
	return 1;
}

int Molecola::Crad_to_COrad(int pos)
{
	if (struttura[pos][1] != 2)
	{
		std::cerr << "ERROR: Crad_to_COrad called on a non C* carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 9;
	return 1;
}

int Molecola::C_to_COrad(int pos)
{
	if (struttura[pos][1] != 1)
	{
		std::cerr << "ERROR: C_to_COrad called on a non C carbon!" << std::endl;
		return 0;
	}
	struttura[pos][1] = 9;
	return 1;
}

bool Molecola::isGroup(int pos)
{
	if (tipo(pos) != 1)
		return true;
	if (isole(pos))
		return true;
	if (isetero(pos))
		return true;
	if (ischeto(pos))
		return true;
	return false;
}

bool Molecola::isNearGroup(int pos)
{
	for (int i = 2; i < 6; i++)
		if (struttura[pos][i] != 0)
			if (isGroup(struttura[pos][i]))
				return true;
	return false;
}

int Molecola::removeAtom(int pos)
{
	//std::cout << *this << std::endl;
	if (pos > size())
		return 0;

	int map[SIZEMAX + 1];
	for (int i = 0; i < SIZEMAX + 1; i++)
		map[i] = 0;
	map[pos] = 0;
	for (int i = 0; i < pos; i++)
		map[i] = i;
	for (int i = pos; i < SIZEMAX; i++)	// shift all the atoms
	{
		for (int j = 0; j < 6; j++)
		{
			map[i + 1] = i;
			struttura[i][j] = struttura[i+1][j];
		}
	}

	for (int i = 0; i < SIZEMAX; i++)	//change the numbering of the atoms
	{
		for (int j = 2; j < 6; j++)
			struttura[i][j] = map[struttura[i][j]];
		//std::cout << ole[i][0] << std::endl;
		ole[i][0] = map[ole[i][0]];
		ole[i][1] = map[ole[i][1]];
		ete[i][0] = map[ete[i][0]];
		ete[i][1] = map[ete[i][1]];
		doppioO[i] = map[doppioO[i]];
	}

	for (int i = 0; i < SIZEMAX; i++)	//remove "holes from bond matrix"
	{
		for (int j = 2; j < 5; j++)
		{
			if(struttura[i][j] == 0)
			{
				struttura[i][j] = struttura[i][j+1];
				struttura[i][j + 1] = 0;
			}
		}
	}


	return 1;
}

int Molecola::ischeto(int pos)
{
	for (int i = 1; i <= numO; i++)
		if (doppioO[i] == pos) return 1;
	return 0;
};

int Molecola::ischeto(void)
{
	if (numO > 0) return doppioO[numO];
	return 0;
};

int Molecola::isKetene(int pos)
{
	if (isole(pos) == true && ischeto(pos) == true)
	{
		return 2;
	}
	else if (isole(pos) == true)
	{
		for (int i = 2; i < 6; i++)
			if (isole(struttura[pos][i])==true && ischeto(struttura[pos][i]) == true)
				return 1;
	}
	return 0;
}

bool Molecola::isInEthRing(int pos)
{
	if (isetero(pos))
		return true;
	for (int i = 0; i < numete; i++)
	{
		int pos1 = ete[i + 1][0];
		int pos2 = ete[i + 1][1];
		int distance = dist(pos1, pos2);
		switch (distance)
		{
		case 1:		// three members ring
			// already satisfied by iseter(pos) check
			break;
		case 2 :	// four members ring
		{
			int found1[SIZEMAX + 1];
			int found2[SIZEMAX + 1];
			scorri(pos1, 1, found1);
			scorri(pos2, 1, found2);
			if (found1[pos] == 1 && found2[pos] == 1)
				return true;
			break;
		}
		case 3:		// five members ring
		{
			int found1[SIZEMAX + 1];
			int found2[SIZEMAX + 1];
			scorri(pos1, 1, found1);
			scorri(pos2, 1, found2);
			if ((found1[pos] == 1 && dist(pos, pos2) == 2)
				|| (found2[pos] == 1 && dist(pos, pos1) == 2))
				return true;
			break;
		}
		default:
			std::cerr << "WARNING: in isInEthRing a ring with " << distance+2 << " members has been found." << std::endl;
		}
	}

	return false;
}

int Molecola::isRO()
{
	if (numCrad() == 1 && numKeto() == 1) return 1;
	return 0;
}

bool Molecola::isAldehyde()			// true if the molecules is an aldehyde (no other functional group present)
{
	if (kindOfSPecies() == CO_)
	{
		int posOxy = posKeto();
		if (numAbstractableH(posOxy) == 1)
			return true;
	}
	return false;
}

bool Molecola::isOleAldehyde()		// true if the molecules is an olefin and aldehyde (no other functional group present  and double bond not involving the carbonyl group)
{
	if (kindOfSPecies() == oleCO_)
	{
		int posOxy = posKeto();
		if (numAbstractableH(posOxy) == 1)
			return true;
	}
	return false;
}


int Molecola::addcheto(int pos)
{
	if (tipoC(pos) == Ct) return 0;
	if (tipoC(pos) == Cq) return 0;
	for (int i = 0; i < numO; i++)
	{
		if (doppioO[i + 1] == pos)
		{
			std::cerr << "WARNING: addcheto used on an already cheto atom!" << std::endl;
			return 0;
		}
	}
	numO++;
	doppioO[numO] = pos;
	return 1;
};

int Molecola::size(void)
{
	int j = 0;
	int i;
	for (i = 1; i <= SIZEMAX; i++)
		if (tipo(i) == 1 || tipo(i) == 2 ||
			tipo(i) == 3 || tipo(i) == 4 ||
			tipo(i) == 8 || tipo(i) == 9)  j++;
	return j;
};

Carbonio Molecola::tipoC(int gruppo)
{
	int trovati = 0;
	for (int i = 1; i <= 4; i++)
		if ((tipo(legato(gruppo, i)) < 5 &&
			tipo(legato(gruppo, i)) > 0) || tipo(legato(gruppo, i)) == 9) trovati++;
	switch (trovati)
	{
	case 1: return Cp;
	case 2: return Cs;
	case 3: return Ct;
	case 4: return Cq;
	case 0: return Cp;   // metano
	};
	cerr << "Errore in tipoC() " << endl;
	exit(100);
};

Idrogeno Molecola::tipoH(int pos)
{
	if (tipo(pos) == 3 || tipo(pos) == 4) return Hcooh;
	switch (tipoC(pos))
	{
	case Cp: return Hp;
	case Cs: return Hs;
	case Ct: return Ht;
	case Cq: cerr << "Chiamata a tipoH su C quaternario!!\n\n";
		exit(0);
	};
};

Radicale Molecola::tipoR(int pos)
{
	if (tipo(pos) != 2)
	{
		cerr << "Chiamata a tipoR non su radicale\n\n";
		exit(0);
	};
	switch (tipoC(pos))
	{
	case Cp: if (size() == 1) { return Rpmet; }
		   else if (size() == 2) { return Rpet; }
		   else { return Rp; };
	case Cs: return Rs;
	case Ct: return Rt;
	case Cq: cerr << "Chiamata a tipoR trova radicale quat !!\n\n";
		exit(0);
	};
};

Radicale Molecola::tipoROO(int pos)
{
	if (tipo(pos) != 3)
	{
		cerr << "Chiamata a tipoROO non su radicale\n\n";
		exit(0);
	};
	switch (tipoC(pos))
	{
	case Cp: if (size() == 1) { return Rpmet; }
		   else if (size() == 2) { return Rpet; }
		   else { return Rp; };
	case Cs: return Rs;
	case Ct: return Rt;
	case Cq: cerr << "Chiamata a tipoROO trova radicale quat !!\n\n";
		exit(0);
	};
};

Radicale Molecola::tipoROOH(int pos)
{
	if (tipo(pos) != 4)
	{
		cerr << "Chiamata a tipoROOH non su OOH\n\n";
		exit(0);
	};
	switch (tipoC(pos))
	{
	case Cp: if (size() == 1) { return Rpmet; }
		   else if (size() == 2) { return Rpet; }
		   else { return Rp; };
	case Cs: return Rs;
	case Ct: return Rt;
	case Cq: cerr << "Chiamata a tipoROOH trova OOH quat !!\n\n";
		exit(0);
	};
};

int Molecola::numH(int gruppo)	// TO DO add all groups
{
	int H = 0;
	if (tipo(gruppo) > 4 && tipo(gruppo) != 8 && tipo(gruppo) != 9)
	{
		cerr << "Chiamata a numH su gruppo " << tipo(gruppo) << " \n\n"; exit(0);
	};
	if (tipo(gruppo) == 1) H = (4 - legami(gruppo));
	if (tipo(gruppo) == 2) H = (3 - legami(gruppo));
	if (tipo(gruppo) == 3) H = (3 - legami(gruppo));
	if (tipo(gruppo) == 4) H = (4 - legami(gruppo));
	if (tipo(gruppo) == 8) H = 0;
	if (tipo(gruppo) == 9) 
		H = (3 - legami(gruppo));
	//else H = (3 - legami(gruppo));
	if (ischeto(gruppo)) H -= 2;
	if (isetero(gruppo)) H -= 1;
	H -= numDoubleBonds(gruppo);
	return H;
};

int Molecola::numAbstractableH(int gruppo)	// TO DO add all groups
{
	int H = 0;
	if (tipo(gruppo) > 4 && tipo(gruppo) != 8 && tipo(gruppo) != 9)
	{
		cerr << "Chiamata a numAbstractableH su gruppo " << tipo(gruppo) << " \n\n"; exit(0);
	};
	if (tipo(gruppo) == 1) H = (4 - legami(gruppo));
	if (tipo(gruppo) == 2) H = (3 - legami(gruppo));
	if (tipo(gruppo) == 3) H = (3 - legami(gruppo));
	if (tipo(gruppo) == 4) H = (3 - legami(gruppo));
	if (tipo(gruppo) == 8) H = 0;
	if (tipo(gruppo) == 9) H = (3 - legami(gruppo));
	//else H = (3 - legami(gruppo));
	if (ischeto(gruppo)) H -= 2;
	if (isetero(gruppo)) H -= 1;
	H -= numDoubleBonds(gruppo);
	return H;
};

int Molecola::numGraphConnections(int pos)
{
	int conn = 0;
	for (int i = 2; i < 6; i++)
		if (struttura[pos][i] != 0)
			conn++;
		
	return conn;
}

int Molecola::numBondedCarbons(int pos)
{
	int C = 0;
	for (int i = 2; i < 6; i++)
	{
		int bondedAtom = struttura[pos][i];
		if (bondedAtom == 0)
			continue;
		if (tipo(bondedAtom) != 6 && tipo(bondedAtom) != 7 && tipo(bondedAtom) != 8)
			C++;
	}
	return C;
}

int Molecola::numSingleBondedCarbons(int pos)
{
	int C = 0;
	for (int i = 2; i < 6; i++)
	{
		int bondedAtom = struttura[pos][i];
		if (bondedAtom != 0 && tipo(bondedAtom) != 8 && isole(bondedAtom) == false && ischeto(bondedAtom) == false)
			C++;
	}
	return C;
}

int Molecola::numDoubleBondedCarbons(int pos)
{
	int C = 0;
	for (int i = 2; i < 6; i++)
	{
		int bondedAtom = struttura[pos][i];
		if (bondedAtom != 0 && tipo(bondedAtom) != 8 && isole(bondedAtom) == true)
			C++;
	}
	return C;
}

int Molecola::numBondedKetoCarbons(int pos)
{
	int C = 0;
	for (int i = 2; i < 6; i++)
	{
		int bondedAtom = struttura[pos][i];
		if (ischeto(bondedAtom) == true)
			C++;
	}
	return C;
}

int Molecola::legami(int gruppo)
{
	int numLeg=0;
	//for (i = 1; i <= 4; i++)
	//	if (legato(gruppo, i) == 0) break;
	//return (i - 1);
	for (int i = 0; i < 4; i++)
		if (legato(gruppo, i + 1) != 0)
			numLeg++;
	return numLeg;
};

int Molecola::trova(int gruppo)
{
	if (gruppo < 5 || gruppo == 8 || gruppo == 9) {		// find the groups saved in the struttura metrix (C C* COO* and COOH)
		int pos = 1;
		for (pos = 1; pos <= SIZEMAX; pos++)
			if (tipo(pos) == gruppo) break;
		if (pos <= SIZEMAX) return pos;
		else return 0;
	}
	else if (gruppo == 7)	//find the group =O			// TO DO: add the possibility of finding more than one group
	{
		if (numO == 0) return 0;	// no group =0 present in the molecule
		else return doppioO[1];		// return the first =O group (we are assuming there is only one)
	}
};

int Molecola::spezza(int g1, int g2, Molecola* a, Molecola* b)
{
	if (g1 == g2)
	{
		std::cerr << "WARNING: function spezza was called with two identical positions. This action has been skipped." << std::endl;
		return 0;
	}
	if (areBonded(g1, g2) == 0)
	{
		std::cerr << "WARNING: function spezza was called with two non bonded atoms. This action has been skipped." << std::endl;
		return 0;
	}
	if (isole(g1, g2)) return 0;
	//if (isetero()) return 0;
	if (trova(2) == g1 || trova(2) == g2) return 0;
	int fatti[SIZEMAX + 1];
	int fatti1[SIZEMAX + 1];
	int i, j;
	for (i = 1; i <= SIZEMAX; fatti[i++] = 0);


	for (i = 1; i <= SIZEMAX; i++)
	{
		for (j = 1; j <= SIZEMAX; fatti1[j++] = 0);
		fatti1[g1] = 2;				// per evitare di scorrere dal lato sbagliato
		scorri_(g2, i, fatti1);		// usa la funzione scorri_
									// salva i risultati in fatti
		for (j = 1; j <= SIZEMAX; j++) if (fatti1[j] == 1) fatti[j] = 1;
	};
	fatti[g2] = 1;
	//-------------- ora in fatti tutti i gruppi a dx di g2 valgono 1
	//-------------- scorriamo fatti e ripartiamo la struttura
	for (i = 1; i <= size(); i++)
	{
		if (fatti[i] == 0)
			for (j = 1; j <= 5; j++)
			{
				a->struttura[i][j] = struttura[i][j];
				b->struttura[i][j] = 0;
			}
		else
			for (j = 1; j <= 5; j++)
			{
				a->struttura[i][j] = 0;
				b->struttura[i][j] = struttura[i][j];
			};  // end else
	};  // end for i
 //  ----------------togliamo il legame che e' stato spezzato
	int r = 2;
	for (j = 2; j <= 5; j++) if (struttura[g1][j] != g2)
	{
		a->struttura[g1][r] = struttura[g1][j];
		r++;
	};
	while (r <= 5) a->struttura[g1][r++] = 0;
	if(a->tipo(g1) != 9)
		a->tipo(g1, 2);

	r = 2;
	for (j = 2; j <= 5; j++) if (struttura[g2][j] != g1)
	{
		b->struttura[g2][r] = struttura[g2][j];
		r++;
	};
	while (r <= 5) b->struttura[g2][r++] = 0;
	if(b->tipo(g2) != 9)
		b->tipo(g2, 2);
	//  ---------------- ripartiamo i doppi legami
	for (i = 1; i <= numole; i++)
	{
		if (fatti[ole[i][0]] == 0)  // vede uno dei gruppi del doppio legame
		{                   // e' nella molecola a
			a->addole(ole[i][0], ole[i][1]);
		}
		else
		{
			b->addole(ole[i][0], ole[i][1]);
		};
	};
	//  ---------------- ripartiamo gli =O
	for (i = 1; i <= numO; i++)
	{
		if (fatti[doppioO[i]] == 0) // =O e' su a
		{
			if(a->ischeto(doppioO[i]) == false)
				a->addcheto(doppioO[i]);
		}
		else if (b->ischeto(doppioO[i]) == false)
		{
			b->addcheto(doppioO[i]);
		}
	};
	// ---------------se la molecola era un radicale bis. fare l' olefina
	// if molcules was a radical add the olefins or transform the linear ether in ketone
	if (trova(2) != 0)
	{
		int p = trova(2);
		if (dist(p, g1) == 1 && tipo(g1) != 8) {
			a->tipo(p, 1);
			a->tipo(g1, 1);
			a->addole(p, g1);
		};
		if (dist(p, g2) == 1 && tipo(g2) != 8) {
			b->tipo(p, 1);
			b->tipo(g2, 1);
			b->addole(p, g2);
		};
		if (dist(p, g1) == 1 && tipo(g1) == 8) {
			a->tipo(p, 1);
			for (j = 1; j <= 5; j++)
				a->struttura[g1][j] = 0;
			a->tipo(g1, 1);
			if(a->ischeto(p) == false)
				a->addcheto(p);
		};
		if (dist(p, g2) == 1 && tipo(g2) == 8) {
			b->tipo(p, 1);
			for (j = 1; j <= 5; j++)
				b->struttura[g2][j] = 0;
			if(b->addcheto(p) == false)
				b->addcheto(p);
		};
	};
	if (a->tipo(g1) == 9)
	{
		a->tipo(g1, 1);
		if(a->ischeto(g1) == false)
			a->addcheto(g1);
	}
	if (b->tipo(g2) == 9)
	{
		b->tipo(g2, 1);
		if(b->ischeto(g2)==false)
			b->addcheto(g2);
	}
	// --------------------- distribute the cyclic -O-
	for (i = 1; i <= numete; i++)
	{
		if (fatti[ete[i][0]] == 0 && fatti[ete[i][1]] == 0)	// cyclyc -O- is all in a
			a->addetero(ete[i][0], ete[i][1]);
		if (fatti[ete[i][0]] == 1 && fatti[ete[i][1]] == 1)	// cyclyc -O- is all in b
			b->addetero(ete[i][0], ete[i][1]);
		if ((fatti[ete[i][0]] == 0 && fatti[ete[i][1]] == 1) || (fatti[ete[i][0]] == 1 && fatti[ete[i][1]] == 0))	// cyclyc -O- is bridging a and b
		{
			Molecola mol;
			int pos = 1;
			for (int k = 1; k <= SIZEMAX; k++)			// merge the two molecules
			{
				if (b->struttura[k][1] == 1 || b->struttura[k][1] == 2 || b->struttura[k][1] == 3 || b->struttura[k][1] == 4 || b->struttura[k][1] == 8 || b->struttura[k][1] == 9)
				{
					for (int l = 1; l <= 5; l++) mol.struttura[k][l] = b->struttura[k][l];
					pos++;
				}
				if (a->struttura[k][1] == 1 || a->struttura[k][1] == 2 || a->struttura[k][1] == 3 || a->struttura[k][1] == 4 || a->struttura[k][1] == 8 || b->struttura[k][1] == 9)
				{
					for (int l = 1; l <= 5; l++) mol.struttura[k][l] = a->struttura[k][l];
					pos++;
				}
			}
			mol.struttura[pos][1] = 8;				// add the -O- atom
			mol.struttura[pos][2] = ete[i][0];
			mol.struttura[pos][3] = ete[i][1];
			for (int l = 1; l <= 5; l++)			// add the bonds to that atom
			{
				if (mol.struttura[ete[i][0]][l] == 0)
				{
					mol.struttura[ete[i][0]][l] = pos;
					break;
				}
			}
			for (int l = 1; l <= 5; l++)			// add the bonds to that atom
			{
				if (mol.struttura[ete[i][1]][l] == 0)
				{
					mol.struttura[ete[i][1]][l] = pos;
					break;
				}
			}
			pos++;
			for (int k = 1; k <= SIZEMAX; k++)			// add the groups
			{
				if (b->struttura[k][1] == 5 || b->struttura[k][1] == 6 || b->struttura[k][1] == 7)
				{
					for (int l = 1; l <= 5; l++) mol.struttura[pos][l] = b->struttura[k][l];
					pos++;
				}
				if (a->struttura[k][1] == 5 || a->struttura[k][1] == 6 || a->struttura[k][1] == 7)
				{
					for (int l = 1; l <= 5; l++) mol.struttura[pos][l] = a->struttura[k][l];
					pos++;
				}
			}
			for (int k = 0; k < b->numole; k++)
				mol.addole(b->ole[k + 1][0], b->ole[k + 1][1]);
			for (int k = 0; k < a->numole; k++)
				mol.addole(a->ole[k + 1][0], a->ole[k + 1][1]);
			for (int k = 0; k < b->numO; k++)
				if(mol.ischeto(b->doppioO[k + 1]) == false)
					mol.addcheto(b->doppioO[k + 1]);
			for (int k = 0; k < b->numO; k++)
				if(mol.ischeto(b->doppioO[k + 1]) == false)
					mol.addcheto(b->doppioO[k + 1]);

			*a = mol;
			for (int k = 1; k <= SIZEMAX; k++)				// erase b
			{
				for (int l = 1; l <= 5; l++) b->struttura[k][l] = 0;
			}
		}
	}

	// ---------------si devono togliere le linee vuote dalla struttura
	a->noempty();
	b->noempty();
	return 1;
};

void Molecola::noempty(void)
{
	int conta[SIZEMAX + 1];
	for (int j = 0; j <= SIZEMAX; conta[j++] = 0);
	int c = 1;
	for (int j = 1; j <= SIZEMAX; j++)
		if (tipo(j) != 0) conta[j] = c++; // crea il vettore conta, che numera
										// i gruppi per le linee non vuote

	for (int r = 1; r <= SIZEMAX; r++)   // scorre tutte le righe
	{
		if (conta[r] == 0) continue;//solo quelle piene vengono copiate indietro
		struttura[conta[r]][1] = struttura[r][1];
		for (int cc = 2; cc <= 5; cc++)
		{
			struttura[conta[r]][cc] = conta[struttura[r][cc]];
		};
	};
	while (c <= SIZEMAX)  // continua con c a cancellare le righe rimaste
	{
		for (int k = 0; k <= 5; struttura[c][k++] = 0);
		c++;
	};
	int i;
	for (i = 1; i <= numole; i++)
	{
		ole[i][0] = conta[ole[i][0]];
		ole[i][1] = conta[ole[i][1]];
	};
	for (i = 1; i <= numete; i++)
	{
		ete[i][0] = conta[ete[i][0]];
		ete[i][1] = conta[ete[i][1]];
	};
	for (i = 1; i <= numO; i++)
		doppioO[i] = conta[doppioO[i]];

	// remove the gaps in the bonds
	for (int k = 1; k <= SIZEMAX; k++)
	{
		for (int l = 2; l <= 5; l++)
		{
			if (struttura[k][l] == 0) {
				for (int g = l; g <= 5; g++)
				{
					if (struttura[k][g] != 0)
					{
						struttura[k][l] = struttura[k][g];
						struttura[k][g] = 0;
					}
				}
			}
		}
	}
};

void Molecola::scorri(int dadove, int dist, int trovati[SIZEMAX + 1])
///////////////////////////////////////////////////////////////////////
//   rende i gruppi a distanza dist nel vettore trovati:
//   per quelli trovati la posizione corrispondente del vettore vale 1
///////////////////////////////////////////////////////////////////////
{
	for (int i = 0; i <= (SIZEMAX); trovati[i++] = 0);
	scorri_(dadove, dist, trovati);
	for (int i = 0; i <= (SIZEMAX); i++)
		if (trovati[i] != 1) trovati[i] = 0;
};

void Molecola::scorri_(int dadove, int dist, int trovati[SIZEMAX + 1])
///////////////////////////////////////////////////////////////////////
//   funzione ricorsiva utilizzata da scorri
///////////////////////////////////////////////////////////////////////
{
	if (dist == 0) trovati[dadove] = 1;
	else if (dist == 1)
	{
		for (int k = 1; k <= legami(dadove); k++)
			if (trovati[legato(dadove, k)] != 2)
				trovati[legato(dadove, k)] = 1;
	} // end if
	else
	{
		trovati[dadove] = 2;
		for (int k = 1; k <= legami(dadove); k++)
			if (trovati[legato(dadove, k)] != 2)
				scorri_(legato(dadove, k), (dist - 1), trovati);
	}; // end else
};   // end scorri_

int Molecola::dist(int pos1, int pos2)
{
	int vicini[SIZEMAX + 1];
	for (int dist = 1; dist <= SIZEMAX; dist++)
	{
		scorri(pos1, dist, vicini);
		if (vicini[pos2] == 1) return dist;
	};
	return -1;
};


/*
int Molecola::operator==(Molecola b)
//
// nota: non fa ancora etero e ole
//       vanno controllate dopo la sovrapposizione o quando?
//
{
	// check if the two molecules are C*-C=C-C and C=C-*C-C, in this case consider them equal
	if (size() == 4 && isole())
	{
		if (alreadyCheckingC4H713 == false)		// to avoid infinite loop that leads to stack over
		{
			if (bothC4H713(b))
			{
				return true;
			}
		}
	}

	if (specialMolecule != 0 || b.specialMolecule != 0)
	{
		if (specialMolecule == b.specialMolecule)
			return 1;
		else
			return 0;
	}

	if (size()         != b.size()        ) return 0;
	if (numberOfC()    != b.numberOfC()   ) return 0;
	if (numberOfH()    != b.numberOfH()   ) return 0;   
	if (numberOfO()    != b.numberOfO()	 ) return 0;
	if (numberOle()    != b.numberOle()	 ) return 0;
	if (numberEthero() != b.numberEthero()) return 0;
	if (numLinEther()  != b.numLinEther() ) return 0;
	if (numKeto()      != b.numKeto()	 ) return 0;
	if (numCrad()      != b.numCrad()	 ) return 0;
	if (numCOORad()    != b.numCOORad()	 ) return 0;
	if (numCORad()     != b.numCORad()	 ) return 0;
	if (numCOOH()      != b.numCOOH()	 ) return 0;

	struct nodeType
	{

		nodeType(int numBonds, std::vector<int> adjBonds)
		{
			type[0] = numBonds;
			std::sort(adjBonds.begin(), adjBonds.end());
			int index = 1;
			for (int i = 0; i < 4; i++)
			{
				if (adjBonds[i] != 0)
				{
					type[index] = adjBonds[i];
					index++;
				}
			}
		}

		bool operator == (nodeType other)
		{
			if (type == other.type)
				return true;
			return false;
		}

		void addNode(int id)
		{
			atomsIDs.push_back(id);
		}


		std::vector<int> type = { 0,0,0,0,0 };
		std::vector<int> atomsIDs = {};

		void print()
		{
			std::cout << type[0] << " " << type[1] << " " << type[2] << " " << type[3] << " " << type[4] << " " << atomsIDs.size() << std::endl;
		}

		int nodePriority()
		{
			return type[0] * 1e4 + type[1] * 1e3 + type[2] * 1e2 + type[3] * 1e1 + type[4];
		}
		
	};

	struct nodesList
	{
		nodesList()
		{

		}

		void addNode(int ID, int numBonds, std::vector<int> adj)
		{
			//std::cout << "ID = " << ID << ", numBonds = " << numBonds << ", adj = " <<
			//	adj[0] << " " << adj[1] << " " << adj[2] << " " << adj[3] << " " << std::endl;
			nodeType newNode(numBonds, adj);
			bool nodeFound = false;
			for (int i = 0; i < nodes.size(); i++)
			{
				if (nodes[i] == newNode)
				{
					nodes[i].addNode(ID);
					nodeFound = true;
					break;
				}
			}
			if (nodeFound == false)
			{
				int index = 0;
				for (int i = 0; i < nodes.size(); i++)
				{
					if (newNode.nodePriority() < nodes[i].nodePriority())
					{
						//std::cout << i << "    ";
						index = i;
						break;
					}
					if (i == nodes.size() - 1)
						index = nodes.size();
				}
				//std::cout << index << std::endl;
				nodes.insert(nodes.begin() + index, newNode);
			}


		}

		bool operator == (nodesList other)
		{
			if (other.nodes.size() != nodes.size())
				return false;

			//for (int i = 0; i < nodes.size(); i++)
			//{
			//	nodes[i].print();
			//	other.nodes[i].print();
			//	std::cout << std::endl;
			//}

			for (int i = 0; i < nodes.size(); i++)
			{
				if (nodes[i].atomsIDs.size() != other.nodes[i].atomsIDs.size())
					return false;
				if (!(nodes[i] == other.nodes[i])) 
				{
					return false;
				}
			}
			return true;
		}

		std::vector<nodeType> nodes = {};
	};

	// build the list of node types for the two molecules
	nodesList nodesList1;
	nodesList nodesList2;

	for (int i = 0; i < size(); i++)
	{
		int numBonds1 = numGraphConnections(i+1);
		int numBonds2 = b.numGraphConnections(i+1);

		std::vector<int> adj1 = { 0,0,0,0 };
		std::vector<int> adj2 = { 0,0,0,0 };

		for (int j = 0; j < 4; j++)
		{
			adj1[j] = numGraphConnections(struttura[i + 1][j + 2]);
			adj2[j] = b.numGraphConnections(b.struttura[i + 1][j + 2]);
		}
		//nodeType node1(numBonds1, adj1);
		//nodeType node2(numBonds2, adj2);

		nodesList1.addNode(i+1, numBonds1, adj1);
		nodesList2.addNode(i+1, numBonds2, adj2);
	}



	// DEBUG
	if (nodesList1 == nodesList2)
		return true;
	
	return false;
	// DEBUG

};   // end ==
*/


int Molecola::operator==(Molecola b)
//
// nota: non fa ancora etero e ole
//       vanno controllate dopo la sovrapposizione o quando?
//
{
	// check if the two molecules are C*-C=C-C and C=C-*C-C, in this case consider them equal
	if (size() == 4 && isole())
	{
		if (alreadyCheckingC4H713 == false)		// to avoid infinite loop that leads to stack over
		{
			if (bothC4H713( b))
			{
				return true;
			}
		}
	}

	if (specialMolecule != 0 || b.specialMolecule != 0)
	{
		if (specialMolecule == b.specialMolecule)
			return 1;
		else
			return 0;
	}


	if (size() != b.size()) return 0;
	int ngruppi = size();
	int fattia[SIZEMAX + 1];
	int fattib[SIZEMAX + 1];
	int corr[SIZEMAX + 1];
	int posa = 1;
	int uguali = 0;
	for (int posb = 1; posb <= ngruppi; posb++)    // cerca di sovrapporre
	{                                           // A partendo da posa=1
		for (int k = 1; k <= SIZEMAX; fattia[k++] = 0);  // e B da posb variabile
		for (int k = 1; k <= SIZEMAX; fattib[k++] = 0);
		for (int k = 1; k <= SIZEMAX; corr[k++] = 0);
		if (confronta(this, &b, posa, posb, fattia, fattib, corr) == 1)
		{
			uguali = 1;
			for (int e = 1; e <= numete; e++)
				if (!b.isetero(corr[ete[e][0]], corr[ete[e][1]]))
				{
					uguali = 0;
					break;    // esce dal ciclo sugli etero
				};
			if (uguali == 1) break;  // e' 1 se ha passato tutti gli etero
		};
	};
	if (uguali == 1) return 1;
	else return 0;
};   // end ==



int Molecola::isIsomer(Molecola* other)
{
	if (numberOfC() != other->numberOfC()) return 0;
	if (numberOfH() != other->numberOfH()) return 0;
	if (numberOle() != other->numberOle()) return 0;
	if (numberEthero() != other->numberEthero()) return 0;
	if (numKeto() != other->numKeto()) return 0;
	if (numCrad() != other->numCrad()) return 0;
	if (numCOORad() != other->numCOORad()) return 0;
	if (numCOOH() != other->numCOOH()) return 0;

	return 1;
}

int Molecola::confronta(Molecola* a, Molecola* b, int posa, int posb,
	int fattia[SIZEMAX + 1], int fattib[SIZEMAX + 1],
	int corr[SIZEMAX + 1])
{
	if (a->tipo(posa) != b->tipo(posb)) return 0;
	if (a->legami(posa) != b->legami(posb)) return 0;
	if (a->ischeto(posa) != b->ischeto(posb)) return 0;
	if (a->isole(posa) != b->isole(posb)) return 0;
	if (a->isetero(posa) != b->isetero(posb)) return 0;

	fattia[posa] = 1;
	fattib[posb] = 1;
	corr[posa] = posb;

	int i;

	for (i = 1; i <= a->legami(posa); i++) if (fattia[a->legato(posa, i)] != 1)
	{
		int trovatoj = 0;
		for (int j = 1; j <= b->legami(posb); j++) if (fattib[b->legato(posb, j)] != 1)
		{
			int nexta = a->legato(posa, i);
			int nextb = b->legato(posb, j);
			if (a->isole(posa, nexta) != b->isole(posb, nextb)) continue;
			if (confronta(a, b, nexta, nextb, fattia, fattib, corr) == 1)
			{
				trovatoj = 1;
				break;
			};
		};
		if (trovatoj != 1) return 0;
	};

	return 1;

};   // end confronta

ostream& operator<<(ostream& stream, Molecola& obj)
//////////////////////////////////////////////////////////////////////////////
//  SCRIVE la molecola su uno stream
//  come matrice o come formula,
//  a secondo di cosa valga obj.output (MATRICE o FORMULA)
//////////////////////////////////////////////////////////////////////////////
{

	if (obj.specialMolecule != 0)
	{
		//1 = N2
		//2 = O2
		//3 = O
		//4 = OH
		//5 = HO2
		//6 = H2O
		//7 = H2O2
		//8 = H
		//9 = H2
		switch (obj.specialMolecule)
		{
		case 1:
			stream << "N2";
			break;
		case 2:
			stream << "O2";
			break;
		case 3:
			stream << "O";
			break;
		case 4:
			stream << "OH";
			break;
		case 5:
			stream << "HO2";
			break;
		case 6:
			stream << "H2O";
			break;
		case 7:
			stream << "H2O2";
			break;
		case 8:
			stream << "H";
			break;
		case 9:
			stream << "H2";
			break;
		case 10:
			stream << "CO";
			break;
		default:
			stream << "ERR";
			cerr << "Special molecule not recognized!";
		}
	}
	else if (obj.output == MATRICE)
		//----------------------------------------------------------------------------
	{
		for (int i = 1; i <= obj.size(); i++)
		{
			for (int j = 1; j <= 5; j++)
			{
				stream << obj.struttura[i][j] << ' ';
			};     // end for j
			stream << endl;
		};        // end for i
		for (int i = 1; i <= obj.numole; i++)
		{
			stream << "Double bond between: " << obj.ole[i][0]
				<< " e " << obj.ole[i][1] << endl;
		};
		for (int i = 1; i <= obj.numete; i++)
		{
			stream << "Cyclic oxygen between: " << obj.ete[i][0]
				<< " e " << obj.ete[i][1] << endl;
		};
		for (int i = 1; i <= obj.numO; i++)
		{
			stream << "Oxygen =O bonded to: " << obj.doppioO[i] << endl;
		};
	}
	else if (obj.output == FORMULA)
		//----------------------------------------------------------------------------
	{
		int fatti[SIZEMAX + 1];
		for (int i = 1; i <= SIZEMAX; fatti[i++] = 0);
		char str[80];
		ostrstream outs(str, sizeof(str));  // prepara la stringa str per
										   // essere riempita come uno stream
										   // da Formula
		obj.Formula(outs, 1, fatti);
		outs << '\0';
		stream << str;
	}
	else   // obj.output FORMULA2
//----------------------------------------------------------------------------
	{
		int fatti[SIZEMAX + 1];
		for (int i = 1; i <= SIZEMAX; fatti[i++] = 0);
		char str[80];
		ostrstream outs(str, sizeof(str));  // prepara la stringa str per
										   // essere riempita come uno stream
										   // da Formula
		obj.Formula(outs, 1, fatti);
		outs << '\0';

		char linea1[80];
		char linea2[80];
		char ch = ' ';
		int j = 0, punti = 0, pass = 0, print = 0;
		while (j <= strlen(str))
		{
			if (str[j] == '^')
				if (pass == 0)
				{
					punti++;
					pass = 1;
					print = 1;
					ch = '-';
					j++;
					continue;
				}
				else
				{
					punti++;
					pass = 0;
					print = 1;
					ch = ' ';
					j++;
					continue;
				};     // end if str==.
			linea1[j - punti] = str[j];
			if (print == 1)
			{
				print = 0;
				if (pass == 1) linea2[j - punti] = '`';
				else linea2[j - punti] = '\'';
			}
			else
				linea2[j - punti] = ch;
			j++;
		};
		linea1[j] = '\0';
		linea2[j] = '\0';
		// qui manca da scrivere l' ossigeno al punto giusto tra i ____
		cout << linea1;
		if (obj.isetero()) { cout << endl << linea2; };
	};
	return stream;
};

void Molecola::Formula(ostream& stream, int pos, int fatti[])
//
// nota: non fa ancora etero
//
{
	//-------------------scrive la sigla del gruppo
	char nome[10];
	switch (tipo(pos))
	{
	case 1: strcpy(nome, "C"); break;
	case 2: strcpy(nome, "*C"); break;
	case 3: strcpy(nome, "COO*"); break;
	case 4: strcpy(nome, "COOH"); break;
	case 8: strcpy(nome, "O"); break;
	case 9: strcpy(nome, "CO*"); break;

	};
	if (isetero(pos)) stream << '^';
	stream << nome;
	if (ischeto(pos)) stream << 'O';
	//-------------------segna fatta la posizione
	fatti[pos] = 1;
	int i, vicini = 0, maxvicino;
	for (i = 1; i <= legami(pos); i++)
		if (fatti[legato(pos, i)] == 0) {
			vicini++;  // conta i vicini da fare per
			maxvicino = i;
		}; //calcolare le parentesi
//-------------------scorre i vicini non fatti
	for (i = 1; i <= legami(pos); i++)
	{
		if (fatti[legato(pos, i)] == 1) continue;
		//------------------------------------- e se ce ne sono...
		if (i != maxvicino)
		{
			int pos2 = legato(pos, i);
			if (isole(pos, pos2)) { stream << "(="; }
			else { stream << "(-"; };
			Formula(stream, legato(pos, i), fatti);
			stream << ")";
		}
		else
		{
			int pos2 = legato(pos, i);
			if (isole(pos, pos2)) { stream << "="; }
			else { stream << "-"; };
			Formula(stream, legato(pos, i), fatti);
		};

	};  // il continue manda qui
};

istream& operator>>(istream& stream, Molecola& obj)
//////////////////////////////////////////////////////////////////////////////
//  Read the molecules from a stream            ITA: LEGGE la molecola da uno stream
//////////////////////////////////////////////////////////////////////////////
{
	int j = 1;   // count the number of rows         ITA: conta le righe
	int num;   // read the number                  ITA: legge il numero
	char ch;   // read until the end of the line   ITA: legge fino fine linea
	while (stream)
	{
		for (int i = 1; i <= 5; i++)          // read the first 5 numbers of the j-th line
		{
			stream >> num;
			obj.struttura[j][i] = num;
		}

		do { stream.get(ch); }              // flush everything on the line after the first 5 numbers
		while (ch != '\n');
		//----------------------check if the group is 5, 6 or 7 before j++     ITA: controlla se il gruppo e' 5 o 6 prima di j++
		if (obj.struttura[j][1] == 5)
		{
			obj.addole(obj.struttura[j][2], obj.struttura[j][3]);
		}
		else if (obj.struttura[j][1] == 6)
		{
			obj.addetero(obj.struttura[j][2], obj.struttura[j][3]);
		}
		else if (obj.struttura[j][1] == 7)
		{
			obj.addcheto(obj.struttura[j][2]);
		}
		else j++;
	}   // end of the file
 // fill the remainig with zeros            ITA: riempie il resto di zeri
	for (int i = j; i <= SIZEMAX; i++)
		for (int k = 1; k <= 5; obj.struttura[i][k++] = 0);
	return stream;
};

int Molecola::numberOfC()
{
	int numberOfCarbons = 0;
	if (isSpecialMolecule() == 0)
	{
		for (int i = 1; i < SIZEMAX + 1; i++)
		{
			if (struttura[i][1] == 1 || struttura[i][1] == 2 || struttura[i][1] == 3 || struttura[i][1] == 4 || struttura[i][1] == 9)
				numberOfCarbons++;
		}
	}
	else
	{
		switch (specialMolecule)
		{
		case 10:
			numberOfCarbons = 1;
			break;
		}
	}
	//int specialMolecule;	// define if the molecule is a special one :	1 = N2
							//												2 = O2
							//												3 = O
							//												4 = OH
							//												5 = HO2
							//												6 = H2O
							//												7 = H2O2
							//												8 = H
							//												9 = H2
							//												10= CO

	return numberOfCarbons;
}

int Molecola::numberOfH() // TO DO: consider all the possible groups
{
	// LIST OF SPECIAL MOLECULES
	//	1 = N2
	//	2 = O2
	//	3 = O
	//	4 = OH
	//	5 = HO2
	//	6 = H2O
	//	7 = H2O2
	//	8 = H
	//	9 = H2
	
	int numberOfHydrogens = 0;
	if (isSpecialMolecule() == 0)
	{
		for (int i = 1; i < SIZEMAX + 1; i++)
		{
			if ((tipo(i) > 4 && tipo(i) != 9) || tipo(i) < 1) continue;
			numberOfHydrogens += numH(i);
		}
	}
	else
	{
		switch (specialMolecule)
		{
		case 4:
			numberOfHydrogens = 1;
			break;
		case 5:
			numberOfHydrogens = 1;
			break;
		case 6:
			numberOfHydrogens = 2;
			break;
		case 7:
			numberOfHydrogens = 2;
			break;
		case 8:
			numberOfHydrogens = 1;
			break;
		case 9:
			numberOfHydrogens = 2;
			break;
		}
	}


	return numberOfHydrogens;
}

int Molecola::numberOfO()				// number of oxygens in the molecule
{

	// LIST OF SPECIAL MOLECULES
	//	1 = N2
	//	2 = O2
	//	3 = O
	//	4 = OH
	//	5 = HO2
	//	6 = H2O
	//	7 = H2O2
	//	8 = H
	//	9 = H2

	int numberOfOxygens = 0;
	if (isSpecialMolecule() == 0)
	{
		for (int i = 1; i < SIZEMAX + 1; i++)
		{
			if (tipo(i) == 3) numberOfOxygens += 2;
			if (tipo(i) == 4) numberOfOxygens += 2;
			if (tipo(i) == 8) numberOfOxygens += 1;
			if (tipo(i) == 9) numberOfOxygens += 1;
		}
		numberOfOxygens += numete;
		numberOfOxygens += numO;
	}
	else
	{
		switch (specialMolecule)
		{
		case 2:
			numberOfOxygens = 2;
			break;
		case 3:
			numberOfOxygens = 1;
			break;
		case 4:
			numberOfOxygens = 1;
			break;
		case 5:
			numberOfOxygens = 2;
			break;
		case 6:
			numberOfOxygens = 1;
			break;
		case 7:
			numberOfOxygens = 2;
			break;
		case 10:
			numberOfOxygens = 1;
			break;
		}
	}
	return numberOfOxygens;
}


int Molecola::numDoubleBonds(int pos)	// return the number of double bonds the atom is involved into
{
	int num = 0;
	for (int i = 1; i <= numole; i++)
	{
		if (ole[i][0] == pos || ole[i][1] == pos) num++;
	};
	return num;
}

int Molecola::numberOle()				// return the number of double bond present in the molecule
{
	return numole;
}

int Molecola::numberEthero()				// retunr the number of ethero present in the molecule
{
	return numete;
}

int Molecola::numKeto()					// return the number of cheto groups in the molecule
{
	return numO;
}

int Molecola::numLinEther()
{
	int num = 0;
	for (int i = 1; i <= SIZEMAX; i++)
		if (tipo(i) == 8)  num++;
	return num;
}

int Molecola::numCrad()					// retunr the number of C* present in the molecule
{
	int num = 0;
	for (int i = 1; i <= SIZEMAX; i++)
		if (tipo(i) == 2)  num++;
	return num;
}

int Molecola::numCOORad()				// retunr the number of COO* present in the molecule
{
	int num = 0;
	for (int i = 1; i <= SIZEMAX; i++)
		if (tipo(i) == 3)  num++;
	return num;
}

int Molecola::numCORad()				// retunr the number of CO* present in the molecule
{
	int num = 0;
	for (int i = 1; i <= SIZEMAX; i++)
		if (tipo(i) == 9)  num++;
	return num;
}

int Molecola::numCOOH()					// retunr the number of COOH present in the molecule
{
	int num = 0;
	for (int i = 1; i <= SIZEMAX; i++)
		if (tipo(i) == 4)  num++;
	return num;
}

bool Molecola::isLinear()
{
	for (int i = 1; i <= SIZEMAX; i++)
		if (legami(i) > 2) return false;
	return true;
}

int Molecola::posCGroup(int group)		// return the position of the specified carbon group (C*, COO* and COOH), return 0 if the group is not present in the moelcule 
{
	int pos = 0;
	for (int i = 1; i <= SIZEMAX; i++)
		if (tipo(i) == group)  pos = i;
	return pos;
}

int Molecola::posCrad()					// return the position of the C* radical, return 0 if no radical is present
{
	return posCGroup(2);
}

int Molecola::posCOOrad()				// return the position of the COO* radical, return 0 if no radical is present
{
	return posCGroup(3);
}

int Molecola::posCOrad()				// return the position of the CO* radical, return 0 if no radical is present
{
	return posCGroup(9);
}

int Molecola::posCOOH()					// return the position of the COOH group, return 0 if no COOH group is present
{
	return posCGroup(4);
}

int Molecola::posOle()
{
	if (numole == 0)
		return 0;
	else
		return std::min(ole[1][0], ole[1][1]);
}

std::vector<int> Molecola::posOle(int num)
{
	std::vector<int> pos = { 0, 0 };
	if (numole == 0)
	{
		std::cerr << "WARNING: posOle called on a non olefin!" << std::endl;
		return pos;
	}
	if (num > numole)
	{
		std::cerr << "ERROR: posOle called num>numberOle!" << std::endl;
		return pos;
	}
	pos[0] = ole[num][0];
	pos[1] = ole[num][1];
	return pos;
}

int Molecola::posKeto()					// return the position of =O, return 0 if no =O is present
{
	if (numO == 0)
		return 0;
	else
		return doppioO[1];
}

std::vector<int> Molecola::posEthero()				// return the two positions of the ether, return [0,0] if no ether is present
{
	std::vector<int> pos(2);
	if (numete == 0)
	{
		pos[0] = 0;
		pos[1] = 0;
		return pos;
	}
	else
	{
		pos[0] = std::min(ete[1][0], ete[1][1]);
		pos[1] = std::max(ete[1][0], ete[1][1]);
		return pos;
	}
}

std::vector<int> Molecola::posOOHinPOOH2()
{
	std::vector<int> pos = { 0, 0 };
	if (numCOOH() != 2 && numCrad() != 0)
	{
		std::cerr << "ERROR: posOOHinPOOH2 called on a non-P(OOH)2 molecule" << std::endl;
		return pos;
	}
	Molecola mol = *this;
	int pos1 = mol.posCOOH();
	mol.removeOOH(pos1);
	int pos2 = mol.posCOOH();
	pos[0] = std::min(pos1, pos2);
	pos[1] = std::max(pos1, pos2);
	return pos;
}

std::vector<int> Molecola::posOOHinROOH2()
{
	std::vector<int> pos = { 0, 0 };
	if (kindOfSPecies() != ROOH2_)
	{
		std::cerr << "ERROR: posOOHinROOH2 called on a non-R(OOH)2 molecule" << std::endl;
		return pos;
	}
	Molecola mol = *this;
	int pos1 = mol.posCOOH();
	mol.removeOOH(pos1);
	int pos2 = mol.posCOOH();
	pos[0] = std::min(pos1, pos2);
	pos[1] = std::max(pos1, pos2);
	return pos;
}

void Molecola::makeCH3()				// make the molecule into CH3
{
	init();
	struttura[1][1] = 2;
	struttura[1][2] = 0;
	struttura[1][3] = 0;
	struttura[1][4] = 0;
	struttura[1][5] = 0;
}
void Molecola::makeCH4()				// make the molecule into CH4
{
	init();
	struttura[1][1] = 1;
	struttura[1][2] = 0;
	struttura[1][3] = 0;
	struttura[1][4] = 0;
	struttura[1][5] = 0;
}
void Molecola::makeC2H5()			// make the molecule into C2H5
{
	init();
	struttura[1][1] = 2;
	struttura[1][2] = 2;
	struttura[1][3] = 0;
	struttura[1][4] = 0;
	struttura[1][5] = 0;
	struttura[2][1] = 1;
	struttura[2][2] = 1;
	struttura[2][3] = 0;
	struttura[2][4] = 0;
	struttura[2][5] = 0;
}
void Molecola::makeC2H6()			// make the molecule into C2H6
{
	init();
	struttura[1][1] = 1;
	struttura[1][2] = 2;
	struttura[1][3] = 0;
	struttura[1][4] = 0;
	struttura[1][5] = 0;
	struttura[2][1] = 1;
	struttura[2][2] = 1;
	struttura[2][3] = 0;
	struttura[2][4] = 0;
	struttura[2][5] = 0;
}

void Molecola::makeCH3OO()
{
	init();
	struttura[1][1] = 3;
	struttura[1][2] = 0;
	struttura[1][3] = 0;
	struttura[1][4] = 0;
	struttura[1][5] = 0;
}

void Molecola::makeCH3OOH()
{
	init();
	struttura[1][1] = 4;
	struttura[1][2] = 0;
	struttura[1][3] = 0;
	struttura[1][4] = 0;
	struttura[1][5] = 0;
}

void Molecola::makeC4H72_1()
{
	init();
	struttura[1][1] = 2; struttura[1][2] = 2; struttura[1][3] = 0;
	struttura[2][1] = 1; struttura[2][2] = 1; struttura[2][3] = 3;
	struttura[3][1] = 1; struttura[3][2] = 2; struttura[3][3] = 4;
	struttura[4][1] = 1; struttura[4][2] = 3; struttura[4][3] = 0;
	addole(2, 3);
}

void Molecola::makeC4H71_3()
{
	struttura[1][1] = 1; struttura[1][2] = 2; struttura[1][3] = 0;
	struttura[2][1] = 1; struttura[2][2] = 1; struttura[2][3] = 3;
	struttura[3][1] = 2; struttura[3][2] = 2; struttura[3][3] = 4;
	struttura[4][1] = 1; struttura[4][2] = 3; struttura[4][3] = 0;
	addole(1, 2);
}

int Molecola::isSpecialMolecule()
{
	return specialMolecule;
}

std::vector<std::vector<int>> Molecola::listOfCCBonds()
{
	std::vector<std::vector<int>> list;

	for (int i = 1; i < SIZEMAX+1; i++)
	{
		if (struttura[i][1] == 1 || struttura[i][1] == 2 || struttura[i][1] == 3 || struttura[i][1] == 4 || struttura[i][1] == 9) // check only for carbon atoms
		{
			for (int j = 2; j < 6; j++)	// check all the possible bonds
			{
				if (struttura[i][j] != 0)	// if there is a bond
				{
					// check if the bond is not already present
					int at = struttura[i][j];
					bool isPresent = false;
					for (int k = 0; k < list.size(); k++)
						if ((list[k][0] == i && list[k][1] == at) || (list[k][0] == at && list[k][1] == i))
							isPresent = true;
					// if it is not present add it
					if (!isPresent)
					{
						std::vector<int> bond = { i, at };
						list.push_back(bond);
					}
				}
			}
		}
	}
	return list;
}

void Molecola::fix()
{
	if (isLinear() && size() > 1)
	{
		// DEBUG
		int asd = size();
		//std::cout << " this molecule " << *this << " has been fixed into ";
		int map[SIZEMAX + 1];	// maps of the position of the atoms in the old structure and in the new one (the index is the position in the old structure
								// and the value is the position in the new structure)
		for (int i = 0; i <= SIZEMAX; i++)
			map[i] = 0;

		int newStruttura[SIZEMAX + 1][5 + 1];
		for (int i = 0; i <= SIZEMAX; i++)
		{
			for (int j = 0; j <= 5; j++) 
				newStruttura[i][j] = 0;
		}

		int alreadyTaken[SIZEMAX + 1];
		for (int i = 0; i <= SIZEMAX; i++)
			alreadyTaken[i] = 0;

		// find the positions of the 2 terminal atoms
		bool isFirstFound = false;
		int firstAtomPos = 0;
		for (int i = 1; i <= size(); i++)
		{
			if (isFirstFound && legami(i) == 1)
			{
				map[i] = size();
				// set this terminal atom as last atom
				for (int j = 1; j <= 5; j++)
					newStruttura[size()][j] = struttura[i][j];
				alreadyTaken[i] = 1;
			}
			if (!isFirstFound && legami(i) == 1)
			{
				map[i] = 1;
				// set this terminal atom as first atom
				for (int j = 1; j <= 5; j++)
					newStruttura[1][j] = struttura[i][j];
				isFirstFound = true;
				alreadyTaken[i] = 1;
				firstAtomPos = i;
			}
			
		}

		// fill the intermediate atoms

		int lastTakenAtom = firstAtomPos;
		for (int i = 2; i <= size()-1; i++) // fill all the atom of the new structure
		{
			int bondedPos = 0;
			// scan for the bonded atom
			for (int k = 1; k <= size(); k++)
			{
				if (alreadyTaken[k] == 0)
				{
					for (int j = 2; j <= 5; j++)
					{
						if (struttura[k][j] == lastTakenAtom)
						{
							bondedPos = k;
						}
					}
				}
			}
			if (bondedPos == 0)
				std::cout << "ERROR: in Molecola.fix() bondedPos = 0!" << std::endl;
			for(int j = 1; j <= 5; j++)
				newStruttura[i][j] = struttura[bondedPos][j];
			map[bondedPos] = i;
			lastTakenAtom = bondedPos;
			alreadyTaken[bondedPos] = 1;
		}

		// change the index of the bonded atoms
		for (int i = 1; i <= size(); i++)
		{
			for (int j = 2; j <= 5; j++)
				newStruttura[i][j] = map[newStruttura[i][j]];
		}
		
		// update the old structure
		for (int i = 1; i <= size(); i++)
		{
			for (int j = 1; j <= 5; j++)
				struttura[i][j] = newStruttura[i][j];
		}

		// update the bonds of the groups
		for (int i = 1; i <= numole; i++)
		{
			ole[i][0] = map[ole[i][0]];
			ole[i][1] = map[ole[i][1]];
		}
		for (int i = 1; i <= numete; i++)
		{
			ete[i][0] = map[ete[i][0]];
			ete[i][1] = map[ete[i][1]];
		}
		for (int i = 1; i <= numO; i++)
		{
			doppioO[i] = map[doppioO[i]];
		}

		// DEBUG
		//std::cout << *this << std::endl;
	}
}

int Molecola::areBonded(int pos1, int pos2)
{
	// check if atom 1 is bonded to atom 2
	bool bond12 = false;
	for (int k = 2; k <= 5; k++)
		if (struttura[pos1][k] == pos2) bond12 = true;

	// check if atom 2 is bonded to atom 1
	bool bond21 = false;
	for (int k = 2; k <= 5; k++)
		if (struttura[pos1][k] == pos2) bond21 = true;

	if (bond12 == false || bond21 == false)
		return 0;
	
	return 1;
}

void Molecola::printStruttura()
{
	for (int i = 0; i < size(); i++)
	{
		for (int j = 0; j < 5; j++)
		{
			std::cout << struttura[i + 1][j + 1] << "  ";
		}
		std::cout << std::endl;
	}
}

void Molecola::printDoppioO()
{
	for (int i = 0; i < numO; i++)
	{
		std::cout << doppioO[i + 1] << "   ";
	}
	std::cout << std::endl;
}

species Molecola::kindOfSPecies()
{
	int numberR = numCrad();
	int numberCOO = numCOORad();
	int numberCOOH = numCOOH();
	int numberEther = numberEthero();
	int numberOLE = numberOle();
	int numberCO = numKeto();
	int specMol = isSpecialMolecule();
	int numLinEth = numLinEther();
	int numCOR = numCORad();


	std::stringstream nameMol;
	nameMol << *this;
	int asdf = 0;
	//if (nameMol.str() == "C=C-CO-C-C")
	//{
	//	printStruttura();
	//	std::cout << std::endl;
	//	printDoppioO();
	//	std::cout << std::endl;
	//	std::cout  
	//		<< "numberR     "<< numCrad()			  << std::endl
	//		<< "numberCOO   "<< numCOORad()		  << std::endl
	//		<< "numberCOOH  "<< numCOOH()			  << std::endl
	//		<< "numberEther "<< numberEthero()	  << std::endl
	//		<< "numberOLE   "<< numberOle()		  << std::endl
	//		<< "numberCO    "<< numKeto()			  << std::endl
	//		<< "specMol     "<< isSpecialMolecule() << std::endl
	//		<< "numLinEth   "<< numLinEther()       << std::endl
	//		<< "numCOR      "<< numCORad() << std::endl << std::endl;
	//
	//}

	if (specMol != 0)																	 // molecule is special molecule
		return special_;
	else if (numberR + numberCOO + numberCOOH + numberEther + numberOLE + numberCO + numLinEth + numCOR == 0) // molecule is alkane
		return fuel_;
	else if (numberR == 1 &&
		numberCOO + numberCOOH + numberEther + numberOLE + numberCO + numLinEth + numCOR == 0) // molecule is R
		return R_;
	else if (numberCOO == 1 &&
		numberR + numberCOOH + numberEther + numberOLE + numberCO + numLinEth + numCOR == 0)	// molecule is ROO
		return ROO_;
	else if (numberCOOH == 1 && numberR == 1 &&
		numberCOO + numberEther + numberOLE + numberCO + numLinEth + numCOR == 0)	// molecule is QOOH
		return QOOH_;
	else if (numberCOO == 1 && numberCOOH == 1 &&
		numberR + numberEther + numberOLE + numberCO + numLinEth + numCOR == 0)	// molecule is OOQOOH
		return OOQOOH_;
	else if (numberOLE == 1 &&
		numberR + numberCOO + numberCOOH + numberEther + numberCO + numLinEth + numCOR == 0)	// molecule is olefin
		return OLE_;
	else if (numberCO == 1 &&
		numberR + numberCOO + numberCOOH + numberEther + numberOLE + numLinEth + numCOR == 0)	// molecule is CO
		return CO_;
	else if (numberEther == 1 &&
		numberR + numberCOO + numberCOOH + numberOLE + numberCO + numLinEth + numCOR == 0)	// molecule is an ether
		return cEth_;
	else if (numberCO == 1 && numberR == 1 &&
		numberCOO + numberCOOH + numberEther + numberOLE + numLinEth + numCOR == 0)	// molecule is an RO
		return RO_;
	else if (numberCO == 1 && numberCOOH == 1 &&
		numberR + numberCOO + numberEther + numberOLE + numLinEth + numCOR == 0)	// molecule is KHP
		return KHP_;
	else if (numberCOOH == 1 &&
		numberR + numberCOO + numberEther + numberOLE + numberCO + numLinEth + numCOR == 0)	// molecule is ROOH
		return ROOH_;
	else if (numberCOOH == 2 &&
		numberR + numberCOO + numberEther + numberOLE + numberCO + numLinEth + numCOR == 0)	// molecule is R(OOH)2
		return ROOH2_;
	else if (numberCOOH == 2 && numberR == 1 &&
		numberCOO + numberEther + numberOLE + numberCO + numLinEth + numCOR == 0)	// molecule is P(OOH)2
		return POOH2_;
	else if (numberCOOH == 1 && numberOLE == 1 &&
		numberR + numberCOO + numberEther + numberCO + numLinEth + numCOR == 0)	// molecule is ole-OOH
		return oleOOH_;
	else if (numberEther == 1 && numberCOOH == 1 &&
		numberR + numberCOO + numberOLE + numberCO + numLinEth + numCOR == 0)	// molecule is an ether-OOH
		return cEthOOH_;
	else if (numberR == 1 && numberOLE == 1 &&
		numberCOO + numberCOOH + numberEther + numberCO + numLinEth + numCOR == 0) // molecule is oleR
		return oleR_;
	else if (numberOLE == 1 && numberCO == 1 &&
		numberR + numberCOO + numberCOOH + numberEther + numLinEth + numCOR == 0)	// molecule is oleCO
		return oleCO_;
	else if (numberEther == 1 && numberR == 1 &&
		numberCOO + numberCOOH + numberOLE + numberCO + numLinEth + numCOR == 0)	// molecule is an etherR
		return cEthR_;
	else if (numberCO == 1 && numberEther == 1 &&
		numberR + numberCOO + numberCOOH + numberOLE + numLinEth + numCOR == 0)	// molecule is etherCO
		return cEthCO_;
	else if (numberCO == 1 && numLinEth == 1 && numberR == 1 &&
		numberCOO + numberCOOH + numberOLE + numberEther + numCOR == 0)	// molecule is linEtherRO
		return lEthRO_;
	else if (numberOLE == 1 && numCOR == 1 &&
		numLinEth + numberCO + numberR + numberCOO + numberCOOH + numberEther == 0)	// molecule is alkenyl RO
		return alkRO_;
	else
		return unidentified_;
}


bool isElementPresent(std::vector<int>* vec, int value)
{
	for (int i = 0; i < vec->size(); i++)
		if ((*vec)[i] == value)
			return true;
	return false;
}

std::vector<int> Molecola::mainChain()
{
	if (size() == 1)
	{
		std::vector<int> a = { 1 };
		return a;
	}
	// prioritize carbons that have a group attached to them or that are nerby
	// a group to avoid putting the group into methyl or ethyl branching
	// ie in this case we want this main chain C-C-C(-C)-CO not this one C-C-C(-C)-CO
	//                                         ^ ^ ^     ^               ^ ^ ^  ^

	std::vector<int> terminalCs(0);		// the list of terminal carbons
	for (int i = 0; i < size(); i++)
		if (legami(i+1) == 1)
			terminalCs.push_back(i+1);

	// find the ends of the longest chain
	int length = 0;
	int firstEnd = 0;
	int secondEnd = 0;
	bool terminalInvolvingGroup = false;	// flag that indicates if a group containing terminal cs is in hte last considered pair of terminal cs
	bool terminalNearGroup = false;

	std::stringstream asdf;
	asdf << *this;
	if (asdf.str() == "C-C(-C)-C(-C)-CO")
		int fdsa = 0;

	for (int i = 0; i < terminalCs.size()-1; i++)
	{
		for (int j = i + 1; j < terminalCs.size(); j++)
		{
			int isgi = isGroup(terminalCs[i]);
			int isgj = isGroup(terminalCs[j]);
			int ingi = isNearGroup(terminalCs[i]);
			int ingj = isNearGroup(terminalCs[j]);
			int d = dist(terminalCs[i], terminalCs[j]) + 1;
			if (d > length)
			{
				length = d;
				firstEnd = terminalCs[i];
				secondEnd = terminalCs[j];
				if (isGroup(terminalCs[i]) || isGroup(terminalCs[j]))
					terminalInvolvingGroup = true;
				else
					terminalInvolvingGroup = false;

				if (isNearGroup(terminalCs[i]) || isNearGroup(terminalCs[j]))
					terminalNearGroup = true;
				else
					terminalNearGroup = false;
			}
			if (d == length && terminalInvolvingGroup == false)
			{
				if (isGroup(terminalCs[i]) || isGroup(terminalCs[j]))
				{
					firstEnd = terminalCs[i];
					secondEnd = terminalCs[j];
					terminalInvolvingGroup = true;
					if (isNearGroup(terminalCs[i]) || isNearGroup(terminalCs[j]))
						terminalNearGroup = true;
				}
				else if (terminalNearGroup == false)
				{
					if (isNearGroup(terminalCs[i]) || isNearGroup(terminalCs[j]))
					{
						firstEnd = terminalCs[i];
						secondEnd = terminalCs[j];
						terminalNearGroup = true;
					}
				}
			}
		}
	}

	// find all the element of the longest chain
	std::vector<int> elements(length,-1);
	bool chainFound = false;
	//std::vector<int> branchingLevel(length, 0);
	std::vector<int> branchingLevel(size(), 0);
	int pos = 0;
	elements[pos] = firstEnd;
	pos++;
	int currentAtom = firstEnd;
	while (chainFound == false)
	{
		int nextAtom = 0;
		while (true)
		{
			nextAtom = struttura[currentAtom][branchingLevel[currentAtom-1] + 2];
			if (isElementPresent(&elements, nextAtom) == false && nextAtom != 0)
				break;
			branchingLevel[currentAtom - 1] ++;
			if (branchingLevel[currentAtom - 1] == 4)
				branchingLevel[currentAtom - 1] = 0;
		}
		 
		if (legami(currentAtom) > 1 &&  branchingLevel[currentAtom - 1] < legami(currentAtom))
		{
			branchingLevel[currentAtom - 1] ++;
			if (branchingLevel[currentAtom - 1] == 4)
				branchingLevel[currentAtom - 1] = 0;
		}
		elements[pos] = nextAtom;
		pos++;
		currentAtom = nextAtom;

		if (legami(currentAtom) == 1)
		{
			//if (pos == length)
			//	break;
			//else if (pos < length )
			//{
			//	pos = 1;
			//	currentAtom = elements[0];
			//	std::fill(elements.begin(), elements.end(), -1);
			//	elements[0] = firstEnd;
			//}
			if (currentAtom == secondEnd)
				break;
			else if (pos <= length)
			{
				pos = 1;
				currentAtom = elements[0];
				std::fill(elements.begin(), elements.end(), -1);
				elements[0] = firstEnd;
			}
			else if (pos > length)
			{
				std::cerr << "ERROR in mainChain(), maximum length exceeded!" << std::endl;
			}
		}
	}

	return elements;
}

std::vector<int> Molecola::listOfMethyl()
{
	std::vector<int> chain = mainChain();
	return listOfMethyl(chain);
}

std::vector<int> Molecola::listOfEthyl()
{
	std::vector<int> chain = mainChain();
	return listOfEthyl(chain);
}

std::vector<int> Molecola::listOfMethyl(std::vector<int> chain)
{
	std::vector<int> methPos(0);

	for (int i = 0; i < chain.size(); i++)
	{
		if (legami(chain[i]) > 2)
		{
			for (int j = 0; j < legami(chain[i]); j++)
			{
				if (isElementPresent(&chain, struttura[chain[i]][j + 2]))
					continue;
				if (legami(struttura[chain[i]][j + 2]) == 1)
					methPos.push_back(i+1);
			}
		}
	}

	return methPos;
}

std::vector<int> Molecola::listOfEthyl(std::vector<int> chain)
{
	std::vector<int> ethPos(0);
	int asd = 0;
	for (int i = 0; i < chain.size(); i++)
	{
		if (legami(chain[i]) > 2)
		{
			for (int j = 0; j < legami(chain[i]); j++)
			{
				if (isElementPresent(&chain, struttura[chain[i]][j + 2]))
					continue;
				if (legami(struttura[chain[i]][j + 2]) == 2)
				{
					int firstAtom = struttura[chain[i]][j + 2];
					for (int k = 0; k < 2; k++)
					{
						if (struttura[firstAtom][k + 2] == chain[i])
							continue;
						if (legami(struttura[firstAtom][k + 2]) == 1)
							ethPos.push_back(i+1);
					}
				}
			}
		}
	}

	return ethPos;
}


void Molecola::orderedListOfMeAndEt(std::vector<int>* chain, std::vector<int>* mePos, std::vector<int>* etPos, std::vector<int>* posToNamingPos)
{
	chain->clear();
	mePos->clear();
	etPos->clear();
	posToNamingPos->clear();

	*chain = mainChain();
	*mePos = listOfMethyl(*chain);
	*etPos = listOfEthyl(*chain);

	int directSum = 0;
	int inverseSum = 0;
	bool needToBeInverted = false;
	
	if (isStructureSymmetric())
	{
		int numberR = numCrad();
		int numberCOO = numCOORad();
		int numberCOOH = numCOOH();
		int numberEther = numberEthero();
		int numberOLE = numberOle();
		int numberCO = numKeto();
		int numLinEth = numLinEther();
		int numCOR = numCORad();

		if (numLinEth > 0)
		{
			if (dist((*chain)[0], trova(8)) > dist((*chain)[chain->size() - 1], trova(8)))
				needToBeInverted = true;
		}
		else if (numberOLE > 0)
		{
			if (dist((*chain)[0], posOle()) > dist((*chain)[chain->size() - 1], posOle()))
				needToBeInverted = true;
		}
		else if (numberR > 0)
		{
			if (dist((*chain)[0], posCrad()) > dist((*chain)[chain->size() - 1], posCrad()))
				needToBeInverted = true;
		}
		else if (numberEther > 0)
		{
			if (dist((*chain)[0], posEthero()[0]) > dist((*chain)[chain->size() - 1], posEthero()[0]))
				needToBeInverted = true;
		}
		else if (numberCO > 0)
		{
			if (dist((*chain)[0], posKeto()) > dist((*chain)[chain->size() - 1], posKeto()))
				needToBeInverted = true;
		}
		else if (numberCOOH > 0)
		{
			if (kindOfSPecies() == ROOH2_)
			{
				std::vector<int> pos = posOOHinROOH2();
				if ((dist((*chain)[0], pos[0]) + (dist((*chain)[0], pos[1]))) > 
					((dist((*chain)[chain->size() - 1], pos[0])) + (dist((*chain)[chain->size() - 1], pos[1]))))
					needToBeInverted = true;
			}
			else
			{
				if (dist((*chain)[0], posCOOH()) > dist((*chain)[chain->size() - 1], posCOOH()))
					needToBeInverted = true;
			}
		}
		else if (numberCOO > 0)
		{
			if (dist((*chain)[0], posCOOrad()) > dist((*chain)[chain->size() - 1], posCOOrad()))
				needToBeInverted = true;
		}
		else if (numCOR > 0)
		{
			if (dist((*chain)[0], posCOrad()) > dist((*chain)[chain->size() - 1], posCOrad()))
				needToBeInverted = true;
		}
	}
	else 
	{
		if (etPos->size() > 0)
		{
			for (int i = 0; i < etPos->size(); i++)
			{
				directSum += (*etPos)[i];
				inverseSum += chain->size() + 1 - (*etPos)[i];
			}
			if (inverseSum < directSum)
				needToBeInverted = true;
		}
		else if (mePos->size() > 0)
		{
			for (int i = 0; i < mePos->size(); i++)
			{
				directSum += (*mePos)[i];
				inverseSum += chain->size() + 1 - (*mePos)[i];
			}
			if (inverseSum < directSum)
				needToBeInverted = true;
		}
	}

	if (needToBeInverted)
	{
		std::vector<int> invertedChain(chain->size());
		for (int i = 0; i < chain->size(); i++)
		{
			invertedChain[i] = (*chain)[chain->size() - i - 1];
		}
		(*chain) = invertedChain;
		for (int i = 0; i < mePos->size(); i++)
			(*mePos)[i] = chain->size() + 1 - (*mePos)[i];
		for (int i = 0; i < etPos->size(); i++)
			(*etPos)[i] = chain->size() + 1 - (*etPos)[i];
	}

	std::vector<bool> isNumbered(size(), false);

	posToNamingPos->resize(size());
	// map the chain
	for (int i = 0; i < chain->size(); i++)
	{
		(*posToNamingPos)[(*chain)[i]-1] = i + 1;
		isNumbered[(*chain)[i] - 1] = true;
	}
	// map all the other carbons
	int j = 1;
	if (needToBeInverted == false)
	{
		for (int i = 0; i < size(); i++)
		{
			if (isNumbered[i] == false)
			{
				(*posToNamingPos)[i] = chain->size() + j;
				j++;
			}
		}
	}
	else
	{
		for (int i = size()-1; i >= 0; i--)
		{
			if (isNumbered[i] == false)
			{
				(*posToNamingPos)[i] = chain->size() + j;
				j++;
			}
		}
	}

	std::sort((*mePos).begin(), (*mePos).end());
	std::sort((*etPos).begin(), (*etPos).end());

}

Molecola Molecola::noRadicalMolecule()
{
	Molecola mol = *this;

	if (posCrad() == 0 && posCOOrad() == 0 && posCOrad() == 0)
		std::cerr << "WARNING: noRadicalMolecule called on a molecule without radicals!" << std::endl;
	else if (posCrad() != 0 && posCOOrad() == 0 && posCOrad() == 0)
		mol.Crad_to_C(mol.posCrad());
	else if (posCrad() == 0 && posCOOrad() != 0 && posCOrad() == 0)
		mol.COOrad_to_COOH(mol.posCOOrad());
	else if (posCrad() == 0 && posCOOrad() == 0 && posCOrad() != 0)
		mol.COrad_to_CO(mol.posCOrad());
	return mol;
}

Molecola Molecola::parentFuel()
{
	Molecola mol = *this;
	for (int i = 0; i < SIZEMAX; i++)
	{
		if (mol.struttura[i + 1][1] != 0 && mol.struttura[i + 1][1] != 8)
		{
			mol.struttura[i + 1][1] = 1;
		}
		mol.ete[i + 1][0] = 0;
		mol.ete[i + 1][1] = 0;
		mol.ole[i + 1][0] = 0;
		mol.ole[i + 1][1] = 0;
		mol.doppioO[i + 1] = 0;
	}
	mol.numole = 0;
	mol.numete = 0;
	mol.numO = 0;

	return mol;
}

bool Molecola::isSymmetric()
{
	std::vector<int> maiChain = mainChain();

	bool isSymmetric = false;

	if (maiChain.size() % 2 == 0)
	{
		int pos1 = maiChain[maiChain.size() / 2 - 1];
		int pos2 = maiChain[maiChain.size() / 2 + 1 - 1];
		Molecola m1, m2;
		spezza(pos1, pos2, &m1, &m2);
		//std::cout << m1 << std::endl << m2 << std::endl;
		//m1.fix();
		//m2.fix();
		if (m1 == m2)
			isSymmetric = true;
	}
	else
	{
		int pos1 = maiChain[(maiChain.size() - 1) / 2 - 1];
		int pos2 = maiChain[(maiChain.size() + 1) / 2 - 1];
		int pos3 = maiChain[(maiChain.size() + 3) / 2 - 1];
		Molecola m1, m2, m3, m4;
		spezza(pos1, pos2, &m1, &m2);
		spezza(pos3, pos2, &m3, &m4);
		//std::cout << m1 << std::endl << m2 << std::endl << m3 << std::endl << m4 << std::endl;
		if (m1 == m3)
			isSymmetric = true;
	}

	return isSymmetric;
}

bool Molecola::isStructureSymmetric()
{
	Molecola structMol = parentFuel();
	return structMol.isSymmetric();
	//std::vector<int> maiChain = structMol.mainChain();
	//
	//bool isSymmetric = false;
	//
	//if (maiChain.size() % 2 == 0)
	//{
	//	int pos1 = maiChain[maiChain.size() / 2 - 1];
	//	int pos2 = maiChain[maiChain.size() / 2 + 1 - 1];
	//	Molecola m1, m2;
	//	structMol.spezza(pos1, pos2, &m1, &m2);
	//	//std::cout << m1 << std::endl << m2 << std::endl;
	//	//m1.fix();
	//	//m2.fix();
	//	if (m1 == m2)
	//		isSymmetric = true;
	//}
	//else 
	//{
	//	int pos1 = maiChain[(maiChain.size()-1) / 2 - 1];
	//	int pos2 = maiChain[(maiChain.size()+1) / 2 - 1];
	//	int pos3 = maiChain[(maiChain.size()+3) / 2 - 1];
	//	Molecola m1, m2, m3, m4;
	//	structMol.spezza(pos1, pos2, &m1, &m2);
	//	structMol.spezza(pos3, pos2, &m3, &m4);
	//	//std::cout << m1 << std::endl << m2 << std::endl << m3 << std::endl << m4 << std::endl;
	//	if (m1 == m3)
	//		isSymmetric = true;
	//}
	//
	//return isSymmetric;
}

bool Molecola::isChiral()
{
	bool isChiral = false;
	for (int i = 1; i < size() + 1; i++)
	{
		//std::cout << "ATOM NUM " << i << std::endl;
		if (ischeto(i) || tipo(i) == 2 || isole(i) || numAbstractableH(i) > 1)
			continue;
		
		Molecola m1, m2;
		std::vector<Molecola> molVec;
		std::vector<int> atomType;
		for (int j = 2; j < 6; j++)
		{
			if (struttura[i][j] == 0)
				continue;
			spezza(i, struttura[i][j], &m1, &m2);
			//std::cout << std::setw(20) << m1 << std::setw(20) << m2 << std::endl << std::endl;
			molVec.push_back(m2);
			atomType.push_back(tipo(struttura[i][j]));
		}

		bool thereAreEqualGroups = false;
		for (int j = 0; j < molVec.size(); j++)
		{
			for (int k = j + 1; k < molVec.size(); k++)
			{
				if (molVec[j] == molVec[k] && atomType[j] == atomType[k])
					thereAreEqualGroups = true;
			}
		}

		if (thereAreEqualGroups == false)
		{
			isChiral = true;
			break;
		}
		//std::cout << std::endl;
	}
	return isChiral;
}



int Molecola::numberOfRotors()
{
	std::vector<std::vector<int>> bondsList = listOfCCBonds();

	int numRot = 0;

	// check number of rotors
	for (int i = 0; i < bondsList.size(); i++)
	{
		int b1 = bondsList[i][0];
		int b2 = bondsList[i][1];
		if(isBondRotor(b1, b2))
			numRot++;
	}

	// add rotors of the functional groups
	for (int i = 1; i < size() + 1; i++)
	{
		switch (tipo(i))
		{
		case 3:	//-COO*
			numRot += 2;
			break;
		case 4:	//-COOH
			numRot += 2;
			break;
		}
	}

	return numRot;
}

bool Molecola::isBondRotor(int b1, int b2)
{
	if (areBonded(b1, b2) == false)
	{
		std::cerr << "ERROR: isBondRotor called on non bonded atoms!" << std::endl;
		return false;
	}
	if (isole(b1) && isole(b2))	// if bond is a double bond it cannot rotate
		return false;

	if ((tipo(b1) == 2 && isole(b2) && isole(b1) == false) ||	// if the bond is a delocalized double bond - radical it cannot rotate
		(tipo(b2) == 2 && isole(b1) && isole(b2) == false))
		return false;

	if (isInEthRing(b1) && isInEthRing(b2)) // if it is in a cyclic ether ring it cannot rotate
		return false;

	return true;
}

bool isBondInVec(int p1, int p2, std::vector<std::vector<int>>* vec)
{
	for (int i = 0; i < vec->size(); i++)
		if ((p1 == (*vec)[i][0] && p2 == (*vec)[i][1]) || (p1 == (*vec)[i][1] && p2 == (*vec)[i][0]))
			return true;
	return false;
}

int Molecola::numberOfSymmetries()
{
	int numSym = 1;
	std::vector<std::vector<int>> rotatedBonds;

	for (int i = 1; i < size() + 1; i++)
	{
		bool alreadyRotated = false;
		if (numAbstractableH(i) == 3)		// CH3
			for (int j = 2; j < 6; j++)
				if (struttura[i][j] != 0)
					if (isBondRotor(i, struttura[i][j]) && isBondInVec(i, struttura[i][j], &rotatedBonds) == false)
					{
						rotatedBonds.push_back(std::vector<int> {i, struttura[i][j]});
						numSym *= 3;
						alreadyRotated = true;
						break;
					}

		if (alreadyRotated)
			continue;
		
		int numBonds = 0;
		for (int j = 2; j < 6; j++)
			if (struttura[i][j] != 0)
				numBonds++;
		if (numBonds == 4)			// check if there are 3 equal groups
		{
			Molecola m1, m2, m3, m4, m5;
			int t1, t2, t3, t4, t5;
			spezza(i, struttura[i][2], &m1, &m2);
			t2 = tipo(struttura[i][2]);
			spezza(i, struttura[i][3], &m1, &m3);
			t3 = tipo(struttura[i][3]);
			spezza(i, struttura[i][4], &m1, &m4);
			t4 = tipo(struttura[i][4]);
			spezza(i, struttura[i][5], &m1, &m5);
			t5 = tipo(struttura[i][5]);
		
			if (m2 == m3 && m3 == m4 && t2 == t3 && t3 == t4)
				if (isBondRotor(i, struttura[i][5]) && isBondInVec(i, struttura[i][5], &rotatedBonds) == false)
				{
					numSym *= 3;
					rotatedBonds.push_back(std::vector<int> {i, struttura[i][5]});
				}
			
			if (m3 == m4 && m4 == m5 && t3 == t4 && t4 == t5)
				if (isBondRotor(i, struttura[i][2]) && isBondInVec(i, struttura[i][2], &rotatedBonds) == false)
				{
					numSym *= 3;
					rotatedBonds.push_back(std::vector<int> {i, struttura[i][2]});
				}
			if (m2 == m4 && m4 == m5 && t2 == t4 && t4 == t5)
				if (isBondRotor(i, struttura[i][3]) && isBondInVec(i, struttura[i][3], &rotatedBonds) == false)
				{
					numSym *= 3;
					rotatedBonds.push_back(std::vector<int> {i, struttura[i][3]});
				}

			if (m2 == m3 && m3 == m5 && t2 == t3 && t3 == t5)
				if (isBondRotor(i, struttura[i][4]) && isBondInVec(i, struttura[i][4], &rotatedBonds) == false)
				{
					numSym *= 3;
					rotatedBonds.push_back(std::vector<int> {i, struttura[i][4]});
				}
		}
	}
	
	if (isSymmetric())
		numSym *= 2;
	
	return numSym;
}

bool Molecola::isCradAllylic(int pos)
{
	if (tipo(pos) != 2)
	{
		std::cerr << "WARNING: isCradAllylic called on a non Crad atom!" << std::endl;
		return false;
	}

	if (isole(pos))
		return false;
	
	int found[SIZEMAX + 1];
	scorri(pos, 1, found);
	for (int i = 1; i < SIZEMAX + 1; i++)
		if (found[i] == 1)
			if (isole(i))
				return true;
	
	return false;
}

bool Molecola::bothC4H713( Molecola mol2)
{
	alreadyCheckingC4H713 = true;
	mol2.alreadyCheckingC4H713 = true;
	Molecola c4h72_1;
	c4h72_1.makeC4H72_1();
	c4h72_1.alreadyCheckingC4H713 = true;
	Molecola c4h71_3;
	c4h71_3.makeC4H71_3();
	c4h71_3.alreadyCheckingC4H713 = true;
	
	if ((*this == c4h72_1 && mol2 == c4h71_3) || (mol2 == c4h72_1 && *this == c4h71_3))
	{
		alreadyCheckingC4H713 = false;
		mol2.alreadyCheckingC4H713 = false;
		return true;
	}

	mol2.alreadyCheckingC4H713 = false;
	alreadyCheckingC4H713 = false;
	return false;
}


std::string Molecola::inchiName()
{
	INCHI_API::inchi_Input molInchi;

	molInchi.num_atoms = numberOfC() + numberOfH() + numberOfO();
	
	INCHI_API::inchi_Atom atoms[1000];

	int indExtraAtom = 0;

	for (int i = 0; i < size(); i++)
	{
		if (tipo(i + 1) != 8)
		{
			INCHI_API::inchi_Atom at;
			char atName[ATOM_EL_LEN] = "C";
			*(at.elname) = *atName;

			at.charge = 0;

			if (tipo(i + 1) == 2)
				at.radical = INCHI_RADICAL_SINGLET;
			else
				at.radical = INCHI_RADICAL_NONE;

			int numBonds = 0;
			// add bonded carbons (or linear cyclic ether)
			for (int j = 0; j < numGraphConnections(i + 1); j++)
			{
				at.neighbor[j] = struttura[i + 1][j + 2] - 1;
				at.bond_stereo[j] = INCHI_BOND_STEREO_NONE;

				if (isole(i + 1, struttura[i + 1][j + 2]))
					at.bond_type[j] = INCHI_BOND_TYPE_DOUBLE;
				else
					at.bond_type[j] = INCHI_BOND_TYPE_SINGLE;
				numBonds++;
			}

			// add bonded hydrogens
			for (int j = 0; j < numAbstractableH(i + 1); j++)
			{
				int IDHydrogen = size() + indExtraAtom;
				at.neighbor[numBonds] = IDHydrogen;
				at.bond_type[numBonds] = INCHI_BOND_TYPE_SINGLE;
				at.bond_stereo[numBonds] = INCHI_BOND_STEREO_NONE;
				numBonds++;
				indExtraAtom++;

				INCHI_API::inchi_Atom hyd;
				char hydName[ATOM_EL_LEN] = "H";
				*(hyd.elname) = *hydName;
				hyd.neighbor[0] = i;
				hyd.bond_type[0] = INCHI_BOND_TYPE_SINGLE;
				hyd.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
				hyd.charge = 0;
				hyd.radical = INCHI_RADICAL_NONE;
				hyd.num_iso_H[0] = -1;
				hyd.num_bonds = 1;
				hyd.isotopic_mass = 0;
				atoms[IDHydrogen] = hyd;
				//std::cout << "ID hyd = " << IDHydrogen << ", bonded to " << hyd.neighbor[0] << std::endl;
			}

			// add COO*
			if (tipo(i + 1) == 3)
			{
				int IDfirstO = size() + indExtraAtom;
				at.neighbor[numBonds] = IDfirstO;
				at.bond_type[numBonds] = INCHI_BOND_TYPE_SINGLE;
				at.bond_stereo[numBonds] = INCHI_BOND_STEREO_NONE;
				numBonds++;
				indExtraAtom++;

				// add first O
				INCHI_API::inchi_Atom firstO;
				char firstOName[ATOM_EL_LEN] = "O";
				*(firstO.elname) = *firstOName;
				int IDsecondO = size() + indExtraAtom;
				indExtraAtom++;
				firstO.neighbor[0] = i;
				firstO.bond_type[0] = INCHI_BOND_TYPE_SINGLE;
				firstO.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
				firstO.neighbor[1] = IDsecondO;
				firstO.bond_type[1] = INCHI_BOND_TYPE_SINGLE;
				firstO.bond_stereo[1] = INCHI_BOND_STEREO_NONE;
				firstO.charge = 0;
				firstO.radical = INCHI_RADICAL_NONE;
				firstO.num_iso_H[0] = -1;
				firstO.num_bonds = 2;
				firstO.isotopic_mass = 0;
				atoms[IDfirstO] = firstO;

				// add second O2
				INCHI_API::inchi_Atom secondO;
				char secondOName[ATOM_EL_LEN] = "O";
				*(secondO.elname) = *secondOName;
				secondO.neighbor[0] = IDfirstO;
				secondO.bond_type[0] = INCHI_BOND_TYPE_SINGLE;
				secondO.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
				secondO.charge = 0;
				secondO.radical = INCHI_RADICAL_SINGLET;
				secondO.num_iso_H[0] = -1;
				secondO.num_bonds = 1;
				secondO.isotopic_mass = 0;
				atoms[IDsecondO] = secondO;
			}

			// add COOH
			if (tipo(i + 1) == 4)
			{
				int IDfirstO = size() + indExtraAtom;
				at.neighbor[numBonds] = IDfirstO;
				at.bond_type[numBonds] = INCHI_BOND_TYPE_SINGLE;
				at.bond_stereo[numBonds] = INCHI_BOND_STEREO_NONE;
				numBonds++;
				indExtraAtom++;

				// add first O
				INCHI_API::inchi_Atom firstO;
				char firstOName[ATOM_EL_LEN] = "O";
				*(firstO.elname) = *firstOName;
				int IDsecondO = size() + indExtraAtom;
				indExtraAtom++;
				firstO.neighbor[0] = i;
				firstO.bond_type[0] = INCHI_BOND_TYPE_SINGLE;
				firstO.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
				firstO.neighbor[1] = IDsecondO;
				firstO.bond_type[1] = INCHI_BOND_TYPE_SINGLE;
				firstO.bond_stereo[1] = INCHI_BOND_STEREO_NONE;
				firstO.charge = 0;
				firstO.radical = INCHI_RADICAL_NONE;
				firstO.num_iso_H[0] = -1;
				firstO.num_bonds = 2;
				firstO.isotopic_mass = 0;
				atoms[IDfirstO] = firstO;

				// add second O
				INCHI_API::inchi_Atom secondO;
				char secondOName[ATOM_EL_LEN] = "O";
				*(secondO.elname) = *secondOName;
				int IDhydr = size() + indExtraAtom;
				indExtraAtom++;
				secondO.neighbor[0] = IDfirstO;
				secondO.bond_type[0] = INCHI_BOND_TYPE_SINGLE;
				secondO.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
				secondO.neighbor[1] = IDhydr;
				secondO.bond_type[1] = INCHI_BOND_TYPE_SINGLE;
				secondO.bond_stereo[1] = INCHI_BOND_STEREO_NONE;
				secondO.charge = 0;
				secondO.radical = INCHI_RADICAL_NONE;
				secondO.num_iso_H[0] = -1;
				secondO.num_bonds = 2;
				secondO.isotopic_mass = 0;
				atoms[IDsecondO] = secondO;

				// add hydrogen
				INCHI_API::inchi_Atom hyd;
				char hydName[ATOM_EL_LEN] = "H";
				*(hyd.elname) = *hydName;
				hyd.neighbor[0] = IDsecondO;
				hyd.bond_type[0] = INCHI_BOND_TYPE_SINGLE;
				hyd.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
				hyd.charge = 0;
				hyd.radical = INCHI_RADICAL_NONE;
				hyd.num_iso_H[0] = -1;
				hyd.num_bonds = 1;
				hyd.isotopic_mass = 0;
				atoms[IDhydr] = hyd;
			}

			// add C=O
			if (ischeto(i + 1))
			{
				int IDOxy = size() + indExtraAtom;
				at.neighbor[numBonds] = IDOxy;
				at.bond_type[numBonds] = INCHI_BOND_TYPE_DOUBLE;
				at.bond_stereo[numBonds] = INCHI_BOND_STEREO_NONE;
				numBonds++;
				indExtraAtom++;

				// add  O
				INCHI_API::inchi_Atom oxy;
				char oxyName[ATOM_EL_LEN] = "O";
				*(oxy.elname) = *oxyName;
				oxy.neighbor[0] = i;
				oxy.bond_type[0] = INCHI_BOND_TYPE_DOUBLE;
				oxy.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
				oxy.charge = 0;
				oxy.radical = INCHI_RADICAL_NONE;
				oxy.num_iso_H[0] = -1;
				oxy.num_bonds = 1;
				oxy.isotopic_mass = 0;
				atoms[IDOxy] = oxy;
			}

			// add C=O*
			if (tipo(i + 1) == 9)
			{
				int IDOxy = size() + indExtraAtom;
				at.neighbor[numBonds] = IDOxy;
				at.bond_type[numBonds] = INCHI_BOND_TYPE_SINGLE;
				at.bond_stereo[numBonds] = INCHI_BOND_STEREO_NONE;
				numBonds++;
				indExtraAtom++;

				// add  O
				INCHI_API::inchi_Atom oxy;
				char oxyName[ATOM_EL_LEN] = "O";
				*(oxy.elname) = *oxyName;
				oxy.neighbor[0] = i;
				oxy.bond_type[0] = INCHI_BOND_TYPE_SINGLE;
				oxy.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
				oxy.charge = 0;
				oxy.radical = INCHI_RADICAL_SINGLET;
				oxy.num_iso_H[0] = -1;
				oxy.num_bonds = 1;
				oxy.isotopic_mass = 0;
				atoms[IDOxy] = oxy;
			}
			//std::cout << "Carbon " << i << ", num bonds: " << numBonds << std::endl;
			at.num_bonds = numBonds;
			at.num_iso_H[0] = -1;
			at.isotopic_mass = 0;
			atoms[i] = at;
		}
		else      // O of linear ether
		{
			INCHI_API::inchi_Atom at;
			char atName[ATOM_EL_LEN] = "O";
			*(at.elname) = *atName;
			at.charge = 0;
			at.radical = INCHI_RADICAL_NONE;
			at.neighbor[0] = struttura[i + 1][2] - 1;
			at.bond_type[0] = INCHI_BOND_TYPE_SINGLE;
			at.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
			at.neighbor[1] = struttura[i + 1][3] - 1;
			at.bond_type[1] = INCHI_BOND_TYPE_SINGLE;
			at.bond_stereo[1] = INCHI_BOND_STEREO_NONE;
			at.num_bonds = 2;
			at.num_iso_H[0] = -1;
			at.isotopic_mass = 0;
			atoms[i] = at;
		}
	}

	// add cyclic ethers
	for (int i = 0; i < numete; i++)
	{
		// update first carbon
		int indFirstC = ete[i + 1][0] - 1;
		int IDOxy = size() + indExtraAtom;
		atoms[indFirstC].neighbor[atoms[indFirstC].num_bonds] = IDOxy;
		atoms[indFirstC].bond_type[atoms[indFirstC].num_bonds] = INCHI_BOND_TYPE_SINGLE;
		atoms[indFirstC].bond_stereo[atoms[indFirstC].num_bonds] = INCHI_BOND_STEREO_NONE;
		atoms[indFirstC].num_bonds++;
		indExtraAtom++;

		// update second carbon
		int indSecondC = ete[i + 1][1] - 1;
		atoms[indSecondC].neighbor[atoms[indSecondC].num_bonds] = IDOxy;
		atoms[indSecondC].bond_type[atoms[indSecondC].num_bonds] = INCHI_BOND_TYPE_SINGLE;
		atoms[indSecondC].bond_stereo[atoms[indSecondC].num_bonds] = INCHI_BOND_STEREO_NONE;
		atoms[indSecondC].num_bonds++;

		// add oxygen
		INCHI_API::inchi_Atom oxy;
		char oxyName[ATOM_EL_LEN] = "O";
		*(oxy.elname) = *oxyName;
		oxy.neighbor[0] = indFirstC;
		oxy.bond_type[0] = INCHI_BOND_TYPE_SINGLE;
		oxy.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
		oxy.neighbor[1] = indSecondC;
		oxy.bond_type[0] = INCHI_BOND_TYPE_SINGLE;
		oxy.bond_stereo[0] = INCHI_BOND_STEREO_NONE;
		oxy.charge = 0;
		oxy.radical = INCHI_RADICAL_SINGLET;
		oxy.num_iso_H[0] = -1;
		oxy.num_bonds = 2;
		oxy.isotopic_mass = 0;
		atoms[IDOxy] = oxy;
	}

	molInchi.atom = atoms;
	INCHI_API::inchi_Output out;
	int flag = INCHI_API::GetINCHI(&molInchi, &out);
	std::string outStr(out.szInChI);
	return outStr;
}

std::vector<int> Molecola::listOfBondedC(int pos)
{
	std::vector<int> list;
	for (int i = 1; i < SIZEMAX + 1; i++)
	{
		if (i == pos)
			continue;
		if (areBonded(pos, i))
			list.push_back(i);
	}
	return list;
}