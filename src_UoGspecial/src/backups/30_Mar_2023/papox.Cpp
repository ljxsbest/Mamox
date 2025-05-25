//**********************************************************************
//**********************************************************************
//**                                                                  **
//**   ######                                                         **
//**   #     #    ##    #####    ####   #    #                        **
//**   #     #   #  #   #    #  #    #   #  #                         **
//**   ######   #    #  #    #  #    #    ##                          **
//**   #        ######  #####   #    #    ##                          **
//**   #        #    #  #       #    #   #  #                         **
//**   #        #    #  #        ####   #    #                        **
//**                                                                  **
//**                                                                  **
//**  Last Release: 15/4/96                                           **
//**  Author:       Enrico Garavaglia,                                **
//**                Via Monti Sabini, 24                              **
//**                20141 Milano Italy                                **
//**                Email: mr..vain@bbs.infosquare.it                 **
//**  Scope:        This program is used to generate and solve        **
//**                an oxidation scheme for whichever ramificated     **
//**                paraffins; product ripartition for nearly null    **
//**                times, the entire detailed scheme used and the    **
//**                branching velocity are obtained after execution.  **
//**                Full help on using the program please             **
//**                   refer to the file PAPOX.TXT or                 **
//**                   type papox.exe -h                              **
//**  Distribution: This work can be used freely (no cost!) and       **
//**                distributed, as long as these file are all        **
//**                present and unmodified:                           **
//**                   papox.cpp                                      **
//**                   molec.hpp                                      **
//**                   molec.cpp                                      **
//**                   kinox.hpp                                      **
//**                   kinox.cpp                                      **
//**                   papox.txt                                      **
//**                                                                  **
//**                                                                  **
//**                                                                  **
//**                                                                  **
//**                                                                  **
//**                                                                  **
//**********************************************************************
//**********************************************************************


#include <stdlib.h>
#include "kinox.h"
#include "molec.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define NOMINMAX		// need to be defined before inclusion of windows.h because otherwise there is a problem with the functions std::min and std::max
#include <windows.h>

#include "chemkinOut.h"
#include "reaction.h"
#include "thermoOut.h"
#include <algorithm>
#include <thread>
#include "chemkinParser.h"
#include "Utilities.h"
#include <chrono>
#include <filesystem>


#include "definitions.h"


struct Pos {
	int ooh;
	int rad;
};

//==============================================================
//-------------------Global Variables---------------------------
//==============================================================

int tipor[SIZEMAX + 1];
int tipoqooh[SIZEMAX + 1][SIZEMAX + 1];

std::vector<double> Temps{ TEMPERATURES };
std::vector<double> Pressures{ PRESSURES };

int Molecola::numOfTemps = Temps.size();
int OH_atoms[3] = { 0, 1, 1 };
////////////////////////////////////////////////////////////////
//  queste matrici coprono come dimensioni tutte le possibili
//  posizioni dei radicali e dei gruppi OOH; il loro contenuto
//  indica a quale isomero esse si riferiscano
//  tipor[i] va da 1 a numr
//  tipoqooh[i][j] va da 1 a numqooh
////////////////////////////////////////////////////////////////

ofstream reaout;    // reactions output file
ofstream lumpout;   // lumped reactions output file
ofstream distrib;	// lumped products distributions

// start the object for outputting the chemkin files
ChemkinOut chemout("OUTPUT/reaction_detailed.inp");

//==============================================================
//----------------------Functions-------------------------------
//==============================================================

int pos_r(int numero);
Pos pos_qooh(int numero);
void print1(int numr, int numqooh, Molecola r[], Molecola roo[],
	Molecola qooh[], Molecola ooqooh[]);
//void print2(int numr, int numqooh, Molecola* r, Molecola* roo,
//	Molecola* qooh, Molecola* ooqooh);
void err(char s[]);
int Add(Molecola nuova, Molecola old[], int* num, int max);
void decomponi(Molecola HC, Kinox& k,
	Molecola* r, Molecola* roo, Molecola* qooh, Molecola* ooqooh,
	int* numr, int* numqooh);
//double makeKini1(Molecola& HC, Kinox& kin);
//double makeKini2(Molecola& HC, Kinox& kin);
double makeKprod(Molecola& HC, Kinox& kin, Molecola* r, int numr);
void help(void);

//--------------------------------------------------------------

int pos_r(int numero)
////////////////////////////////////////////////////////////////
//  Returns the position of the radical in the chain.
//  Parameters:
//      numero: number of the radical (between 1 and numr)
////////////////////////////////////////////////////////////////
{
	// manca controllo che numero sia da 1 a numr
	int pos = 1;
	while (tipor[pos] != numero) pos++;
	return pos;
};

Pos pos_qooh(int numero)
////////////////////////////////////////////////////////////////
//  Return the position of the COOH anc C* in the chain
//  Parameters:
//      numero: number of the COOH group
////////////////////////////////////////////////////////////////
{
	// manca controllo che numero sia da 1 a numqooh
	Pos pos;
	int ooh, r;
	for (ooh = 1; ooh <= SIZEMAX; ooh++)
		for (r = 1; r <= SIZEMAX; r++)
			if (tipoqooh[ooh][r] == numero)
			{
				pos.ooh = ooh;
				pos.rad = r;
				return pos;
			};
	pos.ooh = pos.rad = 0;
	cerr << "pos_qooh() non corretta";  // @@ linea di debug
	return pos;
};

void print1(int numr, int numqooh, Molecola r[], Molecola roo[],
	Molecola qooh[], Molecola ooqooh[])
	////////////////////////////////////////////////////////////////
	//  Generated radicals output               ITA: Output dei radicali generati
	////////////////////////////////////////////////////////////////
{
	std::cout << "Generated radicals:" << endl;
	std::cout << "----------------------------------------------" << endl;
	std::cout << "R*          " << numr << endl;
	std::cout << "ROO*        " << numr << endl;
	std::cout << "*QOOH       " << numqooh << endl;
	std::cout << "*OOQOOH     " << numqooh << endl << endl;

	char buff[80];
	int i;
	std::cout << "R*:" << endl;
	std::cout << "^^^" << endl;
	for (i = 1; i <= numr; i++)
	{
		sprintf(buff, "%2i) r(%2i)         ", i, tipor[i]);
		std::cout << buff << r[i] << endl;
	};
	std::cout << endl;
	std::cout << "ROO*:" << endl;
	std::cout << "^^^^^" << endl;
	for (i = 1; i <= numr; i++)
	{
		sprintf(buff, "%2i) roo(%2i)       ", i, tipor[i]);
		std::cout << buff << roo[i] << endl;
	};
	std::cout << endl;
	std::cout << "*QOOH:" << endl;
	std::cout << "^^^^^^" << endl;
	for (i = 1; i <= numqooh; i++)
	{
		Pos p = pos_qooh(i);
		sprintf(buff, "%2i) qooh(%2i,%2i)     ", i, p.ooh, p.rad);
		std::cout << buff << qooh[i] << endl;
	};
	std::cout << endl;
	std::cout << "*OOQOOH:" << endl;
	std::cout << "^^^^^^^^" << endl;
	for (i = 1; i <= numqooh; i++)
	{
		Pos p = pos_qooh(i);
		sprintf(buff, "%2i) ooqooh(%2i,%2i)   ", i, p.ooh, p.rad);
		std::cout << buff << ooqooh[i] << endl;
	};
	std::cout << endl;
};


void err(char s[])
////////////////////////////////////////////////////////////////
//  Print an error message and exit the code
//  Parameters:
//      s[]: the error message
////////////////////////////////////////////////////////////////
{
	cerr << s << endl << endl;
	exit(100);
};

int Add(Molecola nuova, Molecola old[], int* num, int max)
////////////////////////////////////////////////////////////////
//  Add the moleule "new" to the vector "old" if									Aggiunge la molecola new al vettore old se
//  it is not already present, otherways do not add it.								non esisteva ancora, altrimenti no
//  Returns the position of the molecule in the vector, and evetually				rende la posizione della molecola nel vettore ed eventualmente
//  update the dimension "num".														aggiorna la sua dimensione: num
////////////////////////////////////////////////////////////////
{
	int i;
	for (i = 1; i <= *num; i++)
		if (nuova == old[i]) break;
	if (i > *num)
	{
		if (i > max) {
			cerr << "Vettore insuff. " << nuova << endl;
			exit(300);
		}
		else old[++ * num] = nuova;
	};
	return i;
};

int Add(Molecola newMol, std::vector<Molecola>* vec)   // add newMol to vec
{
	for (int i = 0; i < vec->size(); i++)
		if (newMol == (*vec)[i]) return i;	// if it's present return it's index

	vec->push_back(newMol);					// else add it and return it's index
	return vec->size() - 1;
};

int Add(Reaction newReac, std::vector<Reaction>* vec)   // add newMol to vec
{
	for (int i = 0; i < vec->size(); i++)
		if (newReac == (*vec)[i]) return i;	// if it's present return it's index

	vec->push_back(newReac);					// else add it and return it's index
	return vec->size() - 1;
};

int Add(std::string newString, std::vector<std::string>* vec)   // add newMol to vec
{
	for (int i = 0; i < vec->size(); i++)
		if (newString == (*vec)[i]) return i;	// if it's present return it's index

	vec->push_back(newString);					// else add it and return it's index
	return vec->size() - 1;
};

int AddIsomer(Molecola nuova, Molecola old[], int* num, int max)
////////////////////////////////////////////////////////////////
//  Add the moleule "new" to the vector "old" if									
//  none of its isomers are present, otherways do not add it.		
//  Returns the position of the molecule in the vector, and evetually
//  update the dimension "num".										
////////////////////////////////////////////////////////////////
{
	int i;
	for (i = 1; i <= *num; i++)
		if (nuova.isIsomer(&old[i])) break;
	if (i > *num)
	{
		if (i > max) {
			cerr << "Vettore insuff. " << nuova << endl;
			exit(300);
		}
		else old[++ * num] = nuova;
	};
	return i;
};


void printv(std::vector<float> vect, std::string tag)
{
	std::cout << tag << ": ";
	for (int i = 0; i < vect.size(); i++) std::cout << vect[i] << "  ";
	std::cout << std::endl;
}

void printv(std::vector<double> vect, std::string tag)
{
	std::cout << tag << ": ";
	for (int i = 0; i < vect.size(); i++) std::cout << vect[i] << "  ";
	std::cout << std::endl;
}

void coutDebugFlag()
{
	std::cout << "-DEBUG--DEBUG--DEBUG--DEBUG--DEBUG--DEBUG--DEBUG--DEBUG--DEBUG--DEBUG--DEBUG--DEBUG--DEBUG--DEBUG-" << std::endl;
}

std::string kindToString(species spec)
{
	std::string speClass;
	switch (spec)
	{
	case fuel_:
		speClass = "fuel";
		break;
	case R_:
		speClass = "R";
		break;
	case ROO_:
		speClass = "ROO";
		break;
	case QOOH_:
		speClass = "QOOH";
		break;
	case OOQOOH_:
		speClass = "OOQOOH";
		break;
	case OLE_:
		speClass = "OLE";
		break;
	case CO_:
		speClass = "CO";
		break;
	case cEth_:
		speClass = "cEth";
		break;
	case RO_:
		speClass = "RO";
		break;
	case KHP_:
		speClass = "KHP";
		break;
	case ROOH_:
		speClass = "ROO";
		break;
	case POOH2_:
		speClass = "POOH2";
		break;
	case ROOH2_:
		speClass = "ROOH2";
		break;
	case oleOOH_:
		speClass = "oleOOH";
		break;
	case cEthOOH_:
		speClass = "cEthOOH";
		break;
	case oleR_:
		speClass = "oleR";
		break;
	case oleCO_:
		speClass = "oleCO";
		break;
	case cEthR_:
		speClass = "cEthR";
		break;
	case cEthCO_:
		speClass = "cEthCO";
		break;
	case lEthRO_:
		speClass = "lEthRO";
		break;
	case alkRO_:
		speClass = "alkRO";
		break;
	case special_:
		speClass = "special";
		break;
	case unidentified_:
		speClass = "unidentified";
		break;
	default:
		speClass = "ERROR";
		break;
	}
	return speClass;
}

void printConcHeader(std::string title, int sm)
{
	std::cout << title << endl;
	for (int i = 0; i < title.length(); i++) std::cout << "^";
	std::cout << endl;
	std::cout << std::setw(sm) << "Temperatures:" << "    ";
	std::cout << std::setprecision(6);
	for (int z = 0; z < Temps.size(); z++) std::cout << std::setw(10) << Temps[z];
	std::cout << std::endl;
	std::cout << std::setprecision(3);
}

void addReactionToList(std::vector<Reaction>* vect, Molecola reac, Molecola prod, double kinPar[3], std::string label, std::string comment)
{
	Reaction re(reac, prod, kinPar, label, comment);
	vect->push_back(re);
}

void addReactionToList(std::vector<Reaction>* vect, Molecola reac, std::vector<Molecola> prod, double kinPar[3], std::string label, std::string comment)
{
	Reaction re(reac, prod, kinPar, label, comment);
	vect->push_back(re);
}

void addReactionToList(std::vector<Reaction>* vect, std::vector<Molecola> reac, std::vector<Molecola> prod, double kinPar[3], std::string label, std::string comment)
{
	Reaction re(reac, prod, kinPar, label, comment);
	vect->push_back(re);
}

void addReactionToList(std::vector<Reaction>* vect, std::vector<Molecola> reac, Molecola prod, double kinPar[3], std::string label, std::string comment)
{
	Reaction re(reac, prod, kinPar, label, comment);
	vect->push_back(re);
}

std::string oh_abstraction_correction(Molecola *HC, int posR)
{
	std::vector<int> listBonded = HC->listOfBondedC(posR);
	
	std::string correction = "_none";
	switch (HC->tipoC(posR))
	{
	case Cp:
	{


		switch (HC->tipoC(listBonded[0]))
		{
		case Cs:
			correction = "_1";
			break;
		case Ct:
			correction = "_2";
			break;
		case Cq:
			correction = "_3";
			break;
		}
	}
	break;
	case Cs:
	{
		Carbonio C1 = HC->tipoC(listBonded[0]);
		Carbonio C2 = HC->tipoC(listBonded[1]);
		if ((C1 == Cp && C2 == Cs) || (C1 == Cs && C2 == Cp))
			correction = "_01";
		if ((C1 == Cp && C2 == Ct) || (C1 == Ct && C2 == Cp))
			correction = "_02";
		if ((C1 == Cp && C2 == Cq) || (C1 == Cq && C2 == Cp))
			correction = "_03";
		if ((C1 == Cs && C2 == Cs))
			correction = "_11";
		if ((C1 == Cs && C2 == Ct) || (C1 == Ct && C2 == Cs))
			correction = "_12";
		if ((C1 == Ct && C2 == Cq) || (C1 == Cq && C2 == Ct))
			correction = "_23";
	}
	break;
	case Ct:
	{
		Carbonio C1 = HC->tipoC(listBonded[0]);
		Carbonio C2 = HC->tipoC(listBonded[1]);
		Carbonio C3 = HC->tipoC(listBonded[2]);
		if (C1 == Cp && C2 == Cp && C3 == Cp)
			correction = "_000";
		if ((C1 == Cs && C2 == Cp && C3 == Cp) || (C1 == Cp && C2 == Cs && C3 == Cp) || (C1 == Cp && C2 == Cp && C3 == Cs))
			correction = "_001";
		if ((C1 == Ct && C2 == Cp && C3 == Cp) || (C1 == Cp && C2 == Ct && C3 == Cp) || (C1 == Cp && C2 == Cp && C3 == Ct))
			correction = "_002";
		if ((C1 == Cp && C2 == Cs && C3 == Cs) || (C1 == Cs && C2 == Cp && C3 == Cs) || (C1 == Cs && C2 == Cs && C3 == Cp))
			correction = "_011";
	}
	break;
	default:
		break;
	}
	return correction;
}

void h_abstraction(Kinox k, int numr, Molecola HC, Molecola r[], HAbsRad Hrad, std::vector<Reaction>* reacVec, Molecola HabsMol, Molecola HabsMolProd)
{
	std::string nameHrad = k.nameHAbsRad(Hrad);
	std::string label = "H abstraction by " + nameHrad + ":";

	for (int i = 1; i <= numr; i++) {
		int posR = r[i].trova(2);

		std::string correction = "none";

		if (Hrad == oh)
		{
			correction = oh_abstraction_correction(&HC, posR);
		}
		if (Hrad == ch3oo && HC.tipoC(posR) == Cp)
		{
			correction = oh_abstraction_correction(&HC, posR);
		}

		std::string reacomm = k.v_h_abstraction(Hrad, HC.tipoC(posR), HC.numAbstractableH(posR), r[i].isomeri, correction);
		addReactionToList(reacVec, std::vector<Molecola>{ HC, HabsMol }, std::vector<Molecola>{ r[i], HabsMolProd }, new double[3] { k.A, k.n, k.E }, label, reacomm);
		k.wrirea(reaout, 14 + int(Hrad), k.A, k.n, k.E, HC, r[i]);
	}
}

int addUniqueSpecies(std::vector<Molecola*>* vect, Molecola* mol)
{
	for (int i = 0; i < vect->size(); i++)
	{
		if (*mol == *((*vect)[i]))
		{
			return i;
		}
	}
	vect->push_back(mol);
	return -1;
}

int carbonioToInt(Carbonio c)
{
	//enum Carbonio { Cp, Cs, Ct, Cq };
	switch (c)
	{
	case Cp:
		return 1;
		break;
	case Cs:
		return 2;
		break;
	case Ct:
		return 3;
		break;
	case Cq:
		return 4;
		break;
	default:
		std::cerr << "ERROR: carbonioToInt called on a non recognized Carbonio: " << c << ". Aborted." << std::endl;
		exit(0);
		break;
	}
}

std::vector<Molecola> decomposeRO(Molecola m2)  // decompose the species m2 and return the products
{
	Molecola m3, m4;
	std::vector<Molecola> vec;
	// find first atom of the beta scission
	int found[SIZEMAX + 1];
	int pos2 = m2.posCrad();
	m2.scorri(pos2, 2, found);
	int posBeta1 = 0;
	for (int k = 1; k < SIZEMAX + 1; k++)
	{
		if (found[k] == 1)
		{
			if (m2.dist(m2.posCrad(), k) > 2)
				posBeta1 = k;
		}
	}
	if (posBeta1 != 0) // desired beta scission is possible
	{
		// find second atom of the beta scission
		int found2[SIZEMAX + 1];
		m2.scorri(pos2, 1, found);
		m2.scorri(posBeta1, 1, found2);
		int posBeta2 = 0;
		for (int k = 1; k < SIZEMAX + 1; k++)
			if (found[k] == 1 && found2[k] == 1)
				posBeta2 = k;
		int posRad = m2.posCrad();
		m2.spezza(posBeta1, posBeta2, &m3, &m4);
		vec.push_back(m3);
		vec.push_back(m4);
	}
	else     // desired beta scission is not possible
	{
		m2.scorri(pos2, 1, found);
		int found2[SIZEMAX + 1];
		int posBeta2 = 0;
		for (int k = 0; k < SIZEMAX + 1; k++)
		{
			if (found[k] == 1)
			{
				m2.scorri(k, 1, found2);
				for (int l = 0; l < SIZEMAX + 1; l++)
					if (found2[l] == 1 && l != pos2)
					{
						posBeta1 = k;
						posBeta2 = l;
					}
			}
		}
		if (posBeta1 != 0 && posBeta2 != 0)
		{
			m2.spezza(posBeta1, posBeta2, &m3, &m4);
			vec.push_back(m3);
			vec.push_back(m4);
		}
		else
		{
			vec.push_back(m2);
		}
	}
	//std::cout << m2 << std::endl << m3 << std::endl << m4 << std::endl;

	return vec;
}

std::vector<Molecola> fullyDecomposeRO(Molecola m2) // decompose m2 and its products too (if they are RO) and return the products
{
	std::vector<Molecola> decProds;
	decProds = decomposeRO(m2);
	if (decProds.size() > 1)
	{
		for (int i = 0; i < decProds.size(); i++)
		{
			if (decProds[i].kindOfSPecies() == RO_)
			{
				std::vector<Molecola> secDecProds;
				secDecProds = decomposeRO(decProds[i]);
				decProds.erase(decProds.begin() + i);
				decProds.insert(decProds.end(), secDecProds.begin(), secDecProds.end());
			}
		}
	}
	return decProds;
}

std::vector<Molecola> fullyDecomposeRO(std::vector<Molecola> vec) // decompose all the RO in vec and its products too (if they are RO) and return the products
{
	std::vector<Molecola> decProds;
	for (int j = 0; j < vec.size(); j++)
	{
		Molecola m2 = vec[j];
		if (m2.kindOfSPecies() == RO_)
		{
			decProds = decomposeRO(m2);
			if (decProds.size() > 1)
			{
				for (int i = 0; i < decProds.size(); i++)
				{
					if (decProds[i].kindOfSPecies() == RO_)
					{
						std::vector<Molecola> secDecProds;
						secDecProds = decomposeRO(decProds[i]);
						decProds.erase(decProds.begin() + i);
						decProds.insert(decProds.end(), secDecProds.begin(), secDecProds.end());
					}
				}
			}
			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}
	}

	return vec;
}

std::vector<Molecola> decomposeLinEthRO(Molecola mol)
{
	std::vector<Molecola> vec;
	int dist = mol.dist(mol.posCrad(), mol.trova(8));
	if (dist == 1)
	{
		int found1[SIZEMAX + 1];
		int found2[SIZEMAX + 1];
		mol.scorri(mol.trova(8), 1, found1);
		mol.scorri(mol.posCrad(), 2, found2);
		int pos = 0;
		for (int i = 0; i < SIZEMAX + 1; i++)
			if (found1[i] == 1 && found2[i] == 1)
				pos = i;
		Molecola m1, m2;
		mol.spezza(pos, mol.trova(8), &m1, &m2);
		std::vector<Molecola> prods = decomposeRO(m1);
		vec.push_back(m2);
		for (int i = 0; i < prods.size(); i++)
			vec.push_back(prods[i]);
		//for (int i = 0; i < vec.size(); i++)
		//	std::cout << vec[i] << std::endl;
	}
	else if (dist == 2)
	{
		int found1[SIZEMAX + 1];
		int found2[SIZEMAX + 1];
		mol.scorri(mol.trova(8), 1, found1);
		mol.scorri(mol.posCrad(), 1, found2);
		int pos = 0;
		for (int i = 0; i < SIZEMAX + 1; i++)
			if (found1[i] == 1 && found2[i] == 1)
				pos = i;
		Molecola m1, m2;
		mol.spezza(pos, mol.trova(8), &m1, &m2);
		m2.scorri(m2.posCrad(), 1, found1);
		for (int i = 0; i < SIZEMAX + 1; i++)
			if (found1[i] == 1)
				pos = i;
		m2.removeAtom(m2.posCrad());
		m2.C_to_COrad(pos);
		Molecola m3, m4;
		if(m2.dist(m2.trova(9), m2.trova(7)) == 1)
			m2.spezza(m2.trova(9), m2.trova(7), &m3, &m4);
		if (m2.dist(m2.trova(9), m2.trova(7)) == 2)
		{
			int found3[SIZEMAX + 1];
			m2.scorri(m2.trova(7), 1, found3);
			int pos1 = 0;
			for (int i = 0; i < SIZEMAX + 1; i++)
				if (found3[i] == 1 && m2.dist(m2.trova(9), i) == 1)
					pos1 = i;
			m2.spezza(m2.trova(9), pos1, &m3, &m4);
		}
		//std::cout << mol << std::endl << m1 << std::endl << m3 << std::endl << m4 << std::endl;
		vec.push_back(m1);
		vec.push_back(m3);
		vec.push_back(m4);
	}

	return vec;
}

std::vector<Molecola> decomposeCEthR(Molecola m2)
{
	Molecola m3, m4;
	std::vector<Molecola> vec;
	// decompose the radical cyclic ether
	if (m2.posCrad() == m2.posEthero()[0] || m2.posCrad() == m2.posEthero()[1]) // if the radical is on one of the carbons attached to the cyclic O	
	{
		int pos1 = 0;	// oxygen bond that doubles
		int pos2 = 0;	// oxygen bond that breaks
		if (m2.posCrad() == m2.posEthero()[0])
		{
			pos1 = m2.posEthero()[0];
			pos2 = m2.posEthero()[1];
		}
		else
		{
			pos1 = m2.posEthero()[1];
			pos2 = m2.posEthero()[0];
		}

		//TEST

		m2.removeCycEther(1);
		m2.Crad_to_C(pos1);
		m2.addcheto(pos1);
		m2.C_to_Crad(pos2);
		std::vector<Molecola> prods = decomposeRO(m2);
		for (int k = 0; k < prods.size(); k++)
			vec.push_back(prods[k]);
	}
	//if (m2.posCrad() == m2.posEthero()[0] || m2.posCrad() == m2.posEthero()[0])
	else
	{
		int found[SIZEMAX + 1];
		int pos1 = 0;		// position of the radical
		int pos2 = 0;		// position of the oxygen bond that breaks
		int pos3 = 0;		// position of the oxygen bond that becomes double bond

		m2.scorri(m2.posCrad(), 1, found);
		for (int k = 0; k < SIZEMAX + 1; k++)
		{
			if (found[k] == 1 && k == m2.posEthero()[0])
			{
				pos1 = m2.posCrad();
				pos2 = m2.posEthero()[0];
				pos3 = m2.posEthero()[1];
			}
			if (found[k] == 1 && k == m2.posEthero()[1])
			{
				pos1 = m2.posCrad();
				pos2 = m2.posEthero()[1];
				pos3 = m2.posEthero()[0];
			}
		}
		if (pos1 != 0 && pos2 != 0 && pos3 != 0)	// radical is in alpha with respect of one carbon bonded to the oxygen
		{
			m2.Crad_to_C(pos1);
			m2.removeCycEther(1);
			m2.addole(pos1, pos2);
			m2.addcheto(pos3);
			//Molecola m3, m4;
			m2.scorri(pos3, 1, found);
			int pos4 = 0;
			for (int k = 0; k < SIZEMAX + 1; k++)
			{
				if (found[k] == 1)
				{
					if (pos4 == 0)
						pos4 = k;
					else if (m2.dist(k, pos2) < m2.dist(pos4, pos2))	// select the carbon nearest to the double bond
						pos4 = k;
				}
			}
			m2.spezza(pos3, pos4, &m3, &m4);
			int posRadm3 = m3.posCrad();
			m3.Crad_to_C(posRadm3);
			m3.addcheto(posRadm3);
			vec.push_back(m3);
			vec.push_back(m4);
		}
		else    // radical is in beta with respect of one of the carbons attached to the oxygen
		{
			m2.scorri(m2.posCrad(), 2, found);
			for (int k = 0; k < SIZEMAX + 1; k++)
			{
				if (found[k] == 1 && k == m2.posEthero()[0])
					pos1 = m2.posEthero()[0];
				if (found[k] == 1 && k == m2.posEthero()[1])
					pos1 = m2.posEthero()[1];
			}
			if (pos1 != 0)
			{
				m2.scorri(m2.posCrad(), 1, found);
				int found2[SIZEMAX + 1];
				m2.scorri(pos1, 1, found2);
				for (int k = 0; k < SIZEMAX + 1; k++)
					if (found[k] == 1 && found2[k] == 1)
						pos2 = k;
				//Molecola m3, m4;
				m2.spezza(pos1, pos2, &m3, &m4);
				//std::cout << m3 << std::endl << m4 << std::endl << std::endl;
				if (m4.size() == 0) // m3 is a linear ether 
				{
					Molecola m5, m6;
					m3.scorri(m3.posCrad(), 2, found);
					m3.scorri(m3.trova(8), 1, found2);
					for (int k = 0; k < SIZEMAX + 1; k++)
						if (found[k] == 1 && found2[k] == 1)
							pos1 = k;
					m3.spezza(pos1, m3.trova(8), &m5, &m6);
					//std::cout << m5 << std::endl << m6 << std::endl;
					vec.push_back(m5);
					vec.push_back(m6);
				}
				else
				{
					std::vector<Molecola> prod = decomposeCEthR(m3);
					vec.push_back(m4);
					vec.push_back(prod[0]);
					if (prod.size() == 2)
						vec.push_back(prod[1]);
				}
			}
		}
	}
	return vec;
}


//==============================================================
//------------------------MAIN----------------------------------
//==============================================================

int main(int argc, char* argv[])
{
	//DEBUG
	//{
	//	char nome1[100] = "molecules/test1.hyd";
	//	Molecola mol1(nome1);
	//
	//	char nome2[100] = "molecules/test2.hyd";
	//	Molecola mol2(nome2);
	//
	//	std::cout << mol1 << std::endl;
	//	std::cout << mol2 << std::endl;
	//	if (mol1 == mol2)
	//		std::cout << "they are equal" << std::endl;
	//	else
	//		std::cout << "They are different" << std::endl;
	//
	//	return 0;
	//} 


	//DEBUG
	// setup the console window
	std::string moleculesSubDir = "C7-C16";
	std::vector<std::string> molecsPaths = {
											 "2_2-m-c4",
											 "2_3-m-c4",
											 "2-m-c4",
											 "2-m-c5",
											 "3-m-c5",
											 "neoC5",
											 "nc5",
											 "nc6",
											 "nc7",
											 //"2-m-c7",
											 //"nc8",
											 //"nc9",
											 //"nc10",
											 //"nc11",
											 //"nc12",
											 //"nc13",
											 //"nc14",
											 //"nc15",
											 //"nc16",
											 //"2_2_3-m-c4",
											 //"2_2_4_4-m-c5",
											 //"2_2_4_6_6-m-C7",
											 //"2_2_4-m-c5",
											 //"2_2-m-c5",
											 //"2_3-m-c5",
											 //"2_4-m-c5",
											 //"2_4-m-c7",
											 //"2_5-m-c6",
											 //"2_6-m-c7",
											 //"2_7-m-c8",
											 //"2_9-m-C10",
											 //"2-m-c6",
											 //"2-m-c8",
											 //"2-m-c9",
											 //"2-m-c10",
											 //"2-m-C11",
											 //"3_3-m-c5",
											 //"3-m-c6",
											 //"3-m-c7",
											 //"isocetane"
	};


	HWND console = GetConsoleWindow();
	RECT rec;
	GetWindowRect(console, &rec);
	MoveWindow(console, rec.left, rec.top, 1500, 800, TRUE);

	// check the presence of required resources
	//if (!thereIsDirectory("resources"))
	//{
	//	std::cerr << "ERROR: resources folder is missing!" << std::endl;
	//	return -1;
	//}
	//if (!thereIsDirectory("resources"))
	//{
	//	std::cerr << "ERROR: resources folder is missing!" << std::endl;
	//	return -1;
	//}


	for (int qwer = 0; qwer < molecsPaths.size(); qwer++)
	{

		// generate all the required directories
		_mkdir("OUTPUT");
		_mkdir("resources/temp");

		int tipoVbr = 1;        // type of equation for computing the branching speed
		int errflg = 0;

		//---------------------initialize with zeros the global variables tipor and tipoqooh
		{
			int a, b;
			for (a = 0; a <= SIZEMAX; a++)
			{
				for (b = 0; b <= SIZEMAX; b++) tipoqooh[a][b] = 0;
				tipor[a] = 0;
			};
		};

		//---------------------opens the output file for printing the reactions
		reaout.open("OUTPUT/reactions.out");
		if (!reaout)
		{
			cerr << "Error in opening the file: reactions.out" << endl << endl;
			exit(1);
		};


		////////////////////////////////////////////////////////////////
		//  kinetic data input from file
		////////////////////////////////////////////////////////////////
		//Kinox k((char*)"kinox.kin");
		std::string kinoxInput = "resources/rateRules/kinetic_constants.csv";
		Kinox k(kinoxInput);
		std::cout << "Kinetic constants read form " << kinoxInput << ": \n\n" << k << endl;
		////////////////////////////////////////////////////////////////
		//  molecular data input from file
		////////////////////////////////////////////////////////////////


		//strcpy(nome, "molecules/nc7.hyd");
		//strcpy(nome, "molecules/2_2_3-m-c4.hyd");

		/*
		bool rightAnswer = false;
		bool useBatch = false;
		while (rightAnswer == false)
		{
			std::cout << "Batch or CSTR simulation? (1-batch, 2-CSTR):" << std::endl;
			int typedNumber;
			std::cin >> typedNumber;
			if (typedNumber == 1)
			{
				useBatch = true;
				rightAnswer = true;
			}
			else if (typedNumber == 2)
			{
				useBatch = false;
				rightAnswer = true;
			}
			else
				std::cout << typedNumber << " is not a valid input!" << std::endl;
		}

		char nome[100];
		std::cout << "Type name of molecule:" << std::endl;
		std::string molNameInp;
		std::cin >> molNameInp;
		//molNameInp = "nc6";
		*/

		bool useBatch = true;
		std::string molNameInp = molecsPaths[qwer];
		char nome[100];
		std::string molDirectory = "molecules/";
		if (moleculesSubDir != "")
		{
			molDirectory.append(moleculesSubDir);
			molDirectory.append("/");
		}
		molDirectory.append(molNameInp);
		molDirectory.append(".hyd");
		strcpy(nome, molDirectory.c_str());
		Molecola HC(nome);
		HC.output = MATRICE;
		std::cout << "Molecule read from " << nome << ": \n\n" << HC << endl;
		HC.output = FORMULA;
		std::cout << HC << endl << endl;

		k.sizeReactions = HC.size() * 9;


		chemout.setFuel(&HC);

		reaout << "Primary oxidation reactions for: " << nome << endl;   // ITA: Reazioni primarie per ossidazione di
		reaout << "-------------------------------------------------------" << endl;
		reaout << HC << endl;
		reaout << "-------------------------------------------------------" << endl;
		reaout << endl;

		{
			std::stringstream ssBuff;
			ssBuff << "!" << "Primary oxidation reactions for: " << nome << std::endl;
			ssBuff << "!" << "-------------------------------------------------------" << std::endl;
			ssBuff << "!" << HC << std::endl;
			ssBuff << "!" << "-------------------------------------------------------" << std::endl;
			ssBuff << std::endl;
			chemout.writeHeadingDetailed(&ssBuff);
		}
		{
			std::stringstream ssBuff;
			ssBuff << "!" << "Lumped oxidation reactions for: " << nome << std::endl;
			ssBuff << "!" << "-------------------------------------------------------" << std::endl;
			ssBuff << "!" << HC << std::endl;
			ssBuff << "!" << "-------------------------------------------------------" << std::endl;
			ssBuff << std::endl;
			chemout.writeHeadingLumped(&ssBuff);
		}


		////////////////////////////////////////////////////////////////
		Molecola* r;
		Molecola* roo;
		Molecola* qooh;
		Molecola* ooqooh;
		r = new Molecola[HC.size() + 1]; if (!r) exit(200);
		roo = new Molecola[HC.size() + 1]; if (!roo) exit(201);
		qooh = new Molecola[HC.size() * 6 + 1]; if (!qooh) exit(202);
		ooqooh = new Molecola[HC.size() * 6 + 1]; if (!ooqooh) exit(203);
		int numr = 0, numqooh = 0;

		// compute reactor inputs
		double tau_CSTR = 2.0;
		double molesFuelIN_CSTR = 1.0;
		double equivalenceRatio_CSTR = 1.0;
		double molesO2IN_CSTR = equivalenceRatio_CSTR * molesFuelIN_CSTR * (double(HC.numberOfC()) + double(HC.numberOfH()) / 4.0 - double(HC.numberOfO()) / 2.0);
		double molesN2IN_CSTR = molesO2IN_CSTR / 0.21 * 0.79;
		double totMolesIN_CSTR = molesFuelIN_CSTR + molesN2IN_CSTR + molesO2IN_CSTR;
		double xHC_IN = molesFuelIN_CSTR / totMolesIN_CSTR;

		std::cout << "Moles fuel : " << molesFuelIN_CSTR << std::endl;
		std::cout << "Moles O2   : " << molesO2IN_CSTR << std::endl;
		std::cout << "Moles N2   : " << molesN2IN_CSTR << std::endl;

		// TEST
		//double tau = 0.002;
		//molesFuelIN = 0.01;
		//molesO2IN = 0.08;
		//molesN2IN = 0.91;
		//totMolesIN = molesFuelIN + molesN2IN + molesO2IN;
		//xHC_IN = molesFuelIN / totMolesIN;
		// TEST


		decomponi(HC, k, r, roo, qooh, ooqooh, &numr, &numqooh);

		int numtot = numr * 2 + numqooh * 2;

		////////////////////////////////////////////////////////////////
		print1(numr, numqooh, r, roo, qooh, ooqooh);
		////////////////////////////////////////////////////////////////

		// generate special molecules
		Molecola N2(1);
		Molecola O2(2);
		Molecola O(3);
		Molecola OH(4);
		Molecola HO2(5);
		Molecola H2O(6);
		Molecola H2O2(7);
		Molecola H(8);
		Molecola H2(9);
		Molecola CO(10);
		Molecola CH3;
		CH3.makeCH3();
		Molecola CH4;
		CH4.makeCH4();
		Molecola C2H5;
		C2H5.makeC2H5();
		Molecola C2H6;
		C2H6.makeC2H6();
		Molecola CH3OO;
		CH3OO.makeCH3OO();
		Molecola CH3OOH;
		CH3OOH.makeCH3OOH();


		// list of reactions
		std::vector<Reaction>reactionsList;





		//////////////Filling the non linear systems////////////
		//
		//                          r'
		//                          |
		//          ole + ho2  <--  r  -->  r' + ole
		//                          |
		//                          |
		//              rooh  <--  roo
		//                          |
		//                          |  ,--> beta-dec
		//              rooh  <--  qooh  --> etero + oh
		//                          |  `--> ole + ho2
		//                          |
		//                        ooqooh
		//                          |
		//                          V
		//                      kehy + oh
		//
		////////////////////////////////////////////////////////////////

		std::cout << std::endl << "*****************************************************" << std::endl;
		std::cout << "Filling reaction list:" << std::endl;

		//-------------- initiation reactions (C-C breakage) ---------------
		reaout << "Initiation reactions:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Initiation reactions:");
		{
			Molecola* firstRad;			// array of molecules where the first product is stored
			Molecola* secondRad;		// array of molecules where the second product is stored
			Molecola temp1, temp2;		// molecules where to temporary save the products
			firstRad = new Molecola[HC.size() + 1]; if (!firstRad) exit(200);
			secondRad = new Molecola[HC.size() + 1]; if (!secondRad) exit(201);

			// search for the bonds that can be broken
			std::vector < std::vector<int> > atomsInBond;
			atomsInBond = HC.listOfCCBonds();
			//for (int i = 0; i < atomsInBond.size(); i++)
			//	std::cout << atomsInBond[i][0] << "   " << atomsInBond[i][1] << std::endl;


			int numberOfReactions = 0;
			for (int i = 0; i < atomsInBond.size(); i++)
			{
				bool alreadyPresent = false;			// tells if the an equivalent reaction is already present
				int isPossible = HC.spezza(atomsInBond[i][0], atomsInBond[i][1], &temp1, &temp2);
				if (isPossible == 0) break;

				for (int j = 1; j <= HC.size(); j++)
				{
					if ((temp1 == firstRad[j] && temp2 == secondRad[j]) || (temp2 == firstRad[j] && temp1 == secondRad[j])) // check if an equivalent reaction is already saved
					{
						firstRad[j].isomeri++;				// the number of equivalent reaction is stored in the firsRad.isomeri just for convinience
						alreadyPresent = true;
						break;
					}
					int asd = 1;
				}
				if (!alreadyPresent)
				{
					numberOfReactions++;
					firstRad[numberOfReactions] = temp1;
					secondRad[numberOfReactions] = temp2;
				}
			}

			Molecola* allRad;			// list of all the radical generated by the initiation reactions
			int maxSizeAllRad = HC.size() * 3;
			allRad = new Molecola[maxSizeAllRad]; if (!allRad) exit(200);
			int elementsInAllRad = 0;

			for (int i = 1; i <= numberOfReactions; i++)
			{
				Radicale firstRadicalType = firstRad[i].tipoR(firstRad[i].trova(2));
				Radicale secondRadicalType = secondRad[i].tipoR(secondRad[i].trova(2));

				int positionRad1 = Add(firstRad[i], allRad, &elementsInAllRad, maxSizeAllRad);			// check if firstRad is already present in allRad, if not add it
				int positionRad2 = Add(secondRad[i], allRad, &elementsInAllRad, maxSizeAllRad);			// check if secondRad is already present in allRad, if not add it

				//for (int j = 0; j < Temps.size(); j++)
				//{
				//	double velocity = k.v_initiation(firstRadicalType, secondRadicalType, firstRad[i].isomeri, Temps[j]);
				//	allRad[positionRad1].velo_v[j] += velocity;	// update the formation velocity of radical 1
				//	allRad[positionRad2].velo_v[j] += velocity;	// update the formation velocity of radical 2
				//}
				std::string reacomm = k.v_initiation(firstRadicalType, secondRadicalType, firstRad[i].isomeri);
				addReactionToList(&reactionsList, HC, std::vector<Molecola>{ firstRad[i], secondRad[i] }, new double[3] { k.A, k.n, k.E }, "Initiation reactions:", reacomm);
				k.wrirea(reaout, 2, k.A, k.n, k.E, HC, firstRad[i], secondRad[i]);
				//chemout.wrireaDetailed(2, k.A, k.n, k.E, HC, firstRad[i], secondRad[i]);

			}

			std::cout << "Initiation reactions added." << std::endl;
		}


		//------------ H abstractions------------------------
		reaout << "H abstraction: " << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("H abstraction:");

		//std::vector<Reaction>* reacVec, Molecola HabsMol, Molecola HabsMolProd
		h_abstraction(k, numr, HC, r, o2, &reactionsList, O2, HO2);
		h_abstraction(k, numr, HC, r, oh, &reactionsList, OH, H2O);
		h_abstraction(k, numr, HC, r, h, &reactionsList, H, H2);
		h_abstraction(k, numr, HC, r, o, &reactionsList, O, OH);
		h_abstraction(k, numr, HC, r, ho2, &reactionsList, HO2, H2O2);
		h_abstraction(k, numr, HC, r, ch3, &reactionsList, CH3, CH4);
		h_abstraction(k, numr, HC, r, c2h5, &reactionsList, C2H5, C2H6);
		h_abstraction(k, numr, HC, r, ch3oo, &reactionsList, CH3OO, CH3OOH);
		lumpout << std::endl;

		std::cout << "H-abstraction reactions added." << std::endl;

		//------------O2 addition and inverse reaction------------------------- ITA: Addizione di O2 e inversa
		reaout << "O2 addition:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("O2 addition:");
		for (int i = 1; i <= numr; i++)
		{
			Radicale radicalType = r[i].tipoR(r[i].trova(2));
			std::string reacomm = k.v_o2_add_r(radicalType);
			k.wrirea(reaout, 3, k.A, k.n, k.E, r[i], roo[i]);
			addReactionToList(&reactionsList, std::vector<Molecola>{ r[i], O2 }, roo[i], new double[3] { k.A, k.n, k.E }, "O2 addition to R:", reacomm);
			//chemout.wrireaDetailed( 3, k.A, k.n, k.E, r[i], roo[i]);
		}
		for (int i = 1; i <= numr; i++)
		{
			Radicale radicalType = r[i].tipoR(r[i].trova(2));
			std::string reacomm = k.v_o2_rem_roo(radicalType, r[i].numberOfC());
			addReactionToList(&reactionsList, roo[i], std::vector<Molecola>{ r[i], O2 }, new double[3] { k.A, k.n, k.E }, "O2 elimination from ROO:", reacomm);
			k.wrirea(reaout, 4, k.A, k.n, k.E, roo[i], r[i]);
			//chemout.wrireaDetailed(4, k.A, k.n, k.E, roo[i], r[i]);
		}
		for (int i = 1; i <= numqooh; i++)
		{
			Radicale radicalType = qooh[i].tipoR(qooh[i].trova(2));
			std::string reacomm = k.v_o2_add_qooh(radicalType);
			addReactionToList(&reactionsList, std::vector<Molecola>{ qooh[i], O2 }, ooqooh[i], new double[3] { k.A, k.n, k.E }, "O2 addition to QOOH:", reacomm);
			k.wrirea(reaout, 3, k.A, k.n, k.E, qooh[i], ooqooh[i]);
			//chemout.wrireaDetailed(3, k.A, k.n, k.E, qooh[i], ooqooh[i]);
		}
		for (int i = 1; i <= numqooh; i++)
		{
			Radicale radicalType = qooh[i].tipoR(qooh[i].trova(2));
			std::string reacomm = k.v_o2_rem_ooqooh(radicalType, qooh[i].numberOfC());
			addReactionToList(&reactionsList, ooqooh[i], std::vector<Molecola>{ qooh[i], O2 }, new double[3] { k.A, k.n, k.E }, "O2 elimination from OOQOOH:", reacomm);
			k.wrirea(reaout, 4, k.A, k.n, k.E, ooqooh[i], qooh[i]);
			//chemout.wrireaDetailed( 4, k.A, k.n, k.E, ooqooh[i], qooh[i]);
		}
		reaout << endl;

		std::cout << "O2 addition reactions added." << std::endl;

		//------------Isomerization r -> r'---------------------------

		reaout << "Isomerization R:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Isomerization R:");
		for (int i = 1; i <= numr; i++)   // i scorre i radicali
		{
			int j;
			int pos_rad = pos_r(i);
			int trovati[SIZEMAX + 1];
			//-----------------------------isomerization involving 5 atom ring
			r[i].scorri(pos_rad, 3, trovati);		// flag in the vector "trovati" all the carbons that are at 3 atom distance from the radical
			for (j = 1; j <= SIZEMAX; j++)
				if (trovati[j] == 1 && r[i].numAbstractableH(j) != 0)
				{
					std::string reacomm = k.v_isomerization_r(r[i].tipoR(pos_rad), r[i].tipoH(j), a5, r[i].numAbstractableH(j));
					addReactionToList(&reactionsList, r[i], r[tipor[j]], new double[3] { k.A, k.n, k.E }, "Isomerization R:", reacomm);
					k.wrirea(reaout, 1, k.A, k.n, k.E, r[i], r[tipor[j]]);
					//chemout.wrireaDetailed(1, k.A, k.n, k.E, r[i], r[tipor[j]]);
				};
			//-----------------------------isomerization involving 6 atom ring
			r[i].scorri(pos_rad, 4, trovati);
			for (j = 1; j <= SIZEMAX; j++)
				if (trovati[j] == 1 && r[i].numAbstractableH(j) != 0)
				{
					std::string reacomm = k.v_isomerization_r(r[i].tipoR(pos_rad), r[i].tipoH(j), a6, r[i].numAbstractableH(j));
					addReactionToList(&reactionsList, r[i], r[tipor[j]], new double[3] { k.A, k.n, k.E }, "Isomerization R:", reacomm);
					k.wrirea(reaout, 1, k.A, k.n, k.E, r[i], r[tipor[j]]);
					//chemout.wrireaDetailed(1, k.A, k.n, k.E, r[i], r[tipor[j]]);
				};
			//-----------------------------isomerization involving 7 atom ring
			r[i].scorri(pos_rad, 5, trovati);
			for (j = 1; j <= SIZEMAX; j++)
				if (trovati[j] == 1 && r[i].numAbstractableH(j) != 0)
				{
					std::string reacomm = k.v_isomerization_r(r[i].tipoR(pos_rad), r[i].tipoH(j), a7, r[i].numAbstractableH(j));
					addReactionToList(&reactionsList, r[i], r[tipor[j]], new double[3] { k.A, k.n, k.E }, "Isomerization R:", reacomm);
					k.wrirea(reaout, 1, k.A, k.n, k.E, r[i], r[tipor[j]]);
					//chemout.wrireaDetailed(1, k.A, k.n, k.E, r[i], r[tipor[j]]);
				};
		};
		reaout << endl;
		std::cout << "R isomerization added." << std::endl;


		//------------Isomerizzazione roo -> qooh-----------------------

		reaout << "Isomerization ROO:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Isomerization ROO:");
		for (int i = 1; i <= numr; i++)
		{
			int j;
			int pos_o2 = roo[i].trova(3);
			int trovati[SIZEMAX + 1];

			for (int dist = 1; dist <= 4; dist++)
			{
				Anello ringSize;
				switch (dist)
				{
				case 1:
					ringSize = a5;
					break;
				case 2:
					ringSize = a6;
					break;
				case 3:
					ringSize = a7;
					break;
				case 4:
					ringSize = a8;
					break;
				default:
					break;
				}
				roo[i].scorri(pos_o2, dist, trovati);
				for (j = 1; j <= SIZEMAX; j++)          // j iterates trough the found H
					if (trovati[j] == 1 && roo[i].numAbstractableH(j) != 0)
					{
						std::string reacomm = k.v_isom_roo(roo[i].tipoROO(roo[i].trova(3)), roo[i].tipoH(j), ringSize, roo[i].numAbstractableH(j));
						int pos_ooh = pos_r(i);
						int pos_r = j;
						addReactionToList(&reactionsList, roo[i], qooh[tipoqooh[pos_ooh][pos_r]], new double[3] { k.A, k.n, k.E }, "Isomerization ROO:", reacomm);
						k.wrirea(reaout, 1, k.A, k.n, k.E, roo[i], qooh[tipoqooh[pos_ooh][pos_r]]);
						//chemout.wrireaDetailed(1, k.A, k.n, k.E, roo[i], qooh[tipoqooh[pos_ooh][pos_r]]);
					};
			}
		};  // end for i
		reaout << endl;
		std::cout << "ROO isomerization added." << std::endl;

		//------------Isomerization OOQOOH -> P(OOH)2-----------------------

		reaout << "Isomerization OOQOOH:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Isomerization OOQOOH:");

		std::vector<Molecola> POOH2_isom_reac;		// vectors for printing the inverse reaction later
		std::vector<Molecola> POOH2_isom_prod;
		std::vector<int> POOH2_isom_OOH_pos;
		std::vector<int> POOH2_isom_rad_pos;
		std::vector<Anello> POOH2_isom_ring;

		std::vector<Molecola> pooh2;

		std::vector<Reaction> reactionCount;		// vector in which teh reactions are stored in order to check if there are duplicates
		// in the case a reaction is duplicate it is not saved in the vectors of the reverse reaction
		// because otherwise in the following reactions two equal reactions would be formed.

		// debug
		//for (int i = 1; i <= numqooh; i++)
		//{
		//	std::cout << qooh[i] << "      " << ooqooh[i] << std::endl;
		//}
		//std::cout << std::endl << std::endl;

		for (int i = 1; i <= numqooh; i++)
		{
			int j;
			int pos_o2 = ooqooh[i].trova(3);
			int trovati[SIZEMAX + 1];
			//std::cout << ooqooh[i] << std::endl;

			for (int distance = 1; distance < 5; distance++)
			{
				//std::cout << "  Distance = " << distance << std::endl;
				ooqooh[i].scorri(pos_o2, distance, trovati);
				Anello ring = a5;
				switch (distance)
				{
				case 1:
					ring = a5;
					break;
				case 2:
					ring = a6;
					break;
				case 3:
					ring = a7;
					break;
				case 4:
					ring = a8;
					break;
				default:
					break;
				}
				for (j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H

					if (trovati[j] == 1 && ooqooh[i].numAbstractableH(j) != 0 && ooqooh[i].tipo(j) != 4)
					{
						std::string reacomm = k.v_isom_ooqooh(ooqooh[i].tipoROO(ooqooh[i].trova(3)), ooqooh[i].tipoH(j), ring, ooqooh[i].numAbstractableH(j));
						Molecola mol = ooqooh[i];
						mol.OO_to_OOH(ooqooh[i].trova(3));
						mol.C_to_Crad(j);
						Add(mol, &pooh2);
						//std::cout << "          " << mol << std::endl;
						addReactionToList(&reactionsList, ooqooh[i], mol, new double[3] { k.A, k.n, k.E }, "Isomerization OOQOOH:", reacomm);
						k.wrirea(reaout, 1, k.A, k.n, k.E, ooqooh[i], mol);
						int size = reactionCount.size() - 1;
						if (size != Add(reactionsList[reactionsList.size() - 1], &reactionCount))	// if the new reaction was already present do not add it to the inverse reaction
						{
							POOH2_isom_reac.push_back(mol);
							POOH2_isom_prod.push_back(ooqooh[i]);
							POOH2_isom_OOH_pos.push_back(ooqooh[i].trova(3));
							POOH2_isom_rad_pos.push_back(j);
							POOH2_isom_ring.push_back(ring);

						}
					};
			}
		};  // end for i
		reaout << endl;
		std::cout << "OOQOOH isomerization added." << std::endl;

		//------------Isomerization  P(OOH)2  ->  OOQOOH -----------------------

		reaout << "Isomerization P(OOH)2:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Isomerization P(OOH)2:");
		for (int i = 0; i < POOH2_isom_reac.size(); i++)
		{
			std::string reacomm = k.v_isom_pooh2(POOH2_isom_reac[i].tipoROOH(POOH2_isom_OOH_pos[i]), POOH2_isom_reac[i].tipoR(POOH2_isom_rad_pos[i]), POOH2_isom_ring[i]);
			addReactionToList(&reactionsList, POOH2_isom_reac[i], POOH2_isom_prod[i], new double[3] { k.A, k.n, k.E }, "Isomerization P(OOH)2:", reacomm);
			k.wrirea(reaout, 1, k.A, k.n, k.E, POOH2_isom_reac[i], POOH2_isom_prod[i]);
		}
		reaout << endl;
		std::cout << "P(OOH)2 isomerization added." << std::endl;

		//------------Olefine parenti da   roo -> ole + ho2 -----------
		reaout << "Olefins parent from ROO:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Olefins parent from ROO:");

		std::vector<Molecola> olefins;			// vector to store all the produced olefins (the one of the same size of the fuel)

		for (int i = 1; i <= numr; i++)
		{
			int j;
			int pos_o2 = roo[i].trova(3);
			int trovati[SIZEMAX + 1];

			roo[i].scorri(pos_o2, 1, trovati);	// find the carbon at distance 1 from the carbon with the OO
			for (j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
				if (trovati[j] == 1 && roo[i].numAbstractableH(j) != 0)
				{
					if (roo[i].tipoROO(roo[i].trova(3)) == Rp && roo[i].tipoC(j) == Cp);			// CHECK if it is right to neglect those reactions
					Molecola ole = HC;
					ole.addole(roo[i].trova(3), j);
					std::string reacomm = k.v_ole_par_roo(roo[i].tipoROO(roo[i].trova(3)), roo[i].tipoC(j), roo[i].numAbstractableH(j));
					int pos_ooh = pos_r(i);
					int pos_r = j;
					Add(ole, &olefins);		// add the olefin to the vector
					addReactionToList(&reactionsList, roo[i], std::vector<Molecola>{ ole, HO2 }, new double[3] { k.A, k.n, k.E }, "Olefins parent from ROO:", reacomm);
					k.wrirea(reaout, 6, k.A, k.n, k.E, roo[i], ole);
					//chemout.wrireaDetailed(1, k.A, k.n, k.E, roo[i], qooh[tipoqooh[pos_ooh][pos_r]]);
				};

		};  // end for i
		reaout << endl;
		std::cout << "Olefins parent from ROO added." << std::endl;

		//------------Isomerizzazione qooh -> roo-----------------------
		reaout << "Isomerization QOOH:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Isomerization QOOH:");
		for (int i = 1; i <= numqooh; i++)
		{
			Pos p = pos_qooh(i);
			Anello a;
			if (qooh[i].dist(p.rad, p.ooh) == 1) a = a5;
			else if (qooh[i].dist(p.rad, p.ooh) == 2) a = a6;
			else if (qooh[i].dist(p.rad, p.ooh) == 3) a = a7;
			else continue;
			std::string reacomm = k.v_isom_qooh(qooh[i].tipoROOH(p.ooh), qooh[i].tipoR(qooh[i].trova(2)), a, qooh[i].numberOfC());
			addReactionToList(&reactionsList, qooh[i], roo[tipor[p.ooh]], new double[3] { k.A, k.n, k.E }, "Isomerization QOOH : ", reacomm);
			k.wrirea(reaout, 1, k.A, k.n, k.E, qooh[i], roo[tipor[p.ooh]]);
			//chemout.wrireaDetailed(1, k.A, k.n, k.E, qooh[i], roo[tipor[p.ooh]]);
		};
		reaout << endl;
		std::cout << "QOOH isomerization added." << std::endl;


		int max1 = HC.size() * 2 - 1;
		int max2 = HC.size() - 1 + 1;
		int max3 = numqooh + 1;
		int max4 = numqooh + 1;
		int max5 = HC.size() - 1 + 1;
		int max60 = numqooh + 1; //rco
		int max61 = numqooh + 1;    //ald
		int max7a = HC.size() * 2;              // non sicuro
		int max7b = HC.size() * 4;              // non sicuro

		std::string reactionName;

		Molecola* rbeta;
		rbeta = new Molecola[max1];    if (!rbeta)  err((char*)"No Mem !!");
		Molecola* olebetar;
		olebetar = new Molecola[max1]; if (!olebetar)  err((char*)"No Mem !!");
		Molecola* allProducts;
		allProducts = new Molecola[2 * max1]; if (!olebetar)  err((char*)"No Mem !!");

		Molecola* olepar_r;
		olepar_r = new Molecola[max2]; if (!olepar_r)  err((char*)"No Mem !!");
		double ho2_r = 0;

		Molecola* olepar_qooh;
		olepar_qooh = new Molecola[max5]; if (!olepar_qooh)  err((char*)"No Mem !!");
		double ho2_qooh = 0;

		Molecola* etero;
		etero = new Molecola[max4];        if (!etero)  err((char*)"No Mem !!");
		double oh_etero = 0;

		Molecola* olebetaqooh;
		olebetaqooh = new Molecola[max3];  if (!olebetaqooh)  err((char*)"No Mem !!");
		Molecola* aldbetaqooh;
		aldbetaqooh = new Molecola[max3];  if (!aldbetaqooh)  err((char*)"No Mem !!");
		double oh_betaqooh = 0;

		Molecola* rco_branch;
		rco_branch = new Molecola[max60];   if (!rco_branch)  err((char*)"No Mem !!");
		Molecola* ald_branch;
		ald_branch = new Molecola[max61];   if (!ald_branch)  err((char*)"No Mem !!");
		double oh_branch = 0;

		Molecola* co_omo1;
		co_omo1 = new Molecola[max7a];
		double ho2_omo1 = 0;

		Molecola* co_omo2;
		co_omo2 = new Molecola[max7b];
		double oh_omo2 = 0;

		//------------Beta decomp. radicali  r -> r' + ole -------------

		int elementsInAllProducts = 0;

		reaout << "Radicals beta decomposition" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Radicals beta decomposition:");
		int num10 = 0; int num11 = 0;
		for (int i = 1; i <= numr; i++)
		{
			int pos_rad = pos_r(i);
			int posalfa[SIZEMAX + 1];
			r[i].scorri(pos_rad, 1, posalfa);
			for (int j = 1; j <= SIZEMAX; j++) if (posalfa[j] == 1) // j scorre le pos alfa
			{
				int posbeta[SIZEMAX + 1];
				r[i].scorri(j, 1, posbeta);  // per ogni alfa trova le beta
				for (int z = 1; z <= SIZEMAX; z++) if (posbeta[z] == 1 && z != pos_rad)
				{
					Molecola m1, m2;
					int isPossible = r[i].spezza(j, z, &m1, &m2);
					if (isPossible == 0) continue;  // MODIFICATION

					int a = Add(m1, rbeta, &num10, max1);
					int b = Add(m2, olebetar, &num11, max1);
					int c = Add(m1, allProducts, &elementsInAllProducts, 2 * max1);
					int d = Add(m2, allProducts, &elementsInAllProducts, 2 * max1);
					std::string reacomm = k.v_beta_dec_r(r[i].tipoR(pos_rad), m2.tipoR(m2.trova(2)));
					addReactionToList(&reactionsList, r[i], std::vector<Molecola>{ rbeta[a], olebetar[b] }, new double[3] { k.A, k.n, k.E }, "Radicals beta decomposition", reacomm);
					k.wrirea(reaout, 2, k.A, k.n, k.E, r[i], rbeta[a], olebetar[b]);
					//chemout.wrireaDetailed(2, k.A, k.n, k.E, r[i], rbeta[a], olebetar[b]);
				};		// end pos beta
			};			// end pos alfa
		};
		reaout << endl;
		std::cout << "R beta decomposition added." << std::endl;


		//------------Olefine parenti da  r + o2 -> ole + ho2 ----------
		std::vector<double> velo_olepar_r(Temps.size(), 0.0);

		reaout << "Olefins parent from R + O2:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Olefins parent from R + O2:");
		int num2 = 0;
		for (int i = 1; i <= numr; i++)
		{
			int pos_rad = pos_r(i);
			int posalfa[SIZEMAX + 1];
			r[i].scorri(pos_rad, 1, posalfa);
			for (int j = 1; j <= SIZEMAX; j++) if (posalfa[j] == 1 && r[i].numAbstractableH(j) != 0)
			{                            // j scorre le pos alfa
				Molecola ole = HC;
				ole.conc = 0;
				int isPossible = ole.addole(pos_rad, j);

				int a = Add(ole, olepar_r, &num2, max2);
				Add(ole, &olefins);

				std::string reacomm = k.v_ole_par_r(r[i].numAbstractableH(j));
				addReactionToList(&reactionsList, std::vector<Molecola>{ r[i], O2 }, std::vector<Molecola>{ olepar_r[a], HO2 }, new double[3] { k.A, k.n, k.E }, "Olefins parent from R + O2:", reacomm);
				k.wrirea(reaout, 5, k.A, k.n, k.E, r[i], olepar_r[a]);
				//chemout.wrireaDetailed(5, k.A, k.n, k.E, r[i], olepar_r[a]);
			};
		};
		reaout << endl;
		std::cout << "Olefins from R + O2 added." << std::endl;


		//------------Decomposition gamma-qooh---------------------------------
		Molecola* allProdBetaQooh;
		allProdBetaQooh = new Molecola[2 * max3];
		int numOfProdBetaQooh = 0;

		reaout << "Decomposition gamma-QOOH and delta-QOOH:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Decomposition gamma-QOOH and delta-QOOH:");
		int num30 = 0, num31 = 0;
		for (int i = 1; i <= numqooh; i++)
		{
			Pos p = pos_qooh(i);
			int dist = qooh[i].dist(p.rad, p.ooh);
			if (dist == 2)
			{
				//                  trova la posizione in cui spezzare
				int pos_alfa;
				{
					int alfa_ooh[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					qooh[i].scorri(p.ooh, 1, alfa_ooh);
					qooh[i].scorri(p.rad, 1, alfa_rad);
					int j;
					for (j = 1; j <= SIZEMAX; j++)
						if (alfa_ooh[j] == 1 && alfa_rad[j] == 1) break;
					pos_alfa = j;
				};
				//                         f      crea le molecole spezzate
				Molecola m1, m2;
				int isPossible = qooh[i].spezza(p.ooh, pos_alfa, &m1, &m2);
				if (isPossible == 0) continue;                     // MODIFICATION
				p.ooh = m1.trova(2);
				m1.tipo(p.ooh, 1);
				m1.addcheto(p.ooh);

				int a = Add(m1, aldbetaqooh, &num30, max3);
				int b = Add(m2, olebetaqooh, &num31, max3);
				int c = Add(m1, allProdBetaQooh, &numOfProdBetaQooh, 2 * max3);
				int d = Add(m2, allProdBetaQooh, &numOfProdBetaQooh, 2 * max3);

				std::string reacomm = k.v_ole_from_gamma_qooh(qooh[i].tipoROOH(qooh[i].trova(4)), qooh[i].tipoR(qooh[i].trova(2)));
				addReactionToList(&reactionsList, qooh[i], std::vector<Molecola>{aldbetaqooh[a], olebetaqooh[b], OH}, new double[3] { k.A, k.n, k.E }, "Decomposition gamma-QOOH and delta-QOOH:", reacomm);
				k.wrirea(reaout, 8, k.A, k.n, k.E, qooh[i], aldbetaqooh[a], olebetaqooh[b]);
				//chemout.wrireaDetailed(8, k.A, k.n, k.E, qooh[i], aldbetaqooh[a], olebetaqooh[b]);
			};
		};
		reaout << endl;
		//std::cout << "Decomposition gamma-QOOH added." << std::endl;

		//------------Decomposition delta-qooh---------------------------------
		Molecola* allProdBetaQooh2;
		allProdBetaQooh2 = new Molecola[2 * max3];
		int numOfProdBetaQooh2 = 0;

		//reaout << "Decomposition delta-QOOH:" << endl;
		//chemout.endlineDetailed();
		//chemout.writeDetReacName("Decomposition delta-QOOH:");

		for (int i = 1; i <= numqooh; i++)
		{
			Pos p = pos_qooh(i);
			int dist = qooh[i].dist(p.rad, p.ooh);
			if (dist == 3)
			{
				//                 find the position where to break
				int pos_break1;
				{
					int beta_ooh[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					qooh[i].scorri(p.ooh, 2, beta_ooh);
					qooh[i].scorri(p.rad, 1, alfa_rad);
					int j;
					for (j = 1; j <= SIZEMAX; j++)
						if (beta_ooh[j] == 1 && alfa_rad[j] == 1) break;
					pos_break1 = j;
				};
				int pos_break2;
				{
					int alfa_ooh[SIZEMAX + 1];
					int beta_rad[SIZEMAX + 1];
					qooh[i].scorri(p.ooh, 1, alfa_ooh);
					qooh[i].scorri(p.rad, 2, beta_rad);
					int j;
					for (j = 1; j <= SIZEMAX; j++)
						if (alfa_ooh[j] == 1 && beta_rad[j] == 1) break;
					pos_break2 = j;
				};
				//                               crea le molecole spezzate
				Molecola m1, m2;
				int isPossible = qooh[i].spezza(pos_break1, pos_break2, &m1, &m2);
				if (isPossible == 0) continue;                     // MODIFICATION

				int c = Add(m1, allProdBetaQooh2, &numOfProdBetaQooh2, 2 * max3);
				int d = Add(m2, allProdBetaQooh2, &numOfProdBetaQooh2, 2 * max3);

				std::string reacomm = k.v_ole_from_delta_qooh(qooh[i].tipoROOH(qooh[i].trova(4)), qooh[i].tipoR(qooh[i].trova(2)));
				addReactionToList(&reactionsList, qooh[i], std::vector<Molecola>{m1, m2}, new double[3] { k.A, k.n, k.E }, "Decomposition gamma-QOOH and delta-QOOH:", reacomm);
				k.wrirea(reaout, 2, k.A, k.n, k.E, qooh[i], m1, m2);
				//chemout.wrireaDetailed(8, k.A, k.n, k.E, qooh[i], aldbetaqooh[a], olebetaqooh[b]);
			};
		};
		reaout << endl;
		std::cout << "Decomposition gamma-QOOH and delta-QOOH added." << std::endl;

		//------------Etero da     qooh -> ete + oh---------------------
		std::vector<double> velo_ete_qooh(Temps.size(), 0.0);
		reaout << "Ethers from QOOH:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Ethers from QOOH:");
		std::vector<Molecola> cycEthers;
		int num4 = 0;
		for (int i = 1; i <= numqooh; i++)
		{
			Pos p = pos_qooh(i);
			int dist = qooh[i].dist(p.rad, p.ooh);
			AnelloO a;
			
			if (dist == 1) a = ao3;
			else if (dist == 2) a = ao4;
			else if (dist == 3) a = ao5;
			else if (dist == 4) a = ao6;
			else continue;
			if (qooh[i].tipoROOH(qooh[i].trova(4)) == Rp && qooh[i].tipoR(qooh[i].trova(2)) == Rp && a == ao3) continue;
			Molecola e = HC;
			e.conc = 0;
			e.addetero(p.rad, p.ooh);

			//int aa = Add(e, etero, &num4, max4);
			Add(e, &cycEthers);

			// check if tertiary or quaternary carbons are in the ring
			//int trovati[SIZEMAX + 1];
			//roo[i].scorri(pos_o2, 1, trovati);	// find the carbon at distance 1 from the carbon with the OO
			int found1[SIZEMAX + 1];
			int found2[SIZEMAX + 1];
			qooh[i].scorri(p.rad, dist, found1);
			qooh[i].scorri(p.ooh, dist, found2);
			Carbonio corr = Cs;
			for (int ind = 0; ind < SIZEMAX + 1; ind++)
			{
				if (found1[ind] == 1 && found2[ind] == 1 && qooh[i].tipoC(ind) == Ct)
					corr = Ct;
				if (found1[ind] == 1 && found2[ind] == 1 && qooh[i].tipoC(ind) == Cq)
					corr = Cq;
			}

			std::string reacomm = k.v_ether_from_qooh(qooh[i].tipoROOH(qooh[i].trova(4)), qooh[i].tipoR(qooh[i].trova(2)), a, corr);
			addReactionToList(&reactionsList, qooh[i], std::vector<Molecola>{ e, OH }, new double[3] { k.A, k.n, k.E }, "Ethers from QOOH:", reacomm);
			k.wrirea(reaout, 7, k.A, k.n, k.E, qooh[i], e);
			//chemout.wrireaDetailed(7, k.A, k.n, k.E, qooh[i], etero[aa]);
		};
		reaout << endl;
		std::cout << "Ethers from QOOH added." << std::endl;


		//------------Etero da     P(OOH)2 -> ete-OOH + oh---------------------
		std::vector<double> velo_ete_pooh2(Temps.size(), 0.0);
		reaout << "Ethers-OOH from P(OOH)2:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Ethers-OOH from P(OOH)2:");
		std::vector<Molecola> etherOOH;
		for (int i = 0; i < POOH2_isom_reac.size(); i++)
		{
			AnelloO a;
			int dist = 1;
			switch (POOH2_isom_ring[i])
			{
			case a5:
				a = ao3;
				dist = 1;
				break;
			case a6:
				a = ao4;
				dist = 2;
				break;
			case a7:
				a = ao5;
				dist = 3;
				break;
			case a8:
				a = ao6;
				dist = 4;
				break;
			}
			if (POOH2_isom_reac[i].tipoROOH(POOH2_isom_OOH_pos[i]) == Rp && POOH2_isom_reac[i].tipoR(POOH2_isom_rad_pos[i]) == Rp && a == ao3) continue;
			Molecola mol = POOH2_isom_reac[i];
			mol.removeOOH(POOH2_isom_OOH_pos[i]);
			mol.Crad_to_C(POOH2_isom_rad_pos[i]);
			mol.addetero(POOH2_isom_OOH_pos[i], POOH2_isom_rad_pos[i]);
			Add(mol, &etherOOH);

			int found1[SIZEMAX + 1];
			int found2[SIZEMAX + 1];
			POOH2_isom_reac[i].scorri(POOH2_isom_rad_pos[i], dist, found1);
			POOH2_isom_reac[i].scorri(POOH2_isom_OOH_pos[i], dist, found2);
			Carbonio corr = Cs;
			for (int ind = 0; ind < SIZEMAX + 1; ind++)
			{
				if (found1[ind] == 1 && found2[ind] == 1 && POOH2_isom_reac[i].tipoC(ind) == Ct)
					corr = Ct;
				if (found1[ind] == 1 && found2[ind] == 1 && POOH2_isom_reac[i].tipoC(ind) == Cq)
					corr = Cq;
			}

			std::string reacomm = k.v_ether_from_pooh2(POOH2_isom_reac[i].tipoROOH(POOH2_isom_OOH_pos[i]), POOH2_isom_reac[i].tipoR(POOH2_isom_rad_pos[i]), a, corr);
			addReactionToList(&reactionsList, POOH2_isom_reac[i], std::vector<Molecola>{ mol, OH }, new double[3] { k.A, k.n, k.E }, "Ethers-OOH from P(OOH)2:", reacomm);
			k.wrirea(reaout, 7, k.A, k.n, k.E, POOH2_isom_reac[i], mol);
		}
		reaout << endl;
		std::cout << "Ethers-OOH from P(OOH)2 added." << std::endl;

		//------------- P(OOH)2 decomposition 1: beta-P(OOH)2 -> olefin + HO2 ---------------
		reaout << "P(OOH)2 decomposition 1:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("P(OOH)2 decomposition 1:");
		for (int i = 0; i < POOH2_isom_reac.size(); i++)
		{
			if (POOH2_isom_reac[i].dist(POOH2_isom_OOH_pos[i], POOH2_isom_rad_pos[i]) == 1)
			{
				if (POOH2_isom_reac[i].tipoROOH(POOH2_isom_OOH_pos[i]) == Rp && POOH2_isom_reac[i].tipoR(POOH2_isom_rad_pos[i]) == Rp) continue; // CHECK if it is right to neglect those reactions
				Molecola mol = POOH2_isom_reac[i];
				mol.removeOOH(POOH2_isom_OOH_pos[i]);
				mol.Crad_to_C(POOH2_isom_rad_pos[i]);
				mol.addole(POOH2_isom_OOH_pos[i], POOH2_isom_rad_pos[i]);
				std::string reacomm = k.v_pooh2_dec_1(POOH2_isom_reac[i].tipoROOH(POOH2_isom_OOH_pos[i]), POOH2_isom_reac[i].tipoR(POOH2_isom_rad_pos[i]));
				addReactionToList(&reactionsList, POOH2_isom_reac[i], std::vector<Molecola>{ mol, HO2 }, new double[3] { k.A, k.n, k.E }, "P(OOH)2 decomposition 1:", reacomm);
				k.wrirea(reaout, 6, k.A, k.n, k.E, POOH2_isom_reac[i], mol);
			}
		}
		reaout << endl;
		std::cout << "P(OOH)2 decomposition 1 added" << std::endl;

		//------------- P(OOH)2 decomposition 2: gamma-P(OOH)2 -> olefin + RO + OH ---------------
		reaout << "P(OOH)2 decomposition 2:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("P(OOH)2 decomposition 2:");

		for (int i = 0; i < POOH2_isom_reac.size(); i++)
		{
			if (POOH2_isom_reac[i].dist(POOH2_isom_OOH_pos[i], POOH2_isom_rad_pos[i]) == 2)
			{
				//                  find the position where to break
				int pos_alfa;
				{
					int alfa_ooh[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					POOH2_isom_reac[i].scorri(POOH2_isom_OOH_pos[i], 1, alfa_ooh);
					POOH2_isom_reac[i].scorri(POOH2_isom_rad_pos[i], 1, alfa_rad);
					int j;
					for (j = 1; j <= SIZEMAX; j++)
						if (alfa_ooh[j] == 1 && alfa_rad[j] == 1) break;
					pos_alfa = j;
				};
				//                         f      generate the broken molecules
				if (POOH2_isom_reac[i].tipo(pos_alfa) == 4)
					continue;
				Molecola m1, m2;
				int isPossible = POOH2_isom_reac[i].spezza(POOH2_isom_OOH_pos[i], pos_alfa, &m1, &m2);
				if (isPossible == 0) continue;



				int pos = m1.trova(2);
				m1.tipo(pos, 1);
				m1.addcheto(pos);
				std::string reacomm = k.v_pooh2_dec_2(POOH2_isom_reac[i].tipoROOH(POOH2_isom_OOH_pos[i]), POOH2_isom_reac[i].tipoR(POOH2_isom_rad_pos[i]));
				addReactionToList(&reactionsList, POOH2_isom_reac[i], std::vector<Molecola>{ m1, m2, OH }, new double[3] { k.A, k.n, k.E }, "P(OOH)2 decomposition 2:", reacomm);
				k.wrirea(reaout, 8, k.A, k.n, k.E, POOH2_isom_reac[i], m1, m2);
			}
		}
		reaout << endl;
		std::cout << "P(OOH)2 decomposition 2 added" << std::endl;

		//------------- P(OOH)2 decomposition 3: delta-P(OOH)2 -> olefin + beta-P(OOH)2 ---------------
		reaout << "P(OOH)2 decomposition 3:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("P(OOH)2 decomposition 3:");
		for (int i = 0; i < POOH2_isom_reac.size(); i++)
		{
			if (POOH2_isom_reac[i].dist(POOH2_isom_OOH_pos[i], POOH2_isom_rad_pos[i]) == 3)
			{
				int pos_rad;		// find the position where the radical will be located
				{
					int beta_ooh[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					POOH2_isom_reac[i].scorri(POOH2_isom_OOH_pos[i], 1, beta_ooh);
					POOH2_isom_reac[i].scorri(POOH2_isom_rad_pos[i], 2, alfa_rad);
					int j;
					for (j = 1; j <= SIZEMAX; j++)
						if (beta_ooh[j] == 1 && alfa_rad[j] == 1) break;
					pos_rad = j;
				};
				if (POOH2_isom_reac[i].tipo(pos_rad) == 4)  // if there is a COOH in this position continue
					continue;
				//                  find the position where to break
				int pos_break1;
				{
					int beta_ooh[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					POOH2_isom_reac[i].scorri(POOH2_isom_OOH_pos[i], 2, beta_ooh);
					POOH2_isom_reac[i].scorri(POOH2_isom_rad_pos[i], 1, alfa_rad);
					int j;
					for (j = 1; j <= SIZEMAX; j++)
						if (beta_ooh[j] == 1 && alfa_rad[j] == 1) break;
					pos_break1 = j;
				};
				int pos_break2;
				{
					int alfa_ooh[SIZEMAX + 1];
					int beta_rad[SIZEMAX + 1];
					POOH2_isom_reac[i].scorri(POOH2_isom_OOH_pos[i], 1, alfa_ooh);
					POOH2_isom_reac[i].scorri(POOH2_isom_rad_pos[i], 2, beta_rad);
					int j;
					for (j = 1; j <= SIZEMAX; j++)
						if (alfa_ooh[j] == 1 && beta_rad[j] == 1) break;
					pos_break2 = j;
				};
				if (POOH2_isom_reac[i].tipo(pos_break1) == 4)
					continue;

				//                         f      crea le molecole spezzate
				Molecola m1, m2;
				int isPossible = POOH2_isom_reac[i].spezza(pos_break1, pos_break2, &m1, &m2);
				if (isPossible == 0) continue;

				std::string reacomm = k.v_pooh2_dec_3(POOH2_isom_reac[i].tipoROOH(POOH2_isom_OOH_pos[i]), POOH2_isom_reac[i].tipoR(POOH2_isom_rad_pos[i]));
				addReactionToList(&reactionsList, POOH2_isom_reac[i], std::vector<Molecola>{ m1, m2}, new double[3] { k.A, k.n, k.E }, "P(OOH)2 decomposition 3:", reacomm);
				k.wrirea(reaout, 2, k.A, k.n, k.E, POOH2_isom_reac[i], m1, m2);
			}
		}
		reaout << endl;
		std::cout << "P(OOH)2 decomposition 3 added" << std::endl;

		//------------- ether OOH decomposition ---------------
		reaout << "Ether-OOH decomposition:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Ether-OOH decomposition:");

		for (int i = 0; i < etherOOH.size(); i++)
		{
			int alfa_ooh[SIZEMAX + 1];					// array in which there is 1 if the atom in that position is in alfta to OOH
			etherOOH[i].scorri(etherOOH[i].trova(4), 1, alfa_ooh);
			for (int j = 1; j <= SIZEMAX; j++)
			{
				if (alfa_ooh[j] == 1)
				{
					Molecola m1;
					Molecola m2;
					//int isPossible = oleOOH[i].homoBreak(oleOOH[i].trova(4), j, &m1, &m2);
					int isPossible = etherOOH[i].spezza(etherOOH[i].trova(4), j, &m1, &m2);
					if (isPossible == 0) continue;
					//std::cout << m1 << std::endl << m2 << std::endl << std::endl;
					std::string reacomm = k.v_etherooh_dec(etherOOH[i].tipoROOH(etherOOH[i].trova(4)));
					
					std::vector<Molecola> prods;
					//std::cout << etherOOH[i] << std::endl;
					
					
					if (m2.size() == 0)    // if the ring was opened
					{
						//std::cout << "1 - ";
						m1.Crad_to_keto(etherOOH[i].trova(4));
						prods = decomposeLinEthRO(m1);
					}
					else if (m2.kindOfSPecies() == cEthR_) //do not consider reactions that leads to cyclic ethers CO
					{
						//std::cout << "2 - ";
						m1.Crad_to_keto(m1.trova(2));
						prods = decomposeCEthR(m2);
						prods.push_back(m1);
					}
					else
					{
						continue;
					}
					
					prods = fullyDecomposeRO(prods);

					//std::cout << "TEST TEST TEST" << std::endl;
					//std::cout << etherOOH[i] << std::endl;
					//std::cout << "     " << m1 << std::endl;
					//std::cout << "     " << m2 << std::endl;
					//std::cout << "----------------------------------" << std::endl;
					//for (int asd = 0; asd < prods.size(); asd++)
					//	std::cout << "       " << prods[asd] << std::endl;


					switch (prods.size())
					{
					case 2:
						addReactionToList(&reactionsList, etherOOH[i], std::vector<Molecola>{ prods[0], prods[1], OH }, new double[3] { k.A, k.n, k.E }, "Ether-OOH decomposition:", reacomm);
						k.wrirea(reaout, 8, k.A, k.n, k.E, etherOOH[i], prods[0], prods[1]);
						break;
					case 3:
						addReactionToList(&reactionsList, etherOOH[i], std::vector<Molecola>{ prods[0], prods[1], prods[2], OH }, new double[3] { k.A, k.n, k.E }, "Ether-OOH decomposition:", reacomm);
						k.wrirea(reaout, 12, k.A, k.n, k.E, etherOOH[i], prods[0], prods[1], prods[2]);
						break;
					case 4:
						addReactionToList(&reactionsList, etherOOH[i], std::vector<Molecola>{ prods[0], prods[1], prods[2], prods[3], OH }, new double[3] { k.A, k.n, k.E }, "Ether-OOH decomposition:", reacomm);
						k.wrirea(reaout, 24, k.A, k.n, k.E, etherOOH[i], prods[0], prods[1], prods[2], prods[3]);
						break;
					case 5:
						addReactionToList(&reactionsList, etherOOH[i], std::vector<Molecola>{ prods[0], prods[1], prods[2], prods[3], prods[4], OH }, new double[3] { k.A, k.n, k.E }, "Ether-OOH decomposition:", reacomm);
						k.wrirea(reaout, 27, k.A, k.n, k.E, etherOOH[i], prods[0], prods[1], prods[2], prods[3], prods[4]);
						break;
					default:
						std::cerr << "WARNING: in decomposition of CEthR there is an unexpected number of products! Reaction skipped." << std::endl;
					}
					std::cout << std::endl << std::endl;
				}
			}
		}
		reaout << endl;
		std::cout << "Ether-OOH decomposition added" << std::endl;

		//------------Olefine parenti da   qooh -> ole + ho2 -----------
		std::vector<double> velo_ole_par_qooh(Temps.size(), 0.0);
		reaout << "Decomposition beta-QOOH:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Decomposition beta-QOOH:");
		int num5 = 0;
		for (int i = 1; i <= numqooh; i++)
		{
			Pos p = pos_qooh(i);
			int dist = qooh[i].dist(p.rad, p.ooh);
			if (dist == 1)
			{
				if (qooh[i].tipoR(qooh[i].trova(2)) == Rp && qooh[i].tipoROOH(qooh[i].trova(4)) == Rp) continue;				// CHECK if it is right to skip those reactions
				Molecola ole = HC;
				ole.conc = 0;
				ole.addole(p.rad, p.ooh);

				int a = Add(ole, olepar_qooh, &num5, max5);
				Add(ole, &olefins);
				std::string reacomm = k.v_ole_from_beta_qooh(qooh[i].tipoROOH(qooh[i].trova(4)), qooh[i].tipoR(qooh[i].trova(2)));
				addReactionToList(&reactionsList, qooh[i], std::vector<Molecola>{ olepar_qooh[a], HO2}, new double[3] { k.A, k.n, k.E }, "Decomposition beta-QOOH:", reacomm);
				k.wrirea(reaout, 6, k.A, k.n, k.E, qooh[i], olepar_qooh[a]);
				//chemout.wrireaDetailed(6, k.A, k.n, k.E, qooh[i], olepar_qooh[a]);
			};
		};
		reaout << endl;
		std::cout << "Decomposition beta-QOOH added." << std::endl;

		//------------Olefine from OOQOOH -> ole-OOH + HO2 -----------
		reaout << "Olefins-OOH from OOQOOH:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Olefins-OOH from OOQOOH:");

		std::vector<Molecola> oleOOH;
		for (int i = 1; i <= numqooh; i++)
		{
			int j;
			int pos_o2 = ooqooh[i].trova(3);
			int trovati[SIZEMAX + 1];

			ooqooh[i].scorri(pos_o2, 1, trovati);	// find the carbon at distance 1 from the carbon with the OO
			for (j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
				if (trovati[j] == 1 && ooqooh[i].numAbstractableH(j) != 0 && ooqooh[i].trova(4) != j)
				{
					if (ooqooh[i].tipoROO(ooqooh[i].trova(3)) == Rp && ooqooh[i].tipoC(j) == Cp) continue;   // CHECH if it is correct to neglect those reactions
					Molecola ole = ooqooh[i];
					ole.removeOO(ooqooh[i].trova(3));
					ole.addole(ooqooh[i].trova(3), j);
					Add(ole, &oleOOH);
					std::string reacomm = k.v_ole_par_ooqooh(ooqooh[i].tipoROO(ooqooh[i].trova(3)), ooqooh[i].tipoC(j), ooqooh[i].numAbstractableH(j));
					addReactionToList(&reactionsList, ooqooh[i], std::vector<Molecola>{ ole, HO2 }, new double[3] { k.A, k.n, k.E }, "Olefins-OOH from OOQOOH:", reacomm);
					k.wrirea(reaout, 6, k.A, k.n, k.E, ooqooh[i], ole);
					//chemout.wrireaDetailed(1, k.A, k.n, k.E, roo[i], qooh[tipoqooh[pos_ooh][pos_r]]);
				};

		};  // end for i
		reaout << endl;
		std::cout << "Olefins-OOH from OOQOOH added." << std::endl;

		//------------- ole OOH decomposition ---------------
		reaout << "Ole-OOH decomposition:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Ole-OOH decomposition:");

		for (int i = 0; i < oleOOH.size(); i++)
		{
			int alfa_ooh[SIZEMAX + 1];					// array in which there is 1 if the atom in that position is in alfta to OOH
			oleOOH[i].scorri(oleOOH[i].trova(4), 1, alfa_ooh);
			for (int j = 1; j <= SIZEMAX; j++)
			{
				if (alfa_ooh[j] == 1)
				{
					Molecola m1;
					Molecola m2;
					//int isPossible = oleOOH[i].homoBreak(oleOOH[i].trova(4), j, &m1, &m2);
					int isPossible = oleOOH[i].spezza(oleOOH[i].trova(4), j, &m1, &m2);
					if (isPossible == 0) continue;
					m1.Crad_to_keto(m1.trova(2));
					std::string reacomm = k.v_oleooh_dec(oleOOH[i].tipoROOH(oleOOH[i].trova(4)));
					addReactionToList(&reactionsList, oleOOH[i], std::vector<Molecola>{ m1, m2, OH }, new double[3] { k.A, k.n, k.E }, "Ole-OOH decomposition:", reacomm);
					k.wrirea(reaout, 8, k.A, k.n, k.E, oleOOH[i], m1, m2);
				}
			}
		}
		reaout << endl;
		std::cout << "Ole-OOH decomposition added" << std::endl;

		//------------Cheto da ooqooh -> kehy + oh----------------------
		reaout << "Conversion of OOQOOH to ketohydroperoxide:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Conversion of OOQOOH to ketohydroperoxides:");
		std::vector<double> velo_ooqooh_to_khp(Temps.size(), 0.0);
		Molecola* khp;
		khp = new Molecola[HC.size() * 6 + 1]; if (!khp) exit(203);
		int numkhp = 0;
		for (int i = 1; i <= numqooh; i++)
		{
			int pos_oo = ooqooh[i].trova(3);		// find the position of the oo group
			int pos_ooh = ooqooh[i].trova(4);		// find the position of the ooh group

			int dist = ooqooh[i].dist(pos_oo, pos_ooh);		// find the distance between the oo and ooh group
			// from the distance find the type of ring formed during the reaction
			Anello a;
			if (dist == 1) a = a5;
			else if (dist == 2) a = a6;
			else if (dist == 3) a = a7;
			else if (dist == 4) a = a8;
			else continue;					// if the ring is bigger than 8 atoms skip (the reaction could not happen)

			Molecola temp;
			temp = ooqooh[i];
			temp.tipo(pos_oo, 4);						// change the oo group in an ooh group
			temp.tipo(pos_ooh, 1);						// remove the ooh group and ...
			int isPossible = temp.addcheto(pos_ooh);	// ... replace it with the group =o Check if this step is possible.
			if (isPossible != 0)						// if it is possible save the group
			{
				khp[numkhp + 1] = temp;
				std::string reacomm = k.v_ooqooh_to_khp(ooqooh[i].tipoROO(ooqooh[i].trova(3)), ooqooh[i].tipoROOH(ooqooh[i].trova(4)), a, ooqooh[i].numAbstractableH(ooqooh[i].trova(4)));
				addReactionToList(&reactionsList, ooqooh[i], std::vector<Molecola>{ temp, OH }, new double[3] { k.A, k.n, k.E }, "Conversion of OOQOOH to ketohydroperoxide:", reacomm);
				k.wrirea(reaout, 7, k.A, k.n, k.E, ooqooh[i], temp);
				//chemout.wrireaDetailed(7, k.A, k.n, k.E, ooqooh[i], temp);
				numkhp++;
				//std::cout << khp[i] << std::endl;
			}
		}
		reaout << endl;
		std::cout << "Ketohydroperoxides formation added." << std::endl;

		//for (int i = 1; i <= numkhp; i++)
		//	std::cout << khp[i] << std::endl;

		reaout << "Decomposition of ketohydroperoxides:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Decomposition of ketohydroperoxides:");
		//int num60 = 0, num61 = 0;
		Molecola* prod_khp_dec;
		prod_khp_dec = new Molecola[2 * HC.numberOfC()];
		int number_of_khp_dec_prod = 0;

		for (int i = 1; i <= numkhp; i++)
		{
			int pos_ooh = khp[i].trova(4);		// find the position of the ooh group
			int pos_o = khp[i].trova(7);		// find the position of the =o group
			int dist = khp[i].dist(pos_o, pos_ooh);		// find the distance between the ooh and =o group

			int pos_alfa;
			{
				int alfa_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				khp[i].scorri(pos_o, dist - 1, alfa_ooh);
				khp[i].scorri(pos_ooh, 1, alfa_rad);
				int j;
				for (j = 1; j <= SIZEMAX; j++)
					if (alfa_ooh[j] == 1 && alfa_rad[j] == 1) break;
				pos_alfa = j;
			};

			Molecola m1, m2;

			int isPossible = khp[i].spezza(pos_alfa, pos_ooh, &m1, &m2);
			if (isPossible == 0) continue;                     // MODIFICATION
			int pos_rad_m2 = m2.trova(2);			// find the radical on the molecule m2
			m2.tipo(pos_rad_m2, 1);					// remove the radical and..
			m2.addcheto(pos_rad_m2);				// ... replace it with a cheto group
			int posRO = AddIsomer(m1, prod_khp_dec, &number_of_khp_dec_prod, 2 * HC.numberOfC());
			int posCHO = AddIsomer(m2, prod_khp_dec, &number_of_khp_dec_prod, 2 * HC.numberOfC());

			std::vector<Molecola> RO_dec_prod;
			bool thereIsRO = false;
			//std::cout << std::endl << "DEBUG: " << std::endl;
			//std::cout << "Original reaction:" << std::endl;
			//std::cout << khp[i] << "  =>  " << m1 << "  +  " << m2 << std::endl;

			if (m1.kindOfSPecies() == RO_)
			{
				//std::cout << m1 << " is a RO and decomposed into: ";
				RO_dec_prod = fullyDecomposeRO(m1);
				thereIsRO = true;

			}
			//std::cout << std::endl;
			//std::cout << "distance = " << dist << std::endl;
			std::string reacomm = k.v_khp_decomp(khp[i].tipoROOH(khp[i].trova(4)), dist);

			if (thereIsRO == true && RO_dec_prod.size() == 2)
			{
				addReactionToList(&reactionsList, khp[i], std::vector<Molecola>{ RO_dec_prod[0], RO_dec_prod[1], m2, OH }, new double[3] { k.A, k.n, k.E }, "Decomposition of ketohydroperoxides:", reacomm);
				k.wrirea(reaout, 12, k.A, k.n, k.E, khp[i], RO_dec_prod[0], RO_dec_prod[1], m2);
			}
			else if (thereIsRO == true && RO_dec_prod.size() == 3)
			{
				addReactionToList(&reactionsList, khp[i], std::vector<Molecola>{ RO_dec_prod[0], RO_dec_prod[1], RO_dec_prod[2], m2, OH }, new double[3] { k.A, k.n, k.E }, "Decomposition of ketohydroperoxides:", reacomm);
				k.wrirea(reaout, 24, k.A, k.n, k.E, khp[i], RO_dec_prod[0], RO_dec_prod[1], RO_dec_prod[2], m2);
			}
			else if (thereIsRO == true && RO_dec_prod.size() == 4)
			{
				addReactionToList(&reactionsList, khp[i], std::vector<Molecola>{ RO_dec_prod[0], RO_dec_prod[1], RO_dec_prod[2], RO_dec_prod[3], m2, OH }, new double[3] { k.A, k.n, k.E }, "Decomposition of ketohydroperoxides:", reacomm);
				k.wrirea(reaout, 24, k.A, k.n, k.E, khp[i], RO_dec_prod[0], RO_dec_prod[1], RO_dec_prod[2], RO_dec_prod[3], m2);
			}
			else
			{
				addReactionToList(&reactionsList, khp[i], std::vector<Molecola>{ m1, m2, OH }, new double[3] { k.A, k.n, k.E }, "Decomposition of ketohydroperoxides:", reacomm);
				k.wrirea(reaout, 8, k.A, k.n, k.E, khp[i], m1, m2);
			}
		}
		reaout << endl;
		std::cout << "Ketohydroperoxides decomposition added." << std::endl;


		// CYCLIC ETHERS DECOMPOSITION
		reaout << "Cyclic ethers decomposition:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Cyclic ethers decomposition:");

		for (int i = 0; i < cycEthers.size(); i++)								// iterate through all the molecules in the cyclyc ethers list
		{
			std::vector<int> etherPos = cycEthers[i].posEthero();
			int dist = cycEthers[i].dist(etherPos[0], etherPos[1]);
			for (int j = 0; j < 2; j++)			// iterates through the two carbon the oxygen is bonded with
			{
				int posRad = etherPos[j];
				int posCheto = etherPos[1 - j];
				int numH = cycEthers[i].numH(posCheto);
				int numH2 = cycEthers[i].numH(posRad);

				if (numH == 0 && numH2 == 0) // if there are no abstractable hydrogens in alpha to the oxygen use the secondary path
				{
					// look for the hydrogens to abstract
					// if there are methyls abstract it from them
					Molecola interm = cycEthers[i];
					int alfa_rad[SIZEMAX + 1];
					interm.scorri(posRad, 1, alfa_rad);
					bool normalPathTaken = false;
					int numAbsH = 0;
					int posAbs = 0;
					Carbonio carbType = Cp;
					for (int l = 1; l <= SIZEMAX; l++)
					{
						if (alfa_rad[l] == 1 && !interm.isGroup(l) && interm.tipoC(l) == Cp)
						{
							numAbsH += 3;
							posAbs = l;			// no matter if it is overwritten, we need only one position
						}
					}
					if (posAbs == 0)	// if no methyl was found look for a secondary carbon in alpha
					{
						carbType = Cs;
						for (int l = 1; l <= SIZEMAX; l++)
						{
							if (alfa_rad[l] == 1 && !interm.isGroup(l) && interm.tipoC(l) == Cs)
							{
								numAbsH += 2;
								posAbs = l;			// no matter if it is overwritten, we need only one position
							}
						}
					}

					interm.removeCycEther(posRad, posCheto);
					interm.addole(posAbs, posRad);
					//interm.addcheto(posCheto);

					// find the carbon in alpha to the cheto more close to the double bond
					int alfa_cheto[SIZEMAX + 1];
					interm.scorri(posCheto, 1, alfa_cheto);
					int posAlfaCheto = 0;
					int distAlfaCheto = SIZEMAX;
					for (int l = 1; l <= SIZEMAX; l++)
					{
						if (alfa_cheto[l] == 1)
						{
							if (interm.dist(posAbs, l) < distAlfaCheto)
							{
								distAlfaCheto = interm.dist(posAbs, l);
								posAlfaCheto = l;
							}
						}
					}
					Molecola m1, m2;
					interm.spezza(posCheto, posAlfaCheto, &m1, &m2);
					int posRadM1 = m1.posCrad();
					m1.Crad_to_C(posRadM1);
					m1.addcheto(posRadM1);

					std::string correction = "none";
					correction = oh_abstraction_correction(&HC, posAbs);

					std::string reacomm = k.v_h_abstraction(oh, carbType, numAbsH, 1, correction);
					addReactionToList(&reactionsList, std::vector<Molecola>{ cycEthers[i], OH },
						std::vector<Molecola>{ m1, m2, H2O }, new double[3] { k.A, k.n, k.E }, "Cyclic ethers decomposition:", reacomm);
					k.wrirea(reaout, 23, k.A, k.n, k.E, cycEthers[i], m1, m2);
				}
				else if (numH == 0 && numH2 != 0)
				{
					continue;
				}
				else
				{
					Molecola interm = cycEthers[i];
					interm.removeCycEther(posRad, posCheto);
					interm.addcheto(posCheto);
					interm.C_to_Crad(posRad);
					bool reactionHappened = false;
					Molecola m1, m2;
					switch (dist)
					{
					case 1:				// 3 member ring
					{
						//std::cout << "Molecule " << cycEthers[i] << " is a 3 atom ring cyclic ether." << std::endl;
						// search if there is a carbon in alfa at the cheto and in beta at the radical
						int alfa_cheto[SIZEMAX + 1];
						int beta_rad[SIZEMAX + 1];
						interm.scorri(posCheto, 1, alfa_cheto);
						interm.scorri(posRad, 2, beta_rad);
						bool normalPathTaken = false;
						for (int l = 1; l <= SIZEMAX; l++)				// normal path: 
						{
							if (alfa_cheto[l] == 1 && beta_rad[l] == 1)
							{
								interm.spezza(posCheto, l, &m1, &m2);
								normalPathTaken = true;
								reactionHappened = true;
								break;
							}
						}
						if (!normalPathTaken)
						{
							for (int l = 1; l <= SIZEMAX; l++)
							{
								if (beta_rad[l] == 1)
								{
									for (int h = 1; h <= interm.size(); h++)
									{
										if (interm.dist(h, l) == 1 && interm.dist(h, posRad) == 1)
										{
											interm.spezza(h, l, &m1, &m2);
											reactionHappened = true;
											break;
										}
										if (reactionHappened)
											break;
									}
								}
							}
						}
					}
					break;
					case 2:				// 4 member ring
					{
						//std::cout << "Molecule " << cycEthers[i] << " is a 4 atom ring cyclic ether." << std::endl;
						int alfa_cheto[SIZEMAX + 1];
						int alfa_rad[SIZEMAX + 1];
						interm.scorri(posCheto, 1, alfa_cheto);
						interm.scorri(posRad, 1, alfa_rad);
						for (int l = 1; l <= SIZEMAX; l++)
						{
							if (alfa_cheto[l] == 1 && alfa_rad[l] == 1)
							{
								interm.spezza(posCheto, l, &m1, &m2);
								reactionHappened = true;
								break;
							}
						}
					}
					break;
					case 3:				// 5 member ring
					{
						//std::cout << "Molecule " << cycEthers[i] << " is a 5 atom ring cyclic ether." << std::endl;
						int alfa_cheto[SIZEMAX + 1];
						int beta_cheto[SIZEMAX + 1];
						int alfa_rad[SIZEMAX + 1];
						int beta_rad[SIZEMAX + 1];
						interm.scorri(posCheto, 1, alfa_cheto);
						interm.scorri(posCheto, 2, beta_cheto);
						interm.scorri(posRad, 1, alfa_rad);
						interm.scorri(posRad, 2, beta_rad);
						for (int l = 1; l <= SIZEMAX; l++)
						{
							if (alfa_cheto[l] == 1 && beta_rad[l] == 1)
							{
								for (int h = 1; h <= SIZEMAX; h++)
								{
									if (beta_cheto[h] == 1 && alfa_rad[h] == 1)
									{
										interm.spezza(h, l, &m1, &m2);
										reactionHappened = true;
										break;
									}
								}
							}
						}
					}
					break;
					default:
						std::cout << "Molecule " << cycEthers[i] << " has the carbon bonding with the oxygen at a distance of " << dist << ", this is not a valid cyclic ether!" << std::endl;
						break;
					}

					if (reactionHappened)
					{
						std::vector<Molecola> RO_dec_prod;
						bool ism1RO = false;
						bool ism2RO = false;
						if (m1.kindOfSPecies() == RO_)
						{
							RO_dec_prod = decomposeRO(m1);
							ism1RO = true;
						}

						if (m2.kindOfSPecies() == RO_)
						{
							RO_dec_prod = decomposeRO(m2);
							ism2RO = true;
						}

						if (m1.kindOfSPecies() == RO_ && m2.kindOfSPecies() == RO_)
						{
							std::cerr << "ERROR: in decomposition of cyclic ethers, two RO are formed, this should not happen!" << std::endl;
							std::cerr << "Press any key to continue..." << std::endl;
							std::string input_keyboard;
							std::cin >> input_keyboard;
						}

						std::string reacomm = k.v_cyc_eth_dec();

						if (ism1RO == true && RO_dec_prod.size() == 2)
						{
							addReactionToList(&reactionsList, std::vector<Molecola>{ cycEthers[i], OH },
								std::vector<Molecola>{ RO_dec_prod[0], RO_dec_prod[1], m2, H2O }, new double[3] { k.A, k.n, k.E }, "Cyclic ethers decomposition:", reacomm);
							k.wrirea(reaout, 25, k.A, k.n, k.E, cycEthers[i], RO_dec_prod[0], RO_dec_prod[1], m2);
						}
						else if (ism2RO == true && RO_dec_prod.size() == 2)
						{

							addReactionToList(&reactionsList, std::vector<Molecola>{ cycEthers[i], OH },
								std::vector<Molecola>{ m1, RO_dec_prod[0], RO_dec_prod[1], H2O }, new double[3] { k.A, k.n, k.E }, "Cyclic ethers decomposition:", reacomm);
							k.wrirea(reaout, 25, k.A, k.n, k.E, cycEthers[i], m1, RO_dec_prod[0], RO_dec_prod[1]);
						}
						else
						{

							addReactionToList(&reactionsList, std::vector<Molecola>{ cycEthers[i], OH },
								std::vector<Molecola>{ m1, m2, H2O }, new double[3] { k.A, k.n, k.E }, "Cyclic ethers decomposition:", reacomm);
							k.wrirea(reaout, 23, k.A, k.n, k.E, cycEthers[i], m1, m2);
						}
					}
				}
			}



		}
		reaout << endl;
		std::cout << "Cyclic ethers decomposition added." << std::endl;


		// OLE TO ALLYLIC RADICAL: olefin + OH -> allylic-radical + H2O
		reaout << "Allylic radical formation:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Allylic radical formation:");
		std::vector<Molecola> allylicRads;
		for (int i = 0; i < olefins.size(); i++)								// iterate through all the molecules in the olefins list
		{
			for (int j = 1; j <= olefins[i].numberOle(); j++)					// iterate through all the double bond in the molecule
			{
				std::vector<int> oleCarbons = olefins[i].posOle(j);
				for (int l = 0; l < 2; l++)										// iterate trhough all (two) carbons involved in the double bond
				{
					int posalfa[SIZEMAX + 1];
					olefins[i].scorri(oleCarbons[l], 1, posalfa);
					for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)		// iterate through all the carbon in alfa at the carbon involved in the double bond
					{
						if (olefins[i].tipo(h) != 1) continue;		// the H can be abstracted only on normal carbons
						if (olefins[i].isole(h)) continue;
						int numH = olefins[i].numH(h);
						if (numH == 0) continue;					// if there are not ydrogen to abstract skip the reaction

						Carbonio alfaC = olefins[i].tipoC(h);		// save the type of substitution degree of the two carbons (the carbon we abstract the H from and the carbon invovled in the double
						Carbonio otherOleC;							// bond in beta position to it) in order to decide if the radical isomerizes or not
						otherOleC = olefins[i].tipoC(oleCarbons[1 - l]);

						int posRad, pos1Ole, pos2Ole;

						Molecola product = olefins[i];
						if (carbonioToInt(alfaC) < carbonioToInt(otherOleC)) // radical is going to isomerize
						{
							product.removeOle(oleCarbons[l], oleCarbons[1 - l]);
							product.addole(h, oleCarbons[l]);
							product.C_to_Crad(oleCarbons[1 - l]);
						}
						else												// radical is not going to isomerize
						{
							product.C_to_Crad(h);
						}
						Add(product, &allylicRads);

						//// DEBUG
						//std::cout << "Mol: " << product << ", numC = " << product.numberOfC() << ", numH = " << product.numberOfH() << ", numO = " << product.numberOfO() << std::endl;
						//// DEBUG

						std::string reacomm = k.v_allylic_rad_form(alfaC, numH);
						addReactionToList(&reactionsList, std::vector<Molecola>{ olefins[i], OH },
							std::vector<Molecola>{ product, H2O }, new double[3] { k.A, k.n, k.E },
							"Allylic radical formation:", reacomm);
						k.wrirea(reaout, 15, k.A, k.n, k.E, olefins[i], product);
					}
				}
			}
		}
		reaout << endl;
		std::cout << "Allylic radical formation added." << std::endl;


		// ALLYLIC RADICAL TO ALKENYL RO
		reaout << "Alkenyl RO formation:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Alkenyl RO formation:");

		std::vector<Molecola> alkenylRO;
		for (int i = 0; i < allylicRads.size(); i++)								// iterate through all the molecules in the allylic radicals list
		{

			Molecola product = allylicRads[i];
			product.Crad_to_COrad(product.trova(2));
			std::string reacomm = k.v_alkenyl_ro_form();
			addReactionToList(&reactionsList, std::vector<Molecola>{ allylicRads[i], HO2 },
				std::vector<Molecola>{ product, OH }, new double[3] { k.A, k.n, k.E }, "Alkenyl RO formation:", reacomm);

			Add(product, &alkenylRO);
			//// DEBUG
			//std::cout << "Mol: " << product << ", numC = " << product.numberOfC() << ", numH = " << product.numberOfH() << ", numO = " << product.numberOfO() << std::endl;
			//// DEBUG

			k.wrirea(reaout, 22, k.A, k.n, k.E, allylicRads[i], product);
		}
		reaout << endl;
		std::cout << "Alkenyl RO formation added." << std::endl;


		// ALKENYL RO DECOMPOSITION
		reaout << "Alkenyl RO decomposition:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Alkenyl RO decomposition:");

		for (int i = 0; i < alkenylRO.size(); i++)								// iterate through all the molecules in the alkenyl RO list
		{
			int posRO = alkenylRO[i].trova(9);
			int posalfa[SIZEMAX + 1];
			alkenylRO[i].scorri(posRO, 1, posalfa);		// find the carbons in alpha
			for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)		// iterate through all the carbon in alfa at the CO*
			{
				if (alkenylRO[i].isole(h)) continue;		// if it is the carbon involved in the double bond continue
				Carbonio tipoC = alkenylRO[i].tipoC(h);
				Molecola product = alkenylRO[i];
				Molecola m1, m2;
				int isPossible = alkenylRO[i].spezza(posRO, h, &m1, &m2);
				if (isPossible == 0)
				{
					std::cerr << "WARNING: alkenyl RO decomposition failed in breaking the bond, reaction skipped!" << std::endl;
					continue;
				}

				int posRadM1 = m1.trova(2);
				//m1.Crad_to_C(posRadM1);
				//m1.removeAllKeto();
				m1.addcheto(posRadM1);
				std::string reacomm = k.v_alkenyl_ro_dec(tipoC);
				addReactionToList(&reactionsList, alkenylRO[i],
					std::vector<Molecola>{ m1, m2 }, new double[3] { k.A, k.n, k.E }, "Alkenyl RO decomposition:", reacomm);
				k.wrirea(reaout, 2, k.A, k.n, k.E, alkenylRO[i], m1, m2);
			}
		}
		reaout << endl;
		std::cout << "Alkenyl RO decomposition added." << std::endl;

		// look for produced aldehydes and olefin aldehydes
		std::vector<Molecola*> aldehydes;
		std::vector<Molecola*> oleAldehydes;
		for (int i = 0; i < reactionsList.size(); i++)
		{
			std::vector<Molecola*> prod = reactionsList[i].productList();
			for (int j = 0; j < prod.size(); j++)
			{
				if (prod[j]->isAldehyde())
				{
					int a = addUniqueSpecies(&aldehydes, prod[j]);
				}
				if (prod[j]->isOleAldehyde())
				{
					int a = addUniqueSpecies(&oleAldehydes, prod[j]);
				}
			}
		}

		// DECOMPOSITION ALDEHYDES
		reaout << "Aldehydes decomposition:" << endl;
		chemout.endlineDetailed();
		chemout.writeDetReacName("Aldehydes decomposition:");

		for (int i = 0; i < aldehydes.size(); i++)								// iterate through all the aldehydes
		{
			Molecola ald = *aldehydes[i];
			Molecola m1, m2;
			int posalfa[SIZEMAX + 1];
			ald.scorri(ald.posKeto(), 1, posalfa);
			for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)		// find the carbon in alpha
			{
				ald.spezza(ald.posKeto(), h, &m1, &m2);
			}

			std::string reacomm = k.v_ald_dec();
			addReactionToList(&reactionsList, std::vector<Molecola>{ ald, OH },
				std::vector<Molecola>{ m2, CO, H2O }, new double[3] { k.A, k.n, k.E }, "Aldehydes decomposition:", reacomm);
			k.wrirea(reaout, 26, k.A, k.n, k.E, ald, m2);
		}
		reaout << endl;
		std::cout << "Aldehydes decomposition added." << std::endl;




		//

		//********************************************************************************************
		//********************************************************************************************
		// PROCESSING OF THE REACTION LIST
		//********************************************************************************************
		//********************************************************************************************
		std::cout << std::endl;
		std::cout << "Processing reaction list..." << std::endl;
		// the folowing part of code condenses every equivalent reaction (different paths with the same velocity) 
		// into one with the pre exponential factor that keep in consideration all the reactions
		for (int i = 0; i < reactionsList.size(); i++)
		{
			std::vector<int> duplicatesIndeces;
			for (int j = 0; j < reactionsList.size(); j++)	// check all the reactions to find all the equal reactions
			{
				if (j == i) continue;		// do not check the reactions against itself
				if (reactionsList[i] == reactionsList[j])
					duplicatesIndeces.push_back(j);
			}
			if (duplicatesIndeces.size() > 0)
			{
				std::cout << "-- Those reactions have been condensed toghether: " << std::endl;	
				std::cout << "      " << reactionsList[i] << std::endl;
				for (int j = 0; j < duplicatesIndeces.size(); j++)
				{
					std::cout << "      " << reactionsList[duplicatesIndeces[j]] << std::endl;
					reactionsList[i].setA(reactionsList[i].A() + reactionsList[duplicatesIndeces[j]].A());
				}
				std::cout << "   Into reaction: " << std::endl;
				std::cout << "      " << reactionsList[i] << std::endl << std::endl;

				std::sort(duplicatesIndeces.begin(), duplicatesIndeces.end());					// eliminate reaction by reaction starting from the one with the higher index so that 
				for (int j = duplicatesIndeces.size(); j > 0; j--)								// we do not modify the indexes of the reactino we still need to remove
					reactionsList.erase(reactionsList.begin() + duplicatesIndeces[j - 1]);
			}
			if (duplicatesIndeces.size() > 0)
			{
				if (reactionsList[i].comment.size() > 0)
					reactionsList[i].comment.append("; ");
				reactionsList[i].comment.append("A multiplied by ");
				reactionsList[i].comment.append(std::to_string(duplicatesIndeces.size()+1));
				reactionsList[i].comment.append(" multiple pathways");
			}
		}
		std::cout << "                             ... COMPLETED." << std::endl << std::endl;

		std::cout << "Merging equal species..." << std::endl;
		// the following part of code makes every equivalent molecule the same in order to make simpler to print equation later
		std::vector<Molecola*> speciesList;
		for (int i = 0; i < reactionsList.size(); i++)
		{
			// check for duplicate reactants
			std::vector<Molecola*> reac = reactionsList[i].reactantList();
			for (int j = 0; j < reac.size(); j++)
			{
				int a = addUniqueSpecies(&speciesList, reac[j]);
				if (a != -1) reactionsList[i].updateReactant(j, speciesList[a]);
			}
			// check for duplicate products
			std::vector<Molecola*> prod = reactionsList[i].productList();
			for (int j = 0; j < prod.size(); j++)
			{
				int a = addUniqueSpecies(&speciesList, prod[j]);
				if (a != -1) reactionsList[i].updateProduct(j, speciesList[a]);
			}
		}
		std::cout << "                             ... COMPLETED." << std::endl << std::endl;

		std::cout << "Fixing species structures..." << std::endl;
		// the following part of code fixes all the molecules to avoid linear molecules with fake branching
		for (int i = 0; i < reactionsList.size(); i++)
		{
			// fix reactants
			std::vector<Molecola*> reac = reactionsList[i].reactantList();
			for (int j = 0; j < reac.size(); j++)
			{
				reac[j]->fix();
			}
			// fix products
			std::vector<Molecola*> prod = reactionsList[i].productList();
			for (int j = 0; j < prod.size(); j++)
			{
				prod[j]->fix();
			}
		}
		std::cout << "                             ... COMPLETED." << std::endl << std::endl;

		std::cout << "Adding species to kinetic file..." << std::endl;

		for (int i = 0; i < speciesList.size(); i++)
		{
			//		std::cout << *(speciesList[i]) << std::endl;
			chemout.addSpeciesToKinFile(*(speciesList[i]));
		}

		std::cout << "                             ... COMPLETED." << std::endl << std::endl;

		//**********************************************************************************
		//******************* PRINTING SPECIES MAP FILE ************************************
		std::cout << "Printing species map file ..." << std::endl;
		//ofstream speciesMap("OUTPUT/species_map_" + chemout.molToName(HC) +".csv", std::ofstream::out);
		ofstream speciesMap("OUTPUT/species_map_.csv", std::ofstream::out);
		speciesMap << "Sep=;" << std::endl;
		std::vector<std::string> listOfNames;
		for (int i = 0; i < speciesList.size(); i++)
		{
			// debug
			std::cout << std::setw(30) << *(speciesList[i]) << std::setw(30) << chemout.molToName(*(speciesList[i])) << std::endl;
			// debug

			int size = listOfNames.size();
			//std::cout << i << "    " << * (speciesList[i]) << std::endl;
			if (size != Add(chemout.molToName(*(speciesList[i])), &listOfNames))
			{
				std::cerr << "ERROR: two species have the same name: " << chemout.molToName(*(speciesList[i])) << std::endl;
				return 0;
			}
			speciesMap << chemout.molToName(*(speciesList[i])) << ";" << *(speciesList[i]) << ";";

			if (speciesList[i]->isSpecialMolecule())
				speciesMap << std::endl;
			else
				speciesMap << speciesList[i]->inchiName() << std::endl;
		}

		speciesMap.close();
		std::cout << "                               ... COMPLETED." << std::endl << std::endl;

		//**********************************************************************************
		//******************* PRINTING SPECIES CLASS FILE **********************************
		std::cout << "Printing species class file ..." << std::endl;
		ofstream speciesClass("OUTPUT/species_class.txt", std::ofstream::out);
		for (int i = 0; i < speciesList.size(); i++)
		{
			std::string speClass;
			speClass = kindToString((speciesList[i])->kindOfSPecies());
			speciesClass << std::setw(20) << std::left << chemout.molToName(*(speciesList[i])) << " = " << speClass << std::endl;
		}

		speciesClass.close();
		std::cout << "                               ... COMPLETED." << std::endl << std::endl;


		// print all the reactions on terminal
		std::cout << std::endl;
		std::cout << "REACTION LIST:" << std::endl;
		for (int i = 0; i < reactionsList.size(); i++)
		{
			reactionsList[i].print(cout);
		}

		chemout.writeDetailedReactions(&reactionsList);
		chemout.completeDetailed();


		//**************************************************************************************
		//******************** PRINTINNG DETAILED THERMODYNAMIC FILE ***************************
		//**************************************************************************************
		std::cout << std::endl << std::endl;
		std::cout << "Printing detailed thermodynamic file ..." << std::endl;

		std::ofstream detThermoFile;
		std::string detThermoPath = "OUTPUT/thermo_detailed.inp";
		detThermoFile.open(detThermoPath);
		if (!detThermoFile)
			std::cerr << "ERROR: cannot open " + detThermoPath + "!\n      ... SKIPPED." << std::endl;
		else
		{
			ThermoOut thermout("resources/thermo/groups_contribution.csv",
				"resources/thermo/radicals_corrections.csv", "resources/baseMech/base_thermo.txt", &chemout);
			detThermoFile << "THERMO" << std::endl;
			detThermoFile << "   300.000  1000.000  5000.000" << std::endl;
			for (int i = 0; i < speciesList.size(); i++)
			{
				detThermoFile << thermout.NASAOutput(speciesList[i]);
			}
			detThermoFile << "END" << std::endl;
			detThermoFile.close();
			std::cout << "        ... COMPLETED." << std::endl;
		}

		reaout.close();
		distrib.close();

		_mkdir("ALLOUTPUTS");
		std::string finalPath = "ALLOUTPUTS\\";
		finalPath.append(molNameInp);
		_mkdir(finalPath.c_str());
		std::string command = "xcopy OUTPUT ";
		command.append(finalPath);
		command.append(" /E /Y");
		std::cout << command << std::endl;
		std::system(command.c_str());

		chemout = ChemkinOut("OUTPUT/reaction_detailed.inp");
	}
	return 0;
};

//==============================================================
//==============================================================
//==============================================================



void decomponi(Molecola HC, Kinox& k,
	Molecola* r, Molecola* roo, Molecola* qooh, Molecola* ooqooh,
	int* numr, int* numqooh)
{
	////////////////////////////////////////////////////////////////
	//  Computation of the starting radicals R*                         ITA: calcolo dei radicali di partenza R*
	////////////////////////////////////////////////////////////////
	   //r = new Molecola [HC.size()+1];

	int i;
	for (i = 1; i <= HC.size(); i++)
	{
		if (HC.tipo(i) != 1) continue;
		if (HC.tipoC(i) == Cq) continue;
		// Generate all the radicals for the 1 groups (carbon)        ITA: crea tutti i radicali da gruppi 1 (carbonio)
		// exept the quaternary carbon                                ITA: eccetto carb. quaternario
		Molecola newrad = HC;
		newrad.tipo(i, 2);
		int isomero = 0;
		for (int j = 1; j <= *numr; j++) // search if it already exist     ITA: cerca se esiste gia'
		{
			//std::cout << newrad << "  " << r[j] ;
			std::cout << newrad << "   " << r[j] << std::endl;   // debug
			if (newrad == r[j])
			{
				//std::cout << "   are equal" << std::endl;
				r[j].isomeri++;
				tipor[i] = j;    // the radical in the i position is equal to the one in j position         ITA: il rad in i e' = a quello in j
				isomero = 1;
			}
			//else {
			//	std::cout << "   are NOT equal" << std::endl;
			//}
		}

		if (!isomero)  // if it does not already exist...             ITA: se non esiste gia' ...
		{
			r[++ * numr] = newrad;
			tipor[i] = *numr;
		};
	};  // end for

	////////////////////////////////////////////////////////////////
	//  Computation of the starting radicals ROO*                       ITA: calcolo dei radicali di partenza ROO*
	////////////////////////////////////////////////////////////////
	   //roo = new Molecola [*numr+1];
	for (i = 1; i <= *numr; i++)  // i iterates trough all the r       ITA: i scorre tutti gli r
	{
		roo[i] = r[i];
		int pos_rad = 1;                                // pos_rad iterates trough all the groups           ITA: pos_rad scorre tutti i gruppi
		while (roo[i].tipo(pos_rad) != 2) pos_rad++;  // in order to find the C* (group type 2)           ITA: per trovare il C* (2)
		roo[i].tipo(pos_rad, 3);
	}; //end for i

	////////////////////////////////////////////////////////////////
	//  Computation of the starting *QOOH radicals                      ITA: calcolo dei radicali di partenza *QOOH
	////////////////////////////////////////////////////////////////
	Molecola* temp;
	//qooh = new Molecola [(*numr)*6+1];
	temp = new Molecola; if (!temp) exit(204);
	for (i = 1; i <= *numr; i++)  // i iterates trough the roo          ITA: i scorre gli roo
	{
		int isomero = 0;
		int quat = 0;
		int dist1[SIZEMAX + 1];
		int dist2[SIZEMAX + 1];
		int dist3[SIZEMAX + 1];
		int dist4[SIZEMAX + 1];
		int dist[SIZEMAX + 1];

		for (int k = 0; k <= (SIZEMAX); k++)
			dist1[k] = dist2[k] = dist3[k] = dist[k] = 0;
		//--------------------------------find the position of the i-th COO*          ITA: trova la posizione del COO* i-esimo
		int pos_oo = 1;
		while (roo[i].tipo(pos_oo) != 3) pos_oo++;
		//-----------------find all the postions of all the possible radicals...      ITA: trova tutte le posizioni dei radicali possibili...
		roo[i].scorri(pos_oo, 1, dist1);
		roo[i].scorri(pos_oo, 2, dist2);
		roo[i].scorri(pos_oo, 3, dist3);
		roo[i].scorri(pos_oo, 4, dist4);
		for (int kk = 1; kk <= SIZEMAX; kk++)
			dist[kk] = dist1[kk] + dist2[kk] + dist3[kk] + dist4[kk];
		//------------------------------------------- ...and saves them in dist[]     ITA: ... e le mette in dist[]

		for (int pos_rad = 1; pos_rad <= roo[i].size(); pos_rad++)
			if (dist[pos_rad] == 1)					//  iterates trough this COO*       ITA: scorre per questo COO*
			{										//  all the positions of the        ITA: tutte le posiz. dei
				isomero = 0;						//  found radicals                  ITA: radicali trovati
				quat = 0;
				*temp = roo[i];
				temp->tipo(pos_oo, 4);
				temp->tipo(pos_rad, 2);
				if (temp->tipoC(pos_rad) == Cq) quat = 1;    // if the carbon is quaternary              ITA: se ho fatto un quat
				// the next time it has to be skipped       ITA: dopo va saltato
				for (int j = 1; j <= *numqooh; j++)
					if (qooh[j] == *temp)
					{
						isomero = 1;
						qooh[j].isomeri++;
						tipoqooh[pos_oo][pos_rad] = j;
					};   // end if e for j
				if (!isomero && !quat) {
					qooh[++ * numqooh] = *temp;
					tipoqooh[pos_oo][pos_rad] = *numqooh;
				};
			}; // end if and for pos_rad

	};  // end for i (scorre gli roo)
	delete temp;


	////////////////////////////////////////////////////////////////
	//  Computation of all the starting radical *OOQOOH													ITA: calcolo dei radicali di partenza *OOQOOH
	////////////////////////////////////////////////////////////////
	   //ooqooh = new Molecola [*numqooh+1];
	for (i = 1; i <= *numqooh; i++)  // i iterates trough all the qooh									ITA: i scorre tutti i qooh
	{
		ooqooh[i] = qooh[i];
		int pos_rad = 1;
		while (ooqooh[i].tipo(pos_rad) != 2) pos_rad++; // iterates trough all the groups				ITA: scorre tutti i gruppi
		// in order to find the radical C (group 2)		ITA: trovare il radicale C (2)
		ooqooh[i].tipo(pos_rad, 3);
	}; //end for i

}; // end function

void help(void)
{
	std::cout << "\n";
	//std::cout << "Papox.exe v1.0 (c) Enrico Garavaglia \n";
	//std::cout << "                   email: mr..vain@bbs.infosquare.it \n";
	//std::cout << "                                                     \n";
	// std::cout << "This program is used to generate and solve    \n";
	// std::cout << "an oxidation scheme for whichever ramificated     \n";
	// std::cout << "paraffins; product ripartition for nearly null  \n";
	// std::cout << "times, the entire detailed scheme used and the    \n";
	// std::cout << "branching velocity are obtained after execution. \n";
	// std::cout << "                                                 \n";
	std::cout << "NB: This is a work in progress version of the Papox	\n";
	std::cout << "software.												\n";
	std::cout << "Currently developed by Sirio Brunialti.				\n";
	std::cout << "For any question please send an email to:				\n";
	std::cout << "														\n";
	std::cout << "                          sirio.brunialti@kaust.edu.sa\n";
	std::cout << "														\n";
	std::cout << "Original author: Enrico Garavglia						\n";
	std::cout << "          email: mr..vain@bbs.infosquare.it			\n";
	std::cout << "														\n";
	std::cout << "Syntax:												\n";
	std::cout << "   papox.exe Temp[K] Press[bar] xHC[-] mol_file		\n";
	std::cout << "														\n";
	std::cout << "   -h option = gives help								\n";
	//std::cout << "   -vn option = use metod n to compute Vbranching\n";
	std::cout << "mol_file = description of the paraffin							\n";
	std::cout << "																	\n";
	std::cout << "   EXAMPLE														\n";
	std::cout << "   n-butane is:													\n";
	std::cout << "               1 2 0 0 0											\n";
	std::cout << "               1 1 3 0 0											\n";
	std::cout << "               1 2 4 0 0											\n";
	std::cout << "               1 3 0 0 0											\n";
	std::cout << "																	\n";
	std::cout << "RULES FOR mol_file:												\n";
	std::cout << "																	\n";
	std::cout << "Each carbon atom of the molecules has to be numbered startig from	\n";
	std::cout << "1. This number will be the index used for referring to the atom	\n";
	std::cout << "																	\n";
	std::cout << "The molecules are written in form of a matrix.					\n";
	std::cout << "The matrix must have 5 columns.									\n";
	std::cout << "In each row, after the 5 numbers, anything	can be written (i.e.\n";
	std::cout << "comments).														\n";
	std::cout << "Each row represents a carbon atom. The number of the row (starting\n";
	std::cout << "the count from 1) coincides with the index of the carbon atom the \n";
	std::cout << "carbon atom the row is representing.								\n";
	std::cout << "																	\n";
	std::cout << "The first number represents the type of carbon:					\n";
	std::cout << "    - 1 normal carbon    C										\n";
	std::cout << "    - 2 radical carbon atom    C*									\n";
	std::cout << "    - 3 COO* radical												\n";
	std::cout << "    - 4 COOH														\n";
	std::cout << "    - 7 C=O														\n";
	std::cout << "    - 8 oxygen (non cyclic ether)									\n";
	std::cout << "    - 9 C-O*														\n";
	std::cout << "																	\n";
	std::cout << "The remaining four numbers are used to indicate the index of the	\n";
	std::cout << "carbon atoms bonded with this	carbon atom.						\n";
	std::cout << "If less than 4 carbons are bonded to the represented to the 		\n";
	std::cout << "represented atom, zeros must be added in order to have 5 number	\n";
	std::cout << "for row.															\n";
	std::cout << "																	\n";
	std::cout << "A double bond is defined by appending at the bottom of the matrix	\n";
	std::cout << "a row in this format:												\n";
	std::cout << "    5 x y 0 0														\n";
	std::cout << "Where x and y are the number of the two carbon atom involved in   \n";
	std::cout << "the double bond.													\n";
	std::cout << "																	\n";
	std::cout << "A cyclic oxygen is defined by appending at the bottom of the		\n";
	std::cout << "matrix a row in this format:			 							\n";
	std::cout << "    6 x y 0 0														\n";
	std::cout << "Where x and y are the number of the two carbon atoms bonded with	\n";
	std::cout << "the oxygen atom.													\n";
	std::cout << "																	\n";
	std::cout << "																	\n";
	std::cout << "																	\n";
	std::cout << "Results are reported on:											\n";
	std::cout << "	reactions.out : list of the reactions of the detailed model,	\n";
	std::cout << "			molecules are represented by their structure in Mamox	\n";
	std::cout << "			format.													\n";
	std::cout << "	reaction_detaied.inp : Chemkin compatible kinetic file	for the	\n";
	std::cout << "			detailed model.											\n";
	std::cout << "	species_map.txt : mapping of the Chemkin model names and Mamox	\n";
	std::cout << "			structures for the molecules involved in the detailed	\n";
	std::cout << "			mechanism.												\n";
	std::cout << "	reaction_lumped.inp : Chemkin compatible kinetic file for the	\n";
	std::cout << "			lumped model.											\n";
	std::cout << "	THERMO : folder containing the csv input file for the 			\n";
	std::cout << "			thermodynamic properties calculator of Cloudflame.		\n";
	std::cout << "			The chemkin compatible termodynaic data for the species	\n";
	std::cout << "			included in the Aramcomech 3.0 mechnism are collected	\n";
	std::cout << "			in the file knownThermoData.txt .						\n";
	std::cout << "																	\n";
	std::cout << "																	\n";
	std::cout << "																	\n";
	std::cout << "\n";
	exit(1);

};