#define SIZEMAX 40
#include <iostream>
#include <strstream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include "inchi_api.h"
#include "utilities.h"

//using namespace std;

#ifndef ENUM
#define ENUM
enum Carbonio { Cp, Cs, Ct, Cq };
enum Idrogeno { Hp, Hs, Ht, Hcooh };
enum Radicale { Rpmet, Rpet, Rp, Rs, Rt };
enum Anello { a5, a6, a7, a8};
enum AnelloO { ao3, ao4, ao5, ao6 };
enum Direz { dir, inv };
enum Output { MATRICE, FORMULA, FORMULA2 };
enum HAbsRad { o2, oh, h, o, ho2, ch3, c2h5 , ch3o, ch3o2, ch3oh, ch3o2h}; // JIAXIN
enum species { fuel_, R_, ROO_, QOOH_, OOQOOH_, OLE_, CO_, cEth_, RO_, KHP_, ROOH_, 
	POOH2_, ROOH2_, oleOOH_, cEthOOH_, oleR_, oleCO_, cEthR_, cEthCO_, lEthRO_, alkRO_,
	special_, unidentified_, 
	oleROO_, oleQOOH_, oleO2QOOH_, olePOOH2_, dienes_, olecEth_, olecEthR_, oleKHP_, dieneOOH_,
	ROH_, ROHOO_, QOHOOH_, oleROOH_, oleOH_, OHcEth_, O2QOHOOH_, OHoleOOH_, OHKHP_,
	POHOOH_, OHRO_, OHcEthOOH_, OROOHs_, COHCHO_, oleQO_, olecEthOOH_, OHketo_, OHoleketo_, oleROH_, dienesOH_,
	oleRO_, dienesR_, dienesCO_, diketo_, OHdiketo_, OHoleRketo_}; // JIAXIN
#endif



#ifndef MOLEC_H
#define MOLEC_H

class Molecola
{
	//  i tipi di gruppi sono:
	//        1...carbon C
	//        2...radical carbon C*
	//        3...-COO*
	//        4...-COOH
	//        5...double bond =
	//        6...-O- cyclic oxygen
	//        7...C=O oxygen
	//		  8...-O- non cyclic
	//        9...-CO*
	//        10...-COH   (Jiaxin added)
	//        11...-C*OH
	//        12...-C(OH)OO*
	//        13...-C(OH)OOH

protected:
	int struttura[SIZEMAX + 1][5 + 1];
	int ole[SIZEMAX + 1][2];
	int ete[SIZEMAX + 1][2];  // refers to cyclic ethers
	int doppioO[SIZEMAX + 1]; // =O oxygen connected with a double bond
	int numole, numete, numO;
	int numkyn; // num of alkynes.

	void init(void);                                    // usata da costr
	int  confronta(Molecola* a, Molecola* b, int gruppoa, int gruppob,
		int fattia[SIZEMAX + 1], int fattib[SIZEMAX + 1],
		int   corr[SIZEMAX + 1]);             // usata da ==
	void noempty(void);                                 // usata da spezza
	void scorri_(int, int, int[SIZEMAX + 1]);              // usata da scorri
	void Formula(std::ostream& stream, int pos, int fatti[]); // usata da >>
	int specialMolecule;	// def if the molec is special:	1 = N2
							//								2 = O2
							//								3 = O
							//								4 = OH
							//								5 = HO2
							//								6 = H2O
							//								7 = H2O2
							//								8 = H
							//								9 = H2
							//								10= CO
							//                              11= HCCO
							// JIAXIN						12=CH3O
							// JIAXIN						13=CH3OH
							// JIAXIN						14=CH3O2
							// JIAXIN						15=CH3O2H
																			

public:
	double conc;
	int isomeri;
	enum Output  output;

	//						constructors
	Molecola(void);
	Molecola(int specialMoleculeID);
	Molecola(std::string filePath);
	Molecola(int size_, std::vector<int> structure, int numOle_, std::vector<int> oleBonds, int numEthe_,
		std::vector<int> etheBonds, int numKeto_, std::vector<int> ketoPos);

	//						operators
	int operator==(Molecola);
	friend std::ostream& operator<<(std::ostream&, Molecola&);
	friend std::istream& operator>>(std::istream& stream, Molecola& obj);

	//						functions
	
	void printStruttura();	// print to cout the structure matrix
	void printDoppioO();	// print the list doppioO

	int isIsomer(Molecola* other);
	int addole(int, int);	// add double bond between the two positions
	int isole(int, int);	// true if there is a double bond between the 2 carbons in the position pos1 and pos2
	int isole(int);			// true if the carbon in position int is involved in a double bond
	int isAllylic(int position); // true if the carbon is an allylic carbon
	int isole(void);		// true if there is a double bond in the molecule
	int isCRad(void);		// true if molecule is a C* radical, false otherwise.
	int addetero(int, int); // add cyclic ether
	int isetero(int, int);	// check if there is a cyclic ether between these two positions
	int isetero(int);       // check if this position is connected to a cyclic ether
	int isetero(void);      // check if the molecule contains any cyclic ethers
	int addcheto(int);		// add a =O to this position
	int addcheto1(int);		// JIAXIN

	int ischeto(int);		// check if there is a =O connected to this position
	int ischeto(void);		// check if the molecule contains any =O
	int isKetene(int pos);		// returns:	0 if the atom is not part of a ketene
								//			1 if the atom is the -C= of the ketene
								//			2 if the atom is the =C=O of the ketene
	bool isInEthRing(int pos);	// return true if the atom in position pos is part of a cyclyc ether ring
	int isRO();					// true if the molecule is an RO (has one C* radical and a =O group)
	bool isAldehyde();			// true if the molecules is an aldehyde (no other functional group present)
	bool isOleAldehyde();		// true if the molecules is an olefin and aldehyde (no other functional group present)
	bool isKetone();			// true if the molecules is a ketone
	bool isOleKetone();			// true if the molecules is an olefin and ketone 
	int isSpecialMolecule();	//return the special molecule ID, return 0 if the molecule is not a special one
	int removeOO(int pos);		// convert a COO* group (in position pos) to a C group, if possible return 1, if not possible return 0
	int removeOOH(int pos);		// convert a COOH group (in position pos) to a C group, if possible return 1, if not possible return 0
	int removeOle(int num);		// remove the num-th olefin (order in which they are saved in the ole array), if possible return 1, if not possible return 0
	int removeOle(int pos1, int pos2);		// remove the double bond between the two carbons in pos1 and pos2, if possible return 1, if not possible return 0
	int removeCycEther(int num);// remove the num-th cyclic ether (order in which they are saved in the ete array), if possible return 1, if not possible return 0
	int removeCycEther(int pos1, int pos2);		// remove the cyclic ether between the two carbons in pos1 and pos2, if possible return 1, if not possible return 0
	void removeAllKeto();		// remove all keto from molecule
	int OO_to_OOH(int pos);		// convert a COO* group (in position pos) to a COOH group, return 1 if possible, return 0 if not possible
	int C_to_Crad(int pos);		// convert a C group (in position pos) to a C* group, return 1 if possible, return 0 if not possible
	int Crad_to_C(int pos);		// convert a C* group (in position pos) to a C group, return 1 if possible, return 0 if not possible
	int COHrad_to_C(int pos);
	int Crad_to_COOH(int pos);	// convert a C* group (in position pos) to a COOH group, return 1 if possible, return 0 if not possible
	int Crad_to_ROO(int pos);	// convert a C* group (in position pos) to a COO* group, return 1 if possible, return 0 if not possible
	int Crad_to_ROOH(int pos);	// Jiaxin here. This is for HO2 addition on allyRs.
	int ROOH_to_RO(int pos);    // convert ROOH to RO.
	int rohoo_to_orooh(int pos1, int pos2); // Wadddington.
	
	int COOrad_to_COOH(int pos);	// convert a COOrad group (in position pos) to a COOH group, return 1 if possible, return 0 if not possible
	int COOH_to_COOrad(int pos);	// convert a COOH group (in position pos) to a COOrad group, return 1 if possible, return 0 if not possible
	int COOrad_to_Crad(int pos);	// convert a COOrad group (in position pos) to a C* group, return 1 if possible, return 0 if not possible
	int COrad_to_CO(int pos);	// convert a COrad group (in position pos) to a CO group, return 1 if possible, return 0 if not possible
	int oleR_to_oleROO(int pos); // JIAXIN. convert a ole_R to ole_ROO
	int oleROO_to_oleQOOH(int pos); // JIAXIN. convert a ole_ROO to ole_QOOH
	int Crad_to_keto(int pos);  // convert a C* group (in position pos) to a C=O group, return 1 if possible, return 0 if not possible
	int Crad_to_COrad(int pos);	// convert a C* group (in position pos) to a CO* group, return 1 if possible, return 0 if not possible
	int C_to_COrad(int pos);	// convert a C group (in position pos) to a CO* group, return 1 if possible, return 0 if not possible
	// JIAXIN
	int ole_to_ROH(int pos1, int pos2);    // convert an olefin to ROH radical;

	int ole_to_R(int pos1, int pos2);   // convert an olefin to R_ radical. 
	int ole_to_QOOH(int pos1, int pos2);   // convert an olefin to QOOH_ radical.
	int ole_to_QO(int pos1, int pos2);   // Jiaxin this is for O addition to olefins.

	int addch3(int pos);  // add a ch3 on Rs.

	int removeAtom(int pos);	// remove the atom in position pos, return 1 if possible, return 0 if not possible
	bool isGroup(int pos);		// return true if the atom in pos is involved in any kind of group
	bool isNearGroup(int pos);	// return true if the atom in pos is bonded to a carbon involved in any kind of group
	int addH();					// add hydrogen where it is more stable in the molecule, return 1 if possible, return 1 if not possible
	Molecola isomerizeALLR();	// try to isomerize the alk ro, if it is possible return the isomer, toherwie return the same molecule
	int  tipo(int gruppo) { return struttura[gruppo][1]; };
	void tipo(int pos, int settipo) { struttura[pos][1] = settipo; };
	
	// legami = bonds
	int legami(int);			// return the number of bonds of the atom in the stated position
	int legato(int gruppo, int pos) { return struttura[gruppo][pos + 1]; };

	int numberOfC();				// number of carbon atoms in the molecule
	int numberOfH();				// number of hydrogens in the molecule 
	int numberOfO();				// number of oxygens in the molecule
	int numDoubleBonds(int pos);	// return the number of double bonds the atom is involved into
	int size(void);					// numero dei gruppi
	int numberOle();				// return the number of double bond present in the molecule
	int numberEthero();				// return the number of ethero present in the molecule
	int numLinEther();				// return the number of linear ether in the molecule
	int numKeto();					// return the number of cheto groups in the molecule
	int numCrad();					// return the number of C* present in the molecule
	int numKYN();
	int numCOORad();				// return the number of COO* present in the molecule
	int numCORad();					// return the number of CO* present in the molecule
	int numCOOH();					// return the number of COOH present in the molecule
	int numCOH();       // JIAXIN return the number of C-OH in the molecule.
	bool isLinear();				// return true if the molecule is linear

	int posCGroup(int group);		// return the position of the specified carbon group (C*, COO* and COOH), return 0 if the group is not present in the moelcule 
	int posCrad();					// return the position of the C* radical, return 0 if no radical is present
	int posCOOrad();				// return the position of the COO* radical, return 0 if no radical is present
	int posCOrad();					// return the position of the CO* radical, return 0 if no radical is present
	int posCOOH();					// return the position of the COOH group, return 0 if no COOH group is present
	int posOle();					// return the position of the first carbon involved in the double bond, return 0 if no double bond is present
	int posOles(); // pos in dienes
	std::vector<int> posOle(int num);	// return the position of the two carbons involved in the num-th double bond
	int posKeto();					// return the position of =O, return 0 if no =O is present
	int posKeto1();  //position of the second C=O
	std::vector<int> posEthero();	// return the two positions of the ether, return [0,0] if no ether is present
	int posCOH(); // return the position of -OH group;
	std::vector<int> posOOHinPOOH2();	// return the two positions of the OOHs in a P(OOH) molecule, return [0,0] if molecule is not a P(OOH)2
	std::vector<int> posOOHinROOH2();	// return the two positions of the OOHs in a R(OOH) molecule, return [0,0] if molecule is not a R(OOH)2

	species kindOfSPecies();

	void makeCH3();				// make the molecule into CH3
	void makeCH4();				// make the molecule into CH4
	void makeC2H5();			// make the molecule into C2H5
	void makeC2H6();			// make the molecule into C2H6
	void makeC4H72_1();			// make the molecule into C4H72_1
	void makeC4H71_3();			// make the molecule into C4H71_3
	// JIAXIN
	//void makeCH3O();			// make the molecule into CH3O
	//void makeCH3O2();			// make the molecule into CH3O2
	//void makeCH3OH();			// make the molecule into CH3O
	//void makeCH3O2H();			// make the molecule into CH3O2

	void Allylic_structure(); // convert an allylic to its equivalent structure;



	Carbonio tipoC(int gruppo);	// return the type of carbon in the selected position
	Idrogeno tipoH(int pos);
	Radicale tipoR(int pos);
	Radicale tipoROO(int pos);
	Radicale tipoROOH(int pos);
	int      numH(int);						// return the total number of H in the group (the one attached to the carbon and the ones in the group)
	int      numAbstractableH(int pos);		// return the number of H that can be abstracted (the one bonded with the C)
	int numBondedCarbons(int pos);			// return the number of cabons bonded to the atom in position pos
	int numSingleBondedCarbons(int pos);	// return the number of single bonded cabons bonded to the atom in position pos 
											// (this doesn't mean that the carbon in pos has to be connected with a double bond to the carbon, 
											// (this doesn't mean that the carbon in pos has to be connected with a double bond to the carbon, 
											// this means that the attached carbon is invoved in a double bond, also with another carbon or oxygen).
	int numDoubleBondedCarbons(int pos);	// return the number of double bonded cabons (double bonded with another carbon) bonded to the atom in position pos
	int numBondedKetoCarbons(int pos);		// return the number of keto cabons (double bonded with an oxygen) bonded to the atom in position pos

	int spezza(int pos1, int pos2, Molecola*, Molecola*); // breaks the bond between the atoms in position pos1 and position pos2
	int spezza1(int pos1, int pos2, Molecola*, Molecola*); // breaks the bond between the atoms in position pos1 and position pos2

	int trova(int gruppo);
	void scorri(int, int, int[SIZEMAX + 1]);
	int dist(int pos1, int pos2);
	int areBonded(int pos1, int pos2);		// return 1 if the two atoms in position pos1 and pos2 are bonded

	void fix();		// fix the problems of atoms not being in the right order

	std::vector<int> mainChain();		// return a vector with the positions of the carbons that makes for the main chain of the fuel
	void orderedListOfMeAndEt(std::vector<int>* chain, std::vector<int>* mePos, std::vector<int>* Et, std::vector<int>* posToNamingPos);	// overwrite the 3 input vector with the chain atoms, the position of the methyls groups and the position of the ethyl groups (the vector are ordered as they are supposed to be for the nomenclature of the molecule) the last vector is the map of the atoms as they are supposed to be numbered in the naming (including branching atoms)
	int numGraphConnections(int pos);	// return the number of bonds represented in the structure matrix for the atom in position pos

	bool hasAnyRadical();			// return true if the molecule has any kind of radical
	Molecola parentFuel();			// return the same molecules but without any functional group
	Molecola parentFuel_complete();  // without any group including olefin.
	Molecola parentFuelOH();  // without any group including olefin.
	Molecola noRadicalMolecule();	// return the same molecules but without radical groups

	bool isStructureSymmetric();	// return true if the base structure of the molecule is symmetric
	bool isSymmetric();				// return true if the molecule is symmetric

	bool isChiral();			// return true if the molecule is chiral, otherwise return false
	int numberOfRotors();		// return the number of rotors in the molecule
	int numberOfSymmetries();	// return the number of symmetries in the molecule
	bool isBondRotor(int pos1, int pos2);	// return true if the bond between pos1 and pos2 is a rotor
	std::vector<std::vector<int>> listOfCCBonds();
	
	bool isCradAllylic(int pos); // return true if the radical in position pos is allylic
	
	//JIAXIN
	int isAllyBond(int pos1, int pos2); // return 1 if the bond between Cpos1 and Cpos2 is allylic 

	std::string inchiName();
};



#endif
