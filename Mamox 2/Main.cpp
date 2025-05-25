
/*

888b     d888        d8888 888b     d888  .d88888b. Y88b   d88P  .d8888b.
8888b   d8888       d88888 8888b   d8888 d88P" "Y88b Y88b d88P  d88P  Y88b
88888b.d88888      d88P888 88888b.d88888 888     888  Y88o88P          888
888Y88888P888     d88P 888 888Y88888P888 888     888   Y888P         .d88P
888 Y888P 888    d88P  888 888 Y888P 888 888     888   d888b     .od888P"
888  Y8P  888   d88P   888 888  Y8P  888 888     888  d88888b   d88P"
888   "   888  d8888888888 888   "   888 Y88b. .d88P d88P Y88b  888".....
888       888 d88P     888 888       888  "Y88888P" d88P   Y88b 888888888

Authors:
	- Sirio Brunalti (sirio.brunialti@kaust.edu.sa)
	- Jiaxin Liu (olefin part) (J.liu8@universityofgalway.ie)
	- Enrico Garavaglia
	- Tiziano Faravelli (tiziano.faravelli@polimi.it)
	- Eliseo Ranzi


HOW TO CITE
	If a publication is produced using this software the following references
	must be reported
	- E. Ranzi, A. Frassoldati, S. Granata, T. Faravelli, Wide-Range Kinetic Modeling Study of the Pyrolysis, Partial Oxidation, and Combustion of Heavy n-Alkanes, Ind. Eng. Chem. Res. 44 (2005) 5170-5183
	- E. Ranzi, T. Faravelli, P. Gaffuri, A. Sogaro, Low-Temperature Combustion - Automatic-Generation of Primary Oxidation Reactions and Lumping Procedures, Combustion and Flame 102 (1995) 179-192.
	- S. Brunialti, X. Zhang, T. Faravelli, A. Frassoldati, S.M. Sarathy, Automatically generated detailed and lumped reaction mechanisms for low- and high-temperature oxidation of alkanes, Proceedings of the Combustion Institute 39 (2023) 335-344.

WARNING
	This version of MAMOX2 and its source code are intended for educational purposes.
	Use for commercial purposes is not permitted. For any commercial issue please contact
	Sirio Brunialti (sirio.brunialti@kaust.edu.sa) and Tiziano Faravelli (tiziano.faravelli@polimi.it).

LIMITED WARRANTY
	The Software and related documentation are provided “AS IS” and without
	any warranty of any kind and Seller EXPRESSLY DISCLAIMS ALL WARRANTIES,
	EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
	OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#define NOMINMAX		// need to be defined before inclusion of windows.h because otherwise there is a problem with the functions std::min and std::max
#include <windows.h>
#include <algorithm>
#include <thread>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
#include <direct.h>
#include <iomanip>

#include  "CKMechReader.h"
#include "utilities.h"
#include "readConfig.h"
#include "reaction.h"
#include "kinox.h"
#include "chemkinOut.h"
#include "thermoOut.h"
//#include "cantera/base/ctml.h"
#include "Simulation.h"
#include "LumpedReaction.h"

#define MAX_SIMULATION_TIME 6400

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION DECLARATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Generate the list of all possible radicals that can be formed from the input fuel
*
* @param HC Alkane the child radicals are created for
* @param Rvec pointer to the vector that is going to be filled with the child radicals
* @param cToRMap pointer to the vector that is going to contain the map correlating the
*		 index of the carbon in the structure with the corresponding radical in Rvec
* @return void
*/
void generateR(Molecola HC, std::vector<Molecola>* Rvec, std::vector<int>* cToRMap);

/**
* Generate the list of all possible ROO that can be formed from the input fuel
*
* @param HC Alkane the child ROO are created for
* @param ROOvec pointer to the vector that is going to be filled with the child ROO
* @param cToROOMap pointer to the vector that is going to contain the map correlating the
*		 index of the carbon in the structure with the corresponding ROO in ROOvec
* @return void
*/
void generateROO(Molecola HC, std::vector<Molecola>* ROOvec, std::vector<int>* cToROOMap);

/**
* Generate the list of all possible QOOH that can be formed from the input fuel
*
* @param HC Alkane the child QOOH are created for
* @param QOOHvec pointer to the vector that is going to be filled with the child QOOH
* @param cToQOOHMap pointer to the vector that is going to contain the map correlating the
*		 indexes of OOH and radical (in the order [OOHpos,Rpos] in the structure
*		 with the corresponding QOOH in QOOHvec
* @return void
*/
void generateQOOH(Molecola HC, std::vector<Molecola>* QOOHvec,
	std::vector<std::vector<int>>* cToQOOHMap);

/**
* Generate the list of all possible OOQOOH that can be formed from the input fuel
*
* @param HC Alkane the child OOQOOH are created for
* @param OOQOOHvec pointer to the vector that is going to be filled with the child OOQOOH
* @param cToOOQOOHMap pointer to the vector that is going to contain the map correlating the
*		 indexes of OOH and OO* (in the order [OOHpos,OO*pos] in the structure
*		 with the corresponding OOQOOH in OOQOOHvec
* @return void
*/
void generateOOQOOH(Molecola HC, std::vector<Molecola>* OOQOOHvec,
	std::vector<std::vector<int>>* cToOOQOOHMap);

/**
* Generate all the initiation reactions for the fuel
*
* @param HC the fuel the reaction are generate for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/

std::vector<Reaction> OHAdditionOnOlefins(Molecola HC, Kinox* k); // JIAXIN OH addition declaration.
std::vector<Reaction> O2AdditionToROH(std::vector<Molecola> ROHs, Kinox* k); // JIAXIN O2 addition ROHOO.
std::vector<Reaction> betaROHIsomerization(std::vector<Molecola> ROHs, Kinox* k); // JIAXIN beta-ROH isomerization.
std::vector<Reaction> ROHO2IsomerizationReac(std::vector<Molecola> ROHOOs, Kinox* k);  // JIAXIN ROHO2 isomerization.
std::vector<Reaction> HAdditionToOlefins(Molecola HC, Kinox* k); // JIAXIN H addition on olefin.
std::vector<Reaction> HO2AddToAllyRrecom(std::vector<Molecola> Rs, Kinox* k); // JIAXIN RA+HO2=RAOOH.
std::vector<Reaction> HO2AddToAllyRdecom(std::vector<Molecola> Rs, Kinox* k); // JIAXIN RA+HO2=RAO+OH
std::vector<Reaction> CH3AddToAllyR(std::vector<Molecola> Rs, Kinox* k);      // JIAXIN RA+CH3=RACH3
std::vector<Reaction> HO2AddOlefinToQOOH(Molecola HC, Kinox* k); // JIAXIN RA+CH3=RACH3

std::vector<Reaction> OAtomAddToOlefins(Molecola HC, Kinox* k);      // JIAXIN
std::vector<Reaction> ROOHdecomposition(std::vector<Molecola> ROOHs, Kinox* k);      // JIAXIN RAOOH=RAO+OH;
std::vector<Reaction> alkROdecomposition(std::vector<Molecola> alkROs, Kinox* k);      // JIAXIN RAO to dienes;
//std::vector<Reaction> higherOleconsumption(std::vector<Molecola> higheroles, Kinox* k);      // JIAXIN RAOOH=RAO+OH;
std::vector<Reaction> O2allyRHabstraction(std::vector<Molecola> Rs, Kinox* k);      // JIAXIN RAOOH=RAO+OH;
std::vector<Reaction> ROHO2ToWaddington(std::vector<Molecola> ROHOOs, Kinox* k);      // JIAXIN
std::vector<Reaction> WadOROOHdecomposition(std::vector<Molecola> OROOHs, Kinox* k);      // JIAXIN
std::vector<Reaction> betaROHO2Elimination(std::vector<Molecola> ROHOOs, Kinox* k);      // JIAXIN
std::vector<Reaction> betaOHQOOHScission(std::vector<Molecola> QOHOOHs, Kinox* k);      // JIAXIN
std::vector<Reaction> OHQOOHToCycEther(std::vector<Molecola> QOHOOHs, Kinox* k);      // JIAXIN
std::vector<Reaction> OHQOOHdecompositions(std::vector<Molecola> QOHOOHs, Kinox* k);      // JIAXIN
std::vector<Reaction> O2AdditionToOHQOOH(std::vector<Molecola> QOHOOHs, Kinox* k);      // JIAXIN
std::vector<Reaction> O2QOHOOHElimination(std::vector<Molecola> O2QOHOOHs, Kinox* k);      // JIAXIN
std::vector<Reaction> O2QOHOOHToKHP(std::vector<Molecola> O2QOHOOHs, Kinox* k);      // JIAXIN
std::vector<Reaction> OHKHPdecomposition(std::vector<Molecola> OHKHPs, Kinox* k);      // JIAXIN
std::vector<Reaction> O2QOHOOHIsomerization(std::vector<Molecola> O2QOHOOHs, Kinox* k);      // JIAXIN
std::vector<Reaction> POHOOHdecomposition(std::vector<Molecola> POHOOHs, Kinox* k);      // JIAXIN
std::vector<Reaction> POHOOHToEthers(std::vector<Molecola> POHOOHs, Kinox* k);      // JIAXIN

std::vector<Reaction> initiationReactions(Molecola HC, Kinox* k);
std::vector<Reaction> olefinsFromROOReactions(std::vector<Molecola> ROOs, Kinox* k); //Concerted elimination of HO2 from alkenylperoxy radicals
std::vector<Reaction> QOOHDecompositionReaction(std::vector<Molecola> QOOHs, Kinox* k); // include beta/gamma/also delta?
std::vector<Reaction> QOOHToCEthReactions(std::vector<Molecola> QOOHs, Kinox* k);
std::vector<Reaction> POOH2ToCEthOOHReactions(std::vector<Molecola> POOH2s, Kinox* k);
std::vector<Reaction> POOH2DecompositionReaction(std::vector<Molecola> POOH2s, Kinox* k); // include beta/gamma/delta.
std::vector<Reaction> O2AdditionVinyRs(std::vector<Molecola> Rs, Kinox* k);

std::vector<Reaction> oleOHConsumption(std::vector<Molecola> oleOHs1, Kinox* k, ChemkinOut* chemOut);
std::vector<Reaction> oleOHOOHConsumption(std::vector<Molecola> OHoleOOHs, Kinox* k, ChemkinOut* chemOut);
std::vector<Reaction> dienesOOHConsumption(std::vector<Molecola> dieneOOHs, Kinox* k);
std::vector<Reaction> cyclicEtherOOHDecompositionReactions(std::vector<Molecola> olecEthOOHs, Kinox* k);
std::vector<Reaction> cyclicEthersDecompositionReactions(std::vector<Molecola> olecEths, Kinox* k);
std::vector<Reaction> OHcEthOOHDecompositionReactions(std::vector<Molecola> OHcEthOOHs, Kinox* k);
std::vector<Reaction> OHcEthDecompositionReactions(std::vector<Molecola> OHcEths, Kinox* k);
std::vector<Reaction> alphaQOHOOHToKHP(std::vector<Molecola> QOHOOHs, Kinox* k);


/**
* Generate all the H-abstraction reactions for the fuel and the specific abstracting molecule
*
* @param HC the fuel the reaction are generate for
* @param absR the molecule abstracting the hydrogen
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> hAbstractionReactions(Molecola HC, Molecola absR, Kinox* k);

/**
* Generate all the O2 + R -> ROO reactions for the list of radicals
*
* @param Rs the list of radicals to generate the reactions for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> O2AdditionToRReactions(std::vector<Molecola> Rs, Kinox* k);

/**
* Generate all the ROO -> R + O2 reactions for the list of ROO
*
* @param ROOs the list of ROO to generate the reactions for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> O2EliminationFromROOReactions(std::vector<Molecola> ROOs, Kinox* k);

/**
* Generate all the O2 + QOOH -> OOQOOH reactions for the list QOOHs
*
* @param QOOHs the list of QOOHs to generate the reactions for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> O2AdditionToQOOHReactions(std::vector<Molecola> QOOHs, Kinox* k);

/**
* Generate all the OOQOOH -> QOOH + O2 reactions for the list OOQOOHs
*
* @param OOQOOHs the list of OOQOOHs to generate the reactions for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> O2EliminationFromOOQOOHReactions(std::vector<Molecola> OOQOOHs, Kinox* k);

/**
* Generate all the isomerization reactions for all the radicals in the vector Rs
*
* @param Rs the vector of radicals the reactions are going to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> RIsomerizationReaction(std::vector<Molecola> Rs, Kinox* k);

/**
* Generate all the isomerization reactions for all the ROO in the vector ROOs
*
* @param ROOs the vector of ROO the reactions are going to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> ROOIsomerizationReactions(std::vector<Molecola> ROOs, Kinox* k);

/**
* Generate all the isomerization reactions for all the QOOH in the vector QOOHs
*
* @param QOOHs the vector of QOOH the reactions are going to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> QOOHIsomerizationReactions(std::vector<Molecola> QOOHs, Kinox* k);

/**
* Generate all the isomerization reactions for all the OOQOOH in the vector OOQOOHs
*
* @param OOQOOHs the vector of OOQOOH the reactions are going to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> OOQOOHIsomerizationReactions(std::vector<Molecola> OOQOOHs, Kinox* k);

/**
* Generate all the isomerization reactions for all the P(OOH)2 in the vector POOH2s
*
* @param POOH2s the vector of P(OOH)2 the reactions are going to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> POOH2IsomerizationReactions(std::vector<Molecola> POOH2s, Kinox* k);


/**
* Generate all the reactions of beta decomposition for the provided radicals
*
* @param Rs vector of all the radicals the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> RBetaDecompositioReactions(std::vector<Molecola> Rs, Kinox* k);

/**
* Generate all the reactions of radical decomposition to olefins from H abstraction
* from O2 (R + O2 -> OLE + HO2)
*
* @param Rs vector of all the radicals the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> olefinsFromRPlusO2Reactions(std::vector<Molecola> Rs, Kinox* k);

/**
* Generate all the reactions for the decomposition of OOQOOH species to OLE-OOH
*
* @param OOQOOHs vector of all the OOQOOH the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> OOQOOHToOLEOOHReactions(std::vector<Molecola>OOQOOHs, Kinox* k);

/**
* Generate all the reactions for the decomposition of OLE-OOH species
*
* @param OLEOOHs vector of all the OLEOOH the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> OLEOOHDecompositionReactions(std::vector<Molecola> OLEOOHs, Kinox* k);

/**
* Generate all the reactions for the formation of ketohydroperoxides
*
* @param OOQOOHs vector of all the OOQOOH to decompose to ketohydroperoxides
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> KHPFormationReactions(std::vector<Molecola> OOQOOHs, Kinox* k);

/**
* Generate all the ketohydroperoxides decomposition reactions
*
* @param KHPs vector of all the ketohydroperoxides to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> KHPDecompositionReactions(std::vector<Molecola> KHPs, Kinox* k, ChemkinOut* chemOut);


/**
* Generate all the reactions of olefins' conversion to allylic radicals
*
* @param OLEs vector of all the olefins to convert
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> allylicRadicalsFormationReactions(std::vector<Molecola> OLEs, Kinox* k);

/**
* Generate all the reactions of allylic radicals' conversion to alkenyl RO
*
* @param AllRs vector of all the allylic radicals' to convert
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> alkenylROFormationReactions(std::vector<Molecola> AllRs, Kinox* k);

/**
* Generate all the reactions of alkenyl RO decomposition
*
* @param AlkROs vector of all the alkenyl RO to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> alkenylRODecompositionReactions(std::vector<Molecola> AlkROs, Kinox* k);

/**
* Generate all the aldehydes decomposition reactions
*
* @param ALD vector of all the aldehydes to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> aldehydesDecompositionReactions(std::vector<Molecola> ALDs, Kinox* k);

/**
* Generate all the aldehyde olefins decomposition reactions
*
* @param ALDOLE vector of all the aldehyde olefins to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> aldehydeOlefinsDecompositionReactions(std::vector<Molecola> ALDOLEs, Kinox* k);

/**
* Generate all the ketones decomposition reactions
*
* @param KETOs vector of all the ketones to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> ketonesDecompositionReactions(std::vector<Molecola> KETOs, Kinox* k);

/**
* Generate all the ketones olefins decomposition reactions
*
* @param KETOOLEs vector of all the ketones olefins to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> ketonesOlefinsDecompositionReactions(std::vector<Molecola> KETOOLEs, Kinox* k);


/**
* Decompose the cyclic ether radical
* @param m the cyclic ether radical to decompose
* @return a vector containig all the decomposition products
*/
std::vector<Molecola> decomposeCEthR(Molecola m);

/**
* Decompose the linear ether RO
* @param m the linear ether RO to decompose
* @return a vector containig all the decomposition products
*/
std::vector<Molecola> decomposeLinEthRO(Molecola mol);

/**
* Decompose the RO species
* @param m the RO species to decompose
* @return a vector containig all the decomposition products of the RO in vec
*/
std::vector<Molecola> decomposeRO(Molecola m2);
std::vector<Molecola> decomposeRO1(Molecola m2);
std::vector<Molecola> decomposeketo(Molecola m2); // JIAXIN
std::vector<Molecola> decomposeketo(Molecola m2); // JIAXIN
std::vector<Molecola> DecomposeoleRO(Molecola m2);
std::vector<Molecola> DecomposeDiene(Molecola m2);
/**
* Decompose the RO species and the eventually produced RO products
* @param m the RO species to decompose
* @return a vector containig all the decomposition products of the RO in vec
*/
std::vector<Molecola> fullyDecomposeRO(Molecola vec);


/**
* Decompose all the RO species in the vector and the eventually produced RO products
* @param m vector of species containing the RO to decompose
* @return a vector containig all the non RO species of vec and the decomposition products
*		  of the RO in vec
*/
std::vector<Molecola> fullyDecomposeRO(std::vector<Molecola> vec);
std::vector<Molecola> fullyDecomposeRO1(std::vector<Molecola> vec);
/**
* Get all the products that belong to a specific family of species from all the provided
* reactions. Duplicate species are skipped.
*
* @param reactions vector of reactions the products have to be taken from
* @param kind the kind of species to be taken
* @return a vector containing all the desired species
*/
std::vector<Molecola> getProducts(std::vector<Reaction> reactions, species kind);

/**
* Seek in all the products of the provided reactions for species belonging to the specified
* class of species and add to the molecule vector all the missing species.
*
* @param molVec pointer to the vector containing the list of species to be filled
* @param reactions the vector of reactions the products have to be taken from
* @param kind the kind of species to be taked
* @return integer equal to the number of new elements added to molVec
*/
int getAdditionalProducts(std::vector<Molecola>* molVec, std::vector<Reaction> reactions,
	species kind);

/**
* Convert the enumarator Carbonio to the corresponding number:	primary		-> 1
*																secondary	-> 2
*																tertiary	-> 3
*																quaternary	-> 4
* @param c the Carbonio enumerator to convert
* @return the integer
*/
int carbonioToInt(Carbonio c);

/**
* Process the vector of reactions and deals with duplicate reactions. If two reactions
* with same products and same reactants are found they are deemed to be duplicates. If
* they have also the same rate rule used (and same rate costants) they are merged into
* a single reaction and the reaction comment is modified to take into account the
* multiplication factor. If the reactions don't have the same rate rule they are flagged
* as duplicates.
*
* @param reacVec pointer to the vector to process
* @return int return the number of merged reactions
*/
int processDuplicateReactions(std::vector<Reaction>* reacVec);

/**
* Add all the new species in the reaction vector (both reactants and products) to the
* vector of molecules
*
* @param molVec pointer to the vector of molecules
* @param reacVec pointer to the vector of reactions
* @return void
*/
void addNewSpecies(std::vector<Molecola>* molVec, std::vector<Reaction>* reacVec);

/**
* Print the species in the provided file with the proper formatting
*
* @param outfile pointer to the ofstream in which the species have to be printed
* @param mols the vector of species to print
* @param label the string with the label to print as comment at the beginning
*		 of the list
* @param chemOut pointer to the ChemkinOut object providing the names
* @return void
*/
void printSpeciesInFile(std::ofstream* outfile, std::vector<Molecola> mols,
	std::string label, ChemkinOut* chemOut);

/**
* Print the species names in the provided file with the proper formatting
*
* @param outfile pointer to the ofstream in which the species have to be printed
* @param mols the vector of species' names to print
* @param label the string with the label to print as comment at the beginning
*		 of the list
* @return void
*/
void printSpeciesInFile(std::ofstream* outfile, std::vector<std::string> mols,
	std::string label);

/**
* Print the reaction in Chemkin format
*
* @param outfile pointer to the file the reaction is going to printed in
* @param raec the reaction to print
* @param chemOut pointer to the ChemkinOut object that names the molecules
* @return void
*/
void printReaction(std::ofstream* outfile, Reaction reac, ChemkinOut* chemOut);

/**
* * Print the reactions in Chemkin format
*
* @param outfile pointer to the file the reaction is going to printed in
* @param raecs the vector of reactions to print
* @param chemOut pointer to the ChemkinOut object that names the molecules
* @return void
*/
void printReactions(std::ofstream* outfile, std::vector<Reaction> reacs,
	ChemkinOut* chemOut);

/**
* Find if there is a reaction pathway consuming the species both in the base
* mechanism and in the list of generated reactions
*
* @param spec the molecule to find the consumption pathways for
* @param baseMechReacs the vector of baseReactions of the base mechanism
* @param totReactions the vector of reactions generated
* @param chemOut the pointer to the ChemkinOut object that provides the naming
* @return bool true if there is a consumption pathway otherwise false
*/
bool isThereDecompositionPath(Molecola spec,
	std::vector<baseReaction>* baseMechReacs,
	std::vector<Reaction>* totReactions, ChemkinOut* chemOut);

/**
* Return a string with the name of the kind of species
*
* @param kind the kind of species
* @return the string with the name of the kind of species
*/
std::string speciesToText(species kind);


/**
* Return true if the reaction is included in the list of base reactions.
* A reaction is deemed to be included if there is a reaction with the same
* products and reactants (it is checked also against reverse reactions).
* Kinetic parameters are not taken into account.
*
* @param reac the reaction to check if is included
* @param reacList the list of base reactions in which check if reac is present
* @param chemOut pointer to the ChemkinOut object used for naming the species
* @return bool true if the reaction is present otherwise false
*/
bool isIncluded(Reaction reac, std::vector<baseReaction>* reacList,
	ChemkinOut* chemOut);
bool isSpeciesIncluded(Molecola spec, ChemkinOut* chemOut);

/**
* Return a string with expanded environment variables (Windows)
*
* @param inString the string to expand
* @return the string with expanded environment variables
*/
std::string expandEnvVar(std::string inString);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~ MAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


int find_pos(std::vector<std::string> strList, std::string target) {
	int position = std::find(strList.begin(), strList.end(), target) - strList.begin();
	return (position < strList.size()) ? position : -1;
}



int main(int argc, char* argv[])
{
	// Set window size and position
	HWND console = GetConsoleWindow();
	RECT rec;
	GetWindowRect(console, &rec);
	MoveWindow(console, rec.left, rec.top, 1500, 800, TRUE);

	
	// ###################################################################################
	// ###################### OPEN AND READ CONFIG FILE ##################################
	// ###################################################################################
	// --------------- open file -------------------------------------
	if (argc == 1)
		UTL::fatalError("No input argument has been provided!");
	if (argc > 2)
		UTL::fatalError("Only one arguent mus be provided!");

	std::string configFilePath = argv[1];
	
	std::ifstream configFile(configFilePath);
	if (configFile.is_open() == false)
		UTL::fatalError("Unable to open config File!");
	// --------------- read file --------------------------------------
	// define variables that need to be read
	std::vector<std::string> fuelsNames;  // vector with the names of the mol files
	std::vector<double> Temps;			  // vector of temperatures [K] for the simulations 
	std::vector<double> Press;			  // vector of pressures [atm] for the simulations
	std::string baseMechanismPath;		  // path to base mechanisms' kinetic file
	std::string baseThermoPath;			  // path to base mechanisms' thermo file
	std::string baseNamingPath;			  // path to base mechanisms' glossary
	std::string outFolderPath;			  // path to output folder path
	std::string molFolderPath;			  // path to folder containing the molecules files
	std::string thermoGroupPath;		  // path to the group additivity values csv
	std::string thermoHBIPath;			  // path to the HBI values csv
	std::string rateRulesPath;			  // path to the rate rules csv
	// variables with default value
	int numCores = 2;
	bool reversibeDetailed = false;
	bool generateLumped = false;
	bool useBatch = true;
	double eqRatio = 1.0;
	double tau = 2.0;
	bool lumpWithCoreMech = false;
	bool equilibriumLumping = false;

	if (rdcfg::isKeywordPresent(&configFile, "FUELS"))
		fuelsNames = rdcfg::readStringVec(&configFile, "FUELS");
	else
		UTL::fatalError("In reading config file: No FUELS provided!");

	if (fuelsNames.size() == 0)
		UTL::fatalError("In reading config file: No FUELS provided!");

	if (rdcfg::isKeywordPresent(&configFile, "OUT_FOLDER"))
		outFolderPath = expandEnvVar(rdcfg::readString(&configFile, "OUT_FOLDER"));
	else
		UTL::fatalError("In reading config file: No output folder path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "MOL_FOLDER"))
		molFolderPath = expandEnvVar(rdcfg::readString(&configFile, "MOL_FOLDER"));
	else
		UTL::fatalError("In reading config file: No output folder path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "BASE_MECH"))
		baseMechanismPath = expandEnvVar(rdcfg::readString(&configFile, "BASE_MECH"));
	else
		UTL::fatalError("In reading config file: No base mechanism path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "BASE_THERMO"))
		baseThermoPath = expandEnvVar(rdcfg::readString(&configFile, "BASE_THERMO"));
	else
		UTL::fatalError("In reading config file: No base mechanism's thermodynamic data file path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "BASE_NAMES"))
		baseNamingPath = expandEnvVar(rdcfg::readString(&configFile, "BASE_NAMES"));
	else
		UTL::fatalError("In reading config file: No base mechanism's names file path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "THERMO_GROUP_FILE"))
		thermoGroupPath = expandEnvVar(rdcfg::readString(&configFile, "THERMO_GROUP_FILE"));
	else
		UTL::fatalError("In reading config file: No thermodynamic groups contribution file path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "THERMO_HBI_FILE"))
		thermoHBIPath = expandEnvVar(rdcfg::readString(&configFile, "THERMO_HBI_FILE"));
	else
		UTL::fatalError("In reading config file: No thermodynamic HBI file path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "RATE_RULES_FILE"))
		rateRulesPath = expandEnvVar(rdcfg::readString(&configFile, "RATE_RULES_FILE"));
	else
		UTL::fatalError("In reading config file: No rate rules file path defined!");


	if (rdcfg::isKeywordPresent(&configFile, "REVERSIBLE_DETAILED_MECH"))
	{
		std::string value = rdcfg::readString(&configFile, "REVERSIBLE_DETAILED_MECH");
		if (value == "FALSE")
			reversibeDetailed = false;
		else if (value == "TRUE")
			reversibeDetailed = true;
		else
			UTL::warning("In reading config file: Definition for keyword REVERSIBLE_DETAILED_MECH must be either TRUE or FALSE. Default value FALSE has been set.");
	}

	if (rdcfg::isKeywordPresent(&configFile, "GENERATE_LUMPED"))
	{
		std::string value = rdcfg::readString(&configFile, "GENERATE_LUMPED");
		if (value == "FALSE")
			generateLumped = false;
		else if (value == "TRUE")
			generateLumped = true;
		else
			UTL::warning("In reading config file: Definition for keyword GENERATE_LUMPED must be either TRUE or FALSE. Default value FALSE has been set.");
	}

	if (generateLumped)
	{
		if (rdcfg::isKeywordPresent(&configFile, "NUM_CORES"))
			numCores = rdcfg::readInt(&configFile, "NUM_CORES");

		if (rdcfg::isKeywordPresent(&configFile, "EQ_RATIO"))
			eqRatio = rdcfg::readDouble(&configFile, "EQ_RATIO");

		if (rdcfg::isKeywordPresent(&configFile, "TAU"))
			tau = rdcfg::readDouble(&configFile, "TAU");

		if (rdcfg::isKeywordPresent(&configFile, "LUMPING_T"))
			Temps = rdcfg::readDoubleVec(&configFile, "LUMPING_T");
		else
			UTL::fatalError("In reading config file: If GENERATE_LUMPED is TRUE, LUMPING_T must be provided!");
		if (Temps.size() == 0)
			UTL::fatalError("In reading config file: If GENERATE_LUMPED is TRUE, LUMPING_T must be provided!");


		if (rdcfg::isKeywordPresent(&configFile, "LUMPING_P"))
			Press = rdcfg::readDoubleVec(&configFile, "LUMPING_P");
		else
			UTL::fatalError("In reading config file: If GENERATE_LUMPED is TRUE, LUMPING_P must be provided!");
		if (Press.size() == 0)
			UTL::fatalError("In reading config file: If GENERATE_LUMPED is TRUE, LUMPING_P must be provided!");

		if (rdcfg::isKeywordPresent(&configFile, "SIM_TYPE"))
		{
			std::string value = rdcfg::readString(&configFile, "SIM_TYPE");
			if (value == "BATCH")
			{
				useBatch = true;
				lumpWithCoreMech = true;
			}
			else if (value == "CSTR")
			{
				useBatch = false;
				lumpWithCoreMech = false;
			}
			else if (value == "EQUILIBRIUM")
			{
				equilibriumLumping = true;
			}
			else
				UTL::warning("In reading config file: SIM_TYPE must be either BATCH, CSTR or EQUILIBRIUM. Default value BATCH has been set.");
		}

		if (rdcfg::isKeywordPresent(&configFile, "LUMP_WITH_CORE_MECH"))
		{
			std::string value = rdcfg::readString(&configFile, "LUMP_WITH_CORE_MECH");
			if (value == "FALSE")
			{
				if (useBatch)
				{
					UTL::warning("In reading config file: Keyword LUMP_WITH_CORE_MECH cannot be FALSE when SIM_TYPE is BATCH. TRUE value has been set instead.");
					lumpWithCoreMech = true;
				}
				else
				{
					lumpWithCoreMech = false;
				}
			}
			else if (value == "TRUE")
				lumpWithCoreMech = true;
			else
				UTL::warning("In reading config file: Definition for keyword LUMP_WITH_CORE_MECH must be either TRUE or FALSE. Default value has been set.");

		}
	}


	// print read values to terminal
	UTL::printTitle("VALUES READ FROM CONFIG FILE");
	std::cout << "   FUELS = " << fuelsNames[0] << std::endl;
	for (int i = 1; i < fuelsNames.size(); i++)
		std::cout << "           " << fuelsNames[i] << std::endl;
	std::cout << std::endl;
	std::cout << "   OUT_FOLDER   = " << outFolderPath << std::endl;
	std::cout << std::endl;
	std::cout << "   BASE_MECH   = " << baseMechanismPath << std::endl;
	std::cout << "   BASE_THERMO = " << baseThermoPath << std::endl;
	std::cout << "   BASE_NAMES  = " << baseNamingPath << std::endl;
	std::cout << std::endl;
	std::cout << "   THERMO_GROUP_FILE = " << thermoGroupPath << std::endl;
	std::cout << "   THERMO_HBI_FILE   = " << thermoHBIPath << std::endl;
	std::cout << std::endl;
	std::cout << "   GENRATE_LUMPED = ";
	if (generateLumped)
		std::cout << "TRUE";
	else
		std::cout << "FALSE";
	std::cout << std::endl;
	if (generateLumped)
	{
		std::cout << "   SIM_TYPE  = ";
		if (equilibriumLumping)
		{
			std::cout << "EQUILIBRIUM";
		}
		else
		{
			if (useBatch)
				std::cout << "BATCH";
			else
				std::cout << "CSTR";
		}
		std::cout << std::endl;

		std::cout << "   LUMPING_T = " << Temps[0] << std::endl;
		for (int i = 1; i < Temps.size(); i++)
			std::cout << "               " << Temps[i] << std::endl;
		std::cout << "   LUMPING_P = " << Press[0] << std::endl;
		for (int i = 1; i < Press.size(); i++)
			std::cout << "               " << Press[i] << std::endl;
		std::cout << "   NUM_CORES = " << numCores << std::endl;
	}

	// ############# CHECK IF FILES ARE AVAILABLE ########################################
	if (!UTL::checkFileExistence(baseMechanismPath))
		UTL::fatalError("Base mechanism file is not accessible.");
	if (!UTL::checkFileExistence(baseThermoPath))
		UTL::fatalError("Base mechanism thermodynamic file is not accessible.");
	if (!UTL::checkFileExistence(baseNamingPath))
		UTL::fatalError("Base mechanism naming file is not accessible.");
	if (!UTL::checkFileExistence(thermoGroupPath))
		UTL::fatalError("Thermodynamic group contribution file is not accessible.");
	if (!UTL::checkFileExistence(thermoHBIPath))
		UTL::fatalError("Thermodynamic HBI file is not accessible.");
	if (!UTL::checkFileExistence(rateRulesPath))
		UTL::fatalError("Rate rules file is not accessible.");
	if (!UTL::dirExists(molFolderPath))
		UTL::fatalError("Molecule folder " + molFolderPath + " is not accessible.");

	std::vector<std::string> fuelsPaths(fuelsNames.size());
	for (int i = 0; i < fuelsNames.size(); i++)
	{
		std::string pathToMol = molFolderPath + "\\" + fuelsNames[i] + ".hyd";
		if (!UTL::checkFileExistence(pathToMol))
			UTL::fatalError("Molecule file " + pathToMol + " is not accessible.");
		fuelsPaths[i] = pathToMol;
	}

	// ############# GENERATE REQUIRED FOLDERS ###########################################
	UTL::createDirectory(outFolderPath);

	std::string subMechsFolder = outFolderPath + "\\SUBMECHS";
	UTL::createDirectory(subMechsFolder);

	// ###################################################################################
	// ###################### OPEN AND PARSE FILES #######################################
	// ###################################################################################
	std::cout << std::endl;
	UTL::printTitle("READING FILES");

	// read base mechanism
	CKMechReader baseMech(baseMechanismPath, baseThermoPath, baseNamingPath);
	ChemkinOut chemOut(&baseMech);
	std::cout << "Base mechanism read from : " << baseMechanismPath << std::endl;
	std::cout << "                           " << baseThermoPath << std::endl;
	std::cout << "                           " << baseNamingPath << std::endl << std::endl;

	// read rate rules
	Kinox k(rateRulesPath);
	std::cout << "Rate rules read from : " << rateRulesPath << std::endl << std::endl;

	// read fuels files
	std::cout << "Reading molecules: ";
	std::vector<Molecola> fuels(fuelsNames.size());
	for (int i = 0; i < fuels.size(); i++)
	{
		fuels[i] = Molecola(fuelsPaths[i]);
		std::cout << fuels[i] << "   from " << fuelsPaths[i] << std::endl;
		std::cout << "                   ";
	}
	std::cout << std::endl << std::endl;

	// read thermo files
	ThermoOut thermOut(thermoGroupPath, thermoHBIPath, baseThermoPath, &chemOut);

	// ###################################################################################
	// ###################### INITIALIZE VARIABLES #######################################
	// ###################################################################################
	std::vector<Reaction> totalReactionsList;
	std::vector<Molecola> totalSpeciesList;
	std::vector<LumpedReaction> totalLumpedReactionsList;

	// CREATE SPECIAL MOLECULES
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
	Molecola HCCO(11);
	Molecola CH3O(12);
	Molecola CH3OH(13);
	Molecola CH3O2(14);
	Molecola CH3O2H(15);

	Molecola CH3;
	CH3.makeCH3();
	Molecola CH4;
	CH4.makeCH4();
	Molecola C2H5;
	C2H5.makeC2H5();
	Molecola C2H6;
	C2H6.makeC2H6();
	// JIAXIN
	//Molecola CH3O;
	//CH3O.makeCH3O();
	//Molecola CH3O2;
	//CH3O2.makeCH3O2();
	//Molecola CH3OH;
	//CH3OH.makeCH3OH();
	//Molecola CH3O2H;
	//CH3O2H.makeCH3O2H();

	std::vector<Molecola> specialMols = { N2, O2, O, OH, HO2, H2O, H2O2, H, H2,
		CO, CH3O, CH3OH, CH3O2, CH3O2H, CH3, CH4, C2H5, C2H6};//CH3O2, , CH3O2H};


	// ###################################################################################
	// ####################### GENERATE SUBMECHANISMS ####################################
	// ###################################################################################

	// start time measurement
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	std::ofstream timeFile;
	timeFile.open(outFolderPath+"\\timeStamps.txt");
	std::chrono::steady_clock::time_point currentTime = std::chrono::steady_clock::now();

	for (int fuelIndex = 0; fuelIndex < fuels.size(); fuelIndex++)
	{
		std::cout << std::endl << std::endl;
		UTL::printTitle(fuelsNames[fuelIndex] + " submechanism");
		std::cout << std::endl;

		// ### GENERATE FOLDERS ###
		std::cout << "Creating output folder ..." << std::endl << std::endl;
		std::string subMechFold = subMechsFolder + "\\" + fuelsNames[fuelIndex];
		UTL::createDirectory(subMechFold);


		Molecola HC = fuels[fuelIndex];

		// ### GENERATE RADICALS ####
		//std::cout << "Generating radicals ..." << std::endl;
		std::vector<Molecola> Rs;
		std::vector<int> cToRMap;
		generateR(HC, &Rs, &cToRMap);
		std::cout << "   RADICALS:" << std::endl;
		for (int i = 0; i < Rs.size(); i++)
			std::cout << "   " << i + 1 << ")   " << Rs[i] << std::endl;
		std::cout << std::endl;

		// ### GENERATE ROO ####
		//std::cout << "Generating ROO ..." << std::endl;
		std::vector<Molecola> ROOs;
		std::vector<int> cToROOMap;
		generateROO(HC, &ROOs, &cToROOMap);
		std::cout << "   ROO:" << std::endl;
		for (int i = 0; i < ROOs.size(); i++)
			std::cout << "   " << i + 1 << ")   " << ROOs[i] << std::endl;
		std::cout << std::endl;

		// ### GENERATE QOOH ####
		std::cout << "Generating QOOH ..." << std::endl;
		std::vector<Molecola> QOOHs;
		std::vector<std::vector<int>> cToQOOHMap;
		generateQOOH(HC, &QOOHs, &cToQOOHMap);
		std::cout << "   QOOH:" << std::endl;
		for (int i = 0; i < QOOHs.size(); i++)
			std::cout << "   " << i + 1 << ")   " << QOOHs[i] << std::endl;
		std::cout << std::endl;

		// ### GENERATE OOQOOH ####
		std::cout << "Generating OOQOOH ..." << std::endl;
		std::vector<Molecola> OOQOOHs;
		std::vector<std::vector<int>> cToOOQOOHMap;
		generateOOQOOH(HC, &OOQOOHs, &cToOOQOOHMap);
		std::cout << "   OOQOOH:" << std::endl;
		for (int i = 0; i < OOQOOHs.size(); i++)
			std::cout << "   " << i + 1 << ")   " << OOQOOHs[i] << std::endl;
		std::cout << std::endl;

		// ################ GENERATE REACTIONS ##########################################
		UTL::printEmbeddedString('*', "Generate reactions");

		// ### INITIATION REACTIONS (C-C breakage) ###
		std::vector<Reaction> initReac = initiationReactions(HC, &k);
		//for (auto& rea : initReac)
			//std::cout << rea << std::endl;
		std::cout << "- Class01: Unimolecular decomposition reactions added.       ("
			<< initReac.size() << " reactions)" << std::endl;

		// ### H-ABSTRACTION ###
		std::vector<Reaction> hAbsReac;
		std::vector<Reaction> hAbsReacByH = hAbstractionReactions(HC, H, &k);
		std::vector<Reaction> hAbsReacByOH = hAbstractionReactions(HC, OH, &k);
		std::vector<Reaction> hAbsReacByO = hAbstractionReactions(HC, O, &k);
		std::vector<Reaction> hAbsReacByO2 = hAbstractionReactions(HC, O2, &k);
		std::vector<Reaction> hAbsReacByHO2 = hAbstractionReactions(HC, HO2, &k);
		std::vector<Reaction> hAbsReacByCH3 = hAbstractionReactions(HC, CH3, &k);
		// JIAXIN
		std::vector<Reaction> hAbsReacByCH3O = hAbstractionReactions(HC, CH3O, &k);
		std::vector<Reaction> hAbsReacByCH3O2 = hAbstractionReactions(HC, CH3O2, &k);
		std::vector<Reaction> hAbsReacByC2H5 = hAbstractionReactions(HC, C2H5, &k);

		UTL::concatenate(&hAbsReac, &hAbsReacByH);
		UTL::concatenate(&hAbsReac, &hAbsReacByOH);
		UTL::concatenate(&hAbsReac, &hAbsReacByO);
		UTL::concatenate(&hAbsReac, &hAbsReacByO2);
		UTL::concatenate(&hAbsReac, &hAbsReacByHO2);
		UTL::concatenate(&hAbsReac, &hAbsReacByCH3);
		UTL::concatenate(&hAbsReac, &hAbsReacByCH3O);
		UTL::concatenate(&hAbsReac, &hAbsReacByCH3O2);
		UTL::concatenate(&hAbsReac, &hAbsReacByC2H5);
		//for (auto& rea : hAbsReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class02: H abstraction reactions added.       ("
			<< hAbsReac.size() << " reactions)" << std::endl;
		//
		// 
				// ### RADICALS ISOMERIZATION ###
		std::vector<Reaction> RIsomReac = RIsomerizationReaction(Rs, &k);

		//if (reversibeDetailed)
		bool rev = true;
		if (rev)
		{
			std::vector<Reaction> newVec;
			for (int i = 0; i < RIsomReac.size(); i++)
			{
				bool hasDuplicate = false;
				for (int j = 0; j < i; j++)
				{
					std::vector<Molecola*>reac1 = RIsomReac[i].reactantList();
					std::vector<Molecola*>reac2 = RIsomReac[j].reactantList();
					std::vector<Molecola*>prod1 = RIsomReac[i].productList();
					std::vector<Molecola*>prod2 = RIsomReac[j].productList();
					if (*(reac1[0]) == *(prod2[0]) && *(reac2[0]) == *(prod1[0]))
					{
						hasDuplicate = true;
						break;
					}

				}
				if (hasDuplicate == false)
					newVec.push_back(RIsomReac[i]);
			}
			RIsomReac = newVec;
		}
		//for (auto& rea : RIsomReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class03: Alkenyl radical isomerization added.       ("
			<< RIsomReac.size() << " reactions)" << std::endl;

		// ### RADICALS BETA DECOMPOSITION (R -> R' + OLE) ###
		std::vector<Reaction> RBetaDecReac = RBetaDecompositioReactions(Rs, &k);
		//for (auto& rea : RBetaDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class04: Alkenyl radical beta scissions added.       ("
			<< RBetaDecReac.size() << " reactions)" << std::endl;

		// ### Jiaxin H addition on olefins ###
		std::vector<Reaction> oleHadd = HAdditionToOlefins(HC, &k);
		std::vector<Molecola> alkylRS = getProducts(oleHadd, R_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class05: H addition to olefins reactions added.       ("
			<< oleHadd.size() << " reactions)" << std::endl;

		// ### Jiaxin O atom addition to olefins ###
		std::vector<Reaction> OaddOle = OAtomAddToOlefins(HC, &k);
		//std::vector<Molecola> alkylQOOHs = getProducts(HO2addOle, QOOH_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class06: O addition to olefins reactions added.       ("
			<< OaddOle.size() << " reactions)" << std::endl;

		// ### Jiaxin O2 addition on vinylic radicals ###
		std::vector<Reaction> O2AddVinR = O2AdditionVinyRs(Rs, &k);
		//		std::vector<Molecola> alkROs = getProducts(allyHO2add_dec, alkRO_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class07: O2 addition to vinylic radicals added.       ("
			<< O2AddVinR.size() << " reactions)" << std::endl;

		// ### Jiaxin HO2 addition to allylic radicals ###
		std::vector<Reaction> allyHO2add_rec = HO2AddToAllyRrecom(Rs, &k);
		std::vector<Molecola> ROOHs = getProducts(allyHO2add_rec, oleOOH_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class08: HO2 addition to allyRs (RA+HO2=RAOOH) added.       ("
			<< allyHO2add_rec.size() << " reactions)" << std::endl;

		// ### Jiaxin HO2 add ally producing alkRO_+OH ###
		std::vector<Reaction> allyHO2add_dec = HO2AddToAllyRdecom(Rs, &k);
		std::vector<Molecola> alkROs = getProducts(allyHO2add_dec, alkRO_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class09: HO2 addition to allyRs (RA+HO2=alkROs+OH) added.       ("
			<< allyHO2add_dec.size() << " reactions)" << std::endl;

		// ### Jiaxin RAOOH=RAO+OH ###
		std::vector<Reaction> ROOHdecom = ROOHdecomposition(ROOHs, &k);
		//		std::vector<Molecola> alkROs = getProducts(allyHO2add_dec, alkRO_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class10: RAOOH=RAO+OH reactions added.       ("
			<< ROOHdecom.size() << " reactions)" << std::endl;

		// ### Jiaxin RAO decomposition ###
		std::vector<Reaction> ROdecom = alkROdecomposition(alkROs, &k);
		//		std::vector<Molecola> alkROs = getProducts(allyHO2add_dec, alkRO_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class11: RAO decomposition reactions added.       ("
			<< ROdecom.size() << " reactions)" << std::endl;

		// ### Jiaxin CH3 addition to allylic radicals ###
		std::vector<Reaction> allyCH3add = CH3AddToAllyR(Rs, &k);
		std::vector<Molecola> Higheroles = getProducts(allyCH3add, OLE_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class12: CH3 addition to allyRs (RACH3) reactions added.       ("
			<< allyCH3add.size() << " reactions)" << std::endl;

		// ### Jiaxin HO2 addition to olefins ###
		std::vector<Reaction> HO2addOle = HO2AddOlefinToQOOH(HC, &k);
		std::vector<Molecola> alkylQOOHs = getProducts(HO2addOle, QOOH_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class13: HO2 addition to olefins reactions added.       ("
			<< HO2addOle.size() << " reactions)" << std::endl;

		// ### Jiaxin allylic radicals HAA by O2 ###
		std::vector<Reaction> O2absAllyR = O2allyRHabstraction(Rs, &k);
		std::vector<Molecola> dienes = getProducts(O2absAllyR, dienes_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class14: allyRs HAA by O2 reactions added.       ("
			<< O2absAllyR.size() << " reactions)" << std::endl;

		// ### Jiaxin OH addition on olefins ###
		std::vector<Reaction> oleOHadd = OHAdditionOnOlefins(HC, &k);
		std::vector<Molecola> ROHs = getProducts(oleOHadd, ROH_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class15: OH addition to olefins added.       ("
			<< oleOHadd.size() << " reactions)" << std::endl;

		// ### Jiaxin betaROH isomerization ###

		std::vector<Reaction> ROHisomReac = betaROHIsomerization(ROHs, &k);
		std::vector<Molecola> ROHs1 = getProducts(ROHisomReac, ROH_);
		if (reversibeDetailed)
		{
			std::vector<Reaction> newVec;
			for (int i = 0; i < ROHisomReac.size(); i++)
			{
				bool hasDuplicate = false;
				for (int j = 0; j < i; j++)
				{
					std::vector<Molecola*>reac1 = ROHisomReac[i].reactantList();
					std::vector<Molecola*>reac2 = ROHisomReac[j].reactantList();
					std::vector<Molecola*>prod1 = ROHisomReac[i].productList();
					std::vector<Molecola*>prod2 = ROHisomReac[j].productList();
					if (*(reac1[0]) == *(prod2[0]) && *(reac2[0]) == *(prod1[0]))
					{
						hasDuplicate = true;
						break;
					}

				}
				if (hasDuplicate == false)
					newVec.push_back(ROHisomReac[i]);
			}
			ROHisomReac = newVec;
		}
		//for (auto& rea : RIsomReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class16: betaROH isomerization reactions added.       ("
			<< ROHisomReac.size() << " reactions)" << std::endl;

		// ### Jiaxin O2 addition on ROH ###
		ROHs.insert(ROHs.end(), ROHs1.begin(), ROHs1.end());
		std::vector<Reaction> O2AddROH = O2AdditionToROH(ROHs, &k);
		std::vector<Molecola> ROHOOs = getProducts(O2AddROH, ROHOO_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class17: O2 addition to ROHs reactions added.       ("
			<< O2AddROH.size() << " reactions)" << std::endl;

		// ### ROHOO To Waddington (ROO -> QOOH) ###
		std::vector<Reaction> ROHO2Waddington = ROHO2ToWaddington(ROHOOs, &k);
		std::vector<Molecola> OROOHs = getProducts(ROHO2Waddington, OROOHs_);
		//for (auto& rea : ROHO2Waddington)
		//	std::cout << rea << std::endl;
		std::cout << "- Class18: betaROHO2 Waddington mechanism mechanism added.       ("
			<< ROHO2Waddington.size() << " reactions)" << std::endl;

		// ### Waddington OROOH decomposition ###
		std::vector<Reaction> Waddecom = WadOROOHdecomposition(OROOHs, &k);
		//for (auto& rea : Waddecom)
		//	std::cout << rea << std::endl;
		std::cout << "- Class19: Waddington product OROOH decompositions decompositions added.       ("
			<< Waddecom.size() << " reactions)" << std::endl;

		// ### beta-ROHO2 Elimination ###
		std::vector<Reaction> betaROHO2Elim = betaROHO2Elimination(ROHOOs, &k);
		std::vector<Molecola> oleOHs1 = getProducts(betaROHO2Elim, oleOH_);
		//for (auto& rea : betaROHO2Elim)
		//	std::cout << rea << std::endl;
		std::cout << "- Class20: ROHO2 elimination reactions added added.       ("
			<< betaROHO2Elim.size() << " reactions)" << std::endl;

		// ### ISOMERIZATION ROHOO (ROO -> QOOH) ###
		std::vector<Reaction> ROHO2Isom = ROHO2IsomerizationReac(ROHOOs, &k);
		std::vector<Molecola> QOHOOHs = getProducts(ROHO2Isom, QOHOOH_);
		//for (auto& rea : ROOIsomReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class21: ROHOO isomerization reactions added.       ("
			<< ROHO2Isom.size() << " reactions)" << std::endl;

		// ### beta-QOHOOH Elimination ###
		std::vector<Reaction> betaQOHOOHscis = betaOHQOOHScission(QOHOOHs, &k);
		std::vector<Molecola> oleOHs2 = getProducts(betaQOHOOHscis, oleOH_);
		//for (auto& rea : betaQOHOOHscis)
		//	std::cout << rea << std::endl;
		std::cout << "- Class22: beta-QOHOOH scissions added.       ("
			<< betaQOHOOHscis.size() << " reactions)" << std::endl;


		oleOHs1.insert(oleOHs1.end(), oleOHs2.begin(), oleOHs2.end());

		oleOHs1.erase(std::unique(oleOHs1.begin(), oleOHs1.end()), oleOHs1.end());

		

		// ### QOHOOH decompositions ###
		std::vector<Reaction> QOHOOHdecom = OHQOOHdecompositions(QOHOOHs, &k);
		//std::vector<Molecola> OHcEths = getProducts(QOHOOHdecom, OHcEth_);
		//for (auto& rea : QOHOOHdecom)
		//	std::cout << rea << std::endl;
		std::cout << "- Class22_combine: QOHOOH decomposition reactions added.       ("
			<< QOHOOHdecom.size() << " reactions)" << std::endl;

		// ### CycEther formation from QOHOOH ###
		std::vector<Reaction> QOHOOHToEthers = OHQOOHToCycEther(QOHOOHs, &k);
		std::vector<Molecola> OHcEths = getProducts(QOHOOHToEthers, OHcEth_);
		//for (auto& rea : QOHOOHToEthers)
		//	std::cout << rea << std::endl;
		std::cout << "- Class23: QOHOOH to cyclic ethers reactions added.       ("
			<< QOHOOHToEthers.size() << " reactions)" << std::endl;

		// ### O2 addition to QOHOOH ###
		std::vector<Reaction> O2AddQOHOOH = O2AdditionToOHQOOH(QOHOOHs, &k);
		std::vector<Molecola> O2QOHOOHs = getProducts(O2AddQOHOOH, O2QOHOOH_);
		//for (auto& rea : O2AddQOHOOH)
		//	std::cout << rea << std::endl;
		std::cout << "- Class24: O2 addition to QOHOOH reactions added.       ("
			<< O2AddQOHOOH.size() << " reactions)" << std::endl;

		// ### O2QOHOOH Elimination ###
		std::vector<Reaction> O2QOHOOHElim = O2QOHOOHElimination(O2QOHOOHs, &k);
		std::vector<Molecola> OHoleOOHs = getProducts(O2QOHOOHElim, OHoleOOH_);
		//for (auto& rea : O2QOHOOHElim)
		//	std::cout << rea << std::endl;
		std::cout << "- Class25: O2QOHOOH HO2 elimination reactions added.       ("
			<< O2QOHOOHElim.size() << " reactions)" << std::endl;

		// ### O2QOHOOH to OH-KHP ###
		std::vector<Reaction> O2QOHOOH2KHP = O2QOHOOHToKHP(O2QOHOOHs, &k);
		std::vector<Molecola> OHKHPs = getProducts(O2QOHOOH2KHP, OHKHP_);
		//for (auto& rea : O2QOHOOH2KHP)
		//	std::cout << rea << std::endl;
		std::cout << "- Class26: OH-KHP formation from O2QOHOOH reactions added.       ("
			<< O2QOHOOH2KHP.size() << " reactions)" << std::endl;

		// ### OH-KHP decomposition ###
		std::vector<Reaction> OHKHPdecom = OHKHPdecomposition(OHKHPs, &k);
		//std::vector<Molecola> OHKHPs = getProducts(O2QOHOOH2KHP, OHKHP_);
		//for (auto& rea : OHKHPdecom)
		//	std::cout << rea << std::endl;
		std::cout << "- Class27: OH-KHP decomposition reactions added.       ("
			<< OHKHPdecom.size() << " reactions)" << std::endl;

		// ### O2QOHOOH isomerization ###
		std::vector<Reaction> O2QOHOOHIsom = O2QOHOOHIsomerization(O2QOHOOHs, &k);
		std::vector<Molecola> POHOOHs = getProducts(O2QOHOOHIsom, POHOOH_);
		//for (auto& rea : O2QOHOOHIsom)
		//	std::cout << rea << std::endl;
		std::cout << "- Class28: O2QOHOOH isomerization reactions added.       ("
			<< O2QOHOOHIsom.size() << " reactions)" << std::endl;

		// ### POHOOH decomposition ###
		std::vector<Reaction> POHOOHdecom = POHOOHdecomposition(POHOOHs, &k);
		//std::vector<Molecola> OHPOOH = getProducts(POHOOHdecom, OHPOOH_);
		//for (auto& rea : POHOOHdecom)
		//	std::cout << rea << std::endl;
		std::cout << "- Class29: POHOOH decomposition reactions added.       ("
			<< POHOOHdecom.size() << " reactions)" << std::endl;

		// ### POHOOH to ether-OHOOHs ###
		std::vector<Reaction> POHOOHTocE = POHOOHToEthers(POHOOHs, &k);
		std::vector<Molecola> OHcEthOOHs = getProducts(POHOOHTocE, OHcEthOOH_);
		//for (auto& rea : POHOOHTocE)
		//	std::cout << rea << std::endl;
		std::cout << "- Class30: POHOOH to ether-OHOOHs reactions added.       ("
			<< POHOOHTocE.size() << " reactions)" << std::endl;




		//// ### Jiaxin produced higherOle consumption ###
		//std::vector<Reaction> Higheroleconsum = higherOleconsumption(alkROs, &k);
		////		std::vector<Molecola> alkROs = getProducts(allyHO2add_dec, alkRO_);
		////for (auto& rea : KHPFormReac)
		////	std::cout << rea << std::endl;
		//std::cout << "-Produced higherOle consumption added.       ("
		//	<< Higheroleconsum.size() << " reactions)" << std::endl;
		//
		
		// ### O2 ADDITION TO R (R+O2->ROO) ### Jiaxin
		std::vector<Reaction> O2addToRReac = O2AdditionToRReactions(Rs, &k);
		//for (auto& rea : O2addToRReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class31: O2 addition to oleR reactions added.       ("
			<< O2addToRReac.size() << " reactions)" << std::endl;

		// ### ROO TO OLEFINS (ROO -> ole_ole(dienes) + HO2) ###
		std::vector<Reaction> oleFromROOReac = olefinsFromROOReactions(ROOs, &k);
		std::vector<Molecola> dienes1 = getProducts(oleFromROOReac, dienes_);
		//for (auto& rea : oleFromROOReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class32: oleROO HO2 elimination reactions added.       ("
			<< oleFromROOReac.size() << " reactions)" << std::endl;
		//for (auto& prod : dienes1)
		//	std::cout << prod << std::endl;


		// ### ISOMERIZATION ROO (ROO -> QOOH) ###
		std::vector<Reaction> ROOIsomReac = ROOIsomerizationReactions(ROOs, &k);
		//for (auto& rea : ROOIsomReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class33: oleROO isomerization reactions added.       ("
			<< ROOIsomReac.size() << " reactions)" << std::endl;


		// ### DECOMPOSITION beta/gamma/delta-QOOH ###
		std::vector<Reaction> QOOHDecReac = QOOHDecompositionReaction(QOOHs, &k);
		//getAdditionalProducts(&OLEs, oleFromRPlusO2reac, OLE_);
		//for (auto& rea : betaQOOHDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class34: oleQOOH decomposition reactions added.       ("
			<< QOOHDecReac.size() << " reactions)" << std::endl;

		// ### Cyclic ethers from QOOH (QOOH -> cEth + OH) ###
		std::vector<Reaction> QOOHToCEthReac = QOOHToCEthReactions(QOOHs, &k);
		std::vector<Molecola> olecEths = getProducts(QOOHToCEthReac, olecEth_);
		//for (auto& rea : QOOHToCEthReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class35: oleQOOH to cyclic ethers reactions added.       ("
			<< QOOHToCEthReac.size() << " reactions)" << std::endl;

		std::vector<Reaction> oleOHdecom = oleOHConsumption(oleOHs1, &k, &chemOut);
		//for (auto& rea : oleOHdecom)
		//	std::cout << rea << std::endl;
		std::cout << "- Class36: OleOHs (roho2 elim) consumption reactions added.       ("
			<< oleOHdecom.size() << " reactions)" << std::endl;

		// ### O2 ADDITION TO QOOH (QOOH+O2->OOQOOH) ###
		std::vector<Reaction> O2addToQOOHReac = O2AdditionToQOOHReactions(QOOHs, &k);
		//for (auto& rea : O2addToQOOHReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class37: O2 addition to oleQOOH reactions added.       ("
			<< O2addToQOOHReac.size() << " reactions)" << std::endl;
		
		// ### KETOHYDROPEROXYDES FORMATION (OOQOOH -> KHP + OH) ###
		std::vector<Reaction> KHPFormReac = KHPFormationReactions(OOQOOHs, &k);
		std::vector<Molecola> KHPs = getProducts(KHPFormReac, oleKHP_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class38: oleKHP formation reactions added.       ("
			<< KHPFormReac.size() << " reactions)" << std::endl;

		//// ### KHP DECOMPOSITION (KHP -> DEC.PROD. + OH) ###
		std::vector<Reaction> KHPDecReac = KHPDecompositionReactions(KHPs, &k, &chemOut);
		for (auto& rea : KHPDecReac)
			std::cout << rea << std::endl;
		std::cout << "- Class39: oleKHP decomposition reactions added.       ("
			<< KHPDecReac.size() << " reactions)" << std::endl;

		std::vector<Reaction> oleOHOOHdecom = oleOHOOHConsumption(OHoleOOHs, &k, &chemOut);
		for (auto& rea : oleOHOOHdecom)
			std::cout << rea << std::endl;
		std::cout << "- Class40: OHOleOOHs (o2qohooh elim) consumption reactions added.       ("
			<< oleOHOOHdecom.size() << " reactions)" << std::endl;


		// ### ISOMERIZATION OOQOOH (OOQOOH -> P(OOH)2) ###
		std::vector<Reaction> OOQOOHIsomReac = OOQOOHIsomerizationReactions(OOQOOHs, &k);
		std::vector<Molecola> POOH2s = getProducts(OOQOOHIsomReac, olePOOH2_);
		//for (auto& prod : POOH2s)
		//	std::cout << prod << std::endl;
		//for (auto& rea : OOQOOHIsomReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class41: oleO2QOOH isomerization reactions added.       ("
			<< OOQOOHIsomReac.size() << " reactions)" << std::endl;

		// ### beta/gamma/delta-P(OOH)2 DECOMPOSITION (BETA-P(OOH)2 -> OLE-OOH + HO2) ###
		// ### GAMMA-P(OOH)2 DECOMPOSITION (GAMMA-P(OOH)2 -> OLE' + C''O + OH) ###
		// ### DELTA-P(OOH)2 DECOMPOSITION (DELTAP(OOH)2 -> OLE' + P(OOH)2 + OH) ###
		std::vector<Reaction> POOH2DecReac = POOH2DecompositionReaction(POOH2s, &k);
		std::vector<Molecola> dieneOOHs = getProducts(POOH2DecReac, dieneOOH_);
		//for (auto& rea : betaPOOH2DecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class42: oleP(OOH)2 decomposition added.       ("
			<< POOH2DecReac.size() << " reactions)" << std::endl;
		
		// ### Cyclic ethers-OOH from P(OOH)2 (P(OOH)2 -> cEthOOH + OH) ###
		std::vector<Reaction> POOH2ToCEthOOHReac = POOH2ToCEthOOHReactions(POOH2s, &k);
		std::vector<Molecola> olecEthOOHs = getProducts(POOH2ToCEthOOHReac, olecEthOOH_);
		//for (auto& rea : POOH2ToCEthOOHReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class43: P(OOH)2 to cyclic ethers-OOH reactions added.       ("
			<< POOH2ToCEthOOHReac.size() << " reactions)" << std::endl;


		// ### OOQOOH TO OLE-OOH (OOQOOH -> OLE-OOH + HO2) ###
		std::vector<Reaction> OOQOOHToOLEOOHReac = OOQOOHToOLEOOHReactions(OOQOOHs, &k);
		//getAdditionalProducts(&OLEOOHs, OOQOOHToOLEOOHReac, oleOOH_);
		//for (auto& rea : OOQOOHToOLEOOHReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class44: oleO2QOOH elimination to dienes-OOH reactions added.       ("
			<< OOQOOHToOLEOOHReac.size() << " reactions)" << std::endl;


		// ### CYCLICETHER-OOH DECOMPOSITION ###
		std::vector<Reaction> cEthOOHDecReac = cyclicEtherOOHDecompositionReactions(olecEthOOHs, &k);
		//for (auto& rea : cEthOOHDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class45: olecEthOOHs (from oleP(OOH)2) consumption reactions added.       ("
			<< cEthOOHDecReac.size() << " reactions)" << std::endl;
		
		// ### CYCLIC ETHERS DECOMPOSITION (CETH + OH -> DEC.PRO. + H2O) ###
		std::vector<Reaction> cEthDecReac = cyclicEthersDecompositionReactions(olecEths, &k);
		//for (auto& rea : cEthDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Class46: olecEths (from oleQOOH) consumption reactions added.       ("
			<< cEthDecReac.size() << " reactions)" << std::endl;

		std::vector<Reaction> OHcEthOOHdec = OHcEthOOHDecompositionReactions(OHcEthOOHs, &k);
		//for (auto& rea : OHcEthOOHdec)
		//	std::cout << rea << std::endl;
		std::cout << "- Class47: OHcEthOOHs (from POHOOH2) consumption reactions added.       ("
			<< OHcEthOOHdec.size() << " reactions)" << std::endl;

		std::vector<Reaction> OHcEthDec = OHcEthDecompositionReactions(OHcEths, &k);
		//for (auto& rea : OHcEthDec)
		//	std::cout << rea << std::endl;
		std::cout << "- Class48: OHcEthDec (from QOHOOH) consumption reactions added.       ("
			<< OHcEthDec.size() << " reactions)" << std::endl;

		std::vector<Reaction> dieneOOHdecom = dienesOOHConsumption(dieneOOHs, &k);
		//for (auto& rea : dieneOOHdecom)
		//	std::cout << rea << std::endl;
		std::cout << "- Class49: dieneOOHdecom (oleO2QOOH elim) consumption reactions added.       ("
			<< dieneOOHdecom.size() << " reactions)" << std::endl;

		std::vector<Reaction> alphaqohooh = alphaQOHOOHToKHP(QOHOOHs, &k);
		//for (auto& rea : alphaqohooh)
		//	std::cout << rea << std::endl;
		std::cout << "- Class50: alpha-QOHOOH To KHPs reactions added.       ("
			<< alphaqohooh.size() << " reactions)" << std::endl;


		//// ### OLE-OOH DECOMPOSITION (OLE-OOH -> DEC.PROD. + OH) ###
		//std::vector<Reaction> OLEOOHDecReac = OLEOOHDecompositionReactions(OLEOOHs, &k);
		////for (auto& rea : OLEOOHDecReac)
		////	std::cout << rea << std::endl;
		//std::cout << "- OLE-OOH decomposition reactions added.       ("
		//	<< OLEOOHDecReac.size() << " reactions)" << std::endl;
		//


		
		//// ### ALLYLIC RADICAL FORMATION (OLE + OH -> ALLR + H2O) ###
		//std::vector<Reaction> allRadFormRac = allylicRadicalsFormationReactions(OLEs, &k);
		//std::vector<Molecola> AllRs = getProducts(allRadFormRac, oleR_);
		////for (auto& prod : AllRs)
		////	std::cout << prod << std::endl;
		////for (auto& rea : allRadFormRac)
		////	std::cout << rea << std::endl;
		//std::cout << "- Allylic radicals formation reactions added.       ("
		//	<< allRadFormRac.size() << " reactions)" << std::endl;
		//
		//// ### ALKENYL RO FORMATION (ALLR + HO2 -> ALKRO OH) ###
		//std::vector<Reaction> alkROFormReac = alkenylROFormationReactions(AllRs, &k);
		//std::vector<Molecola> AlkROs = getProducts(alkROFormReac, alkRO_);
		////for (auto& prod : AlkROs)
		////	std::cout << prod << std::endl;
		////for (auto& rea : alkROFormReac)
		////	std::cout << rea << std::endl;
		//std::cout << "- AlkenylRO formation reactions added.       ("
		//	<< alkROFormReac.size() << " reactions)" << std::endl;
		//
		//// ### ALKENYL RO DECOMPOSITION (ALKRO -> DEC.PROD.) ###
		//std::vector<Reaction> alkRODecReac = alkenylRODecompositionReactions(AlkROs, &k);
		////for (auto& rea : alkRODecReac)
		////	std::cout << rea << std::endl;
		//std::cout << "- AlkenylRO decomposition reactions added.       ("
		//	<< alkRODecReac.size() << " reactions)" << std::endl;


		// #########################################################################
		// ############### PROCESS SPECIES AND REACTIONS ###########################
		// #########################################################################

		std::cout << std::endl;
		UTL::printEmbeddedString('*', "Processing species and reactions");
		std::vector<Molecola> allSpecies;
		UTL::concatenate(&allSpecies, &specialMols);
		addNewSpecies(&allSpecies, &initReac);
		addNewSpecies(&allSpecies, &hAbsReac);
		addNewSpecies(&allSpecies, &O2addToRReac);
		addNewSpecies(&allSpecies, &O2addToQOOHReac);
		//if (reversibeDetailed == false)
		//	addNewSpecies(&allSpecies, &O2ElimFromOOQOOHReac);
		addNewSpecies(&allSpecies, &RIsomReac);
		addNewSpecies(&allSpecies, &ROOIsomReac);
		//if (reversibeDetailed == false)
		//	addNewSpecies(&allSpecies, &QOOHIsomReac);
		addNewSpecies(&allSpecies, &OOQOOHIsomReac);
		//if (reversibeDetailed == false)
		//	addNewSpecies(&allSpecies, &POOH2IsomReac);
		addNewSpecies(&allSpecies, &oleFromROOReac);
		addNewSpecies(&allSpecies, &RBetaDecReac);
		addNewSpecies(&allSpecies, &QOOHDecReac);
		addNewSpecies(&allSpecies, &QOOHToCEthReac);
		addNewSpecies(&allSpecies, &POOH2ToCEthOOHReac);
		addNewSpecies(&allSpecies, &POOH2DecReac);
		//addNewSpecies(&allSpecies, &cEthOOHDecReac);
		addNewSpecies(&allSpecies, &OOQOOHToOLEOOHReac);
		//addNewSpecies(&allSpecies, &OLEOOHDecReac);
		addNewSpecies(&allSpecies, &KHPFormReac);
		addNewSpecies(&allSpecies, &KHPDecReac);
		//addNewSpecies(&allSpecies, &cEthDecReac);
		//addNewSpecies(&allSpecies, &allRadFormRac);
		//addNewSpecies(&allSpecies, &alkROFormReac);
		//addNewSpecies(&allSpecies, &alkRODecReac);
		addNewSpecies(&allSpecies, &oleOHadd);
		addNewSpecies(&allSpecies, &O2AddROH);
		addNewSpecies(&allSpecies, &ROHisomReac);
		addNewSpecies(&allSpecies, &ROHO2Isom);
		addNewSpecies(&allSpecies, &oleHadd);
		addNewSpecies(&allSpecies, &allyHO2add_rec);
		addNewSpecies(&allSpecies, &allyHO2add_dec);
		addNewSpecies(&allSpecies, &allyCH3add);
		addNewSpecies(&allSpecies, &HO2addOle);

		addNewSpecies(&allSpecies, &OaddOle);
		addNewSpecies(&allSpecies, &ROOHdecom);
		addNewSpecies(&allSpecies, &ROdecom);
		//addNewSpecies(&allSpecies, &Higheroleconsum);
		addNewSpecies(&allSpecies, &O2AddVinR);
		addNewSpecies(&allSpecies, &O2absAllyR);
		addNewSpecies(&allSpecies, &ROHO2Waddington);
		addNewSpecies(&allSpecies, &Waddecom);
		addNewSpecies(&allSpecies, &betaROHO2Elim);
		addNewSpecies(&allSpecies, &betaQOHOOHscis);
		addNewSpecies(&allSpecies, &QOHOOHToEthers);
		addNewSpecies(&allSpecies, &QOHOOHdecom);
		addNewSpecies(&allSpecies, &O2AddQOHOOH);
		addNewSpecies(&allSpecies, &O2QOHOOHElim);
		addNewSpecies(&allSpecies, &O2QOHOOH2KHP);
		addNewSpecies(&allSpecies, &OHKHPdecom);
		addNewSpecies(&allSpecies, &O2QOHOOHIsom);
		addNewSpecies(&allSpecies, &POHOOHdecom);
		addNewSpecies(&allSpecies, &POHOOHTocE);

		addNewSpecies(&allSpecies, &oleOHdecom);
		addNewSpecies(&allSpecies, &oleOHOOHdecom);
		addNewSpecies(&allSpecies, &dieneOOHdecom);
		addNewSpecies(&allSpecies, &cEthOOHDecReac);
		addNewSpecies(&allSpecies, &cEthDecReac);
		addNewSpecies(&allSpecies, &OHcEthOOHdec);
		addNewSpecies(&allSpecies, &OHcEthDec);
		addNewSpecies(&allSpecies, &alphaqohooh);

		for (auto& mol : allSpecies)
			mol.fix();

		std::vector<std::string> allNames(allSpecies.size());
		std::vector<std::string> SpeciesInChINameList(allSpecies.size());

		for (int i = 0; i < allSpecies.size(); i++)
		{

			if (!allSpecies[i].isSpecialMolecule())
			{
				//std::cout << allSpecies[i] << std::endl;
				SpeciesInChINameList[i] = allSpecies[i].inchiName();
				//std::cout << SpeciesInChINameList[i] << std::endl;
				allNames[i] = chemOut.molToName(allSpecies[i]);
				
			}
			else 
			{
				SpeciesInChINameList[i] = std::to_string(i);
				allNames[i] = chemOut.molToName(allSpecies[i]);
			}
		}

		//for (int i = 0; i < allSpecies.size(); i++)
		//{
		//	SpeciesInChINameList[i] = allSpecies[i].inchiName();
		//	std::cout << "allNames[j]: " << allSpecies[i] << std::endl;
		//	allNames[i] = chemOut.molToName(allSpecies[i]);
		//	//if (i == allSpecies.size() -2)
		//	//{
		//	//}
		//}
		// JIAXIN. make sure the species with the same InChIs has only one name.
		std::unordered_map<std::string, size_t> firstpos;
		for (size_t i = 0; i < SpeciesInChINameList.size(); ++i) {
			const std::string& species = SpeciesInChINameList[i];

			if (firstpos.find(species) != firstpos.end()) 
			{
				// index of first appear if possible;
				size_t firstIndex = firstpos[species];
				// allNames[firstIndex] to allNames[i]
				allNames[i] = allNames[firstIndex];
			}
			else {
				// index of first appear
				firstpos[species] = i;
			}
		}
		//for (size_t i = 0; i < allNames.size(); ++i) {
		//	std::cout << "Species: " << SpeciesInChINameList[i]
		//		<< " -> Name: " << allNames[i] << std::endl;
		//}



		for (int i = 0; i < allNames.size(); i++)
			for (int j = i + 1; j < allNames.size(); j++)
				if(allNames[i] == allNames[j])
				{
					UTL::error("Two same names found");
					std::cout << allNames[i] << "  " << allSpecies[i] << "   " << allSpecies[j] << std::endl;
				}


		std::cout << " - List of all species compiled: " << allSpecies.size()
			<< " species." << std::endl;

		int dupReac = 0;
		dupReac += processDuplicateReactions(&initReac);
		dupReac += processDuplicateReactions(&hAbsReac);
		dupReac += processDuplicateReactions(&O2addToRReac);
		//dupReac += processDuplicateReactions(&O2ElimFromROOReac);
		dupReac += processDuplicateReactions(&O2addToQOOHReac);
		//dupReac += processDuplicateReactions(&O2ElimFromOOQOOHReac);
		dupReac += processDuplicateReactions(&RIsomReac);
		dupReac += processDuplicateReactions(&ROOIsomReac);
		//dupReac += processDuplicateReactions(&QOOHIsomReac);
		dupReac += processDuplicateReactions(&OOQOOHIsomReac);
		//dupReac += processDuplicateReactions(&POOH2IsomReac);
		dupReac += processDuplicateReactions(&oleFromROOReac);
		dupReac += processDuplicateReactions(&RBetaDecReac);
		//dupReac += processDuplicateReactions(&oleFromRPlusO2reac);
		dupReac += processDuplicateReactions(&QOOHDecReac);
		dupReac += processDuplicateReactions(&QOOHToCEthReac);
		dupReac += processDuplicateReactions(&POOH2ToCEthOOHReac);
		dupReac += processDuplicateReactions(&POOH2DecReac);
		//dupReac += processDuplicateReactions(&cEthOOHDecReac);
		dupReac += processDuplicateReactions(&OOQOOHToOLEOOHReac);
		//dupReac += processDuplicateReactions(&OLEOOHDecReac);
		dupReac += processDuplicateReactions(&KHPFormReac);
		dupReac += processDuplicateReactions(&KHPDecReac);
		//dupReac += processDuplicateReactions(&cEthDecReac);
		//dupReac += processDuplicateReactions(&allRadFormRac);
		//dupReac += processDuplicateReactions(&alkROFormReac);
		//dupReac += processDuplicateReactions(&alkRODecReac);
		// JIAXIN
		dupReac += processDuplicateReactions(&oleOHadd);
		dupReac += processDuplicateReactions(&O2AddROH);
		dupReac += processDuplicateReactions(&ROHisomReac);
		dupReac += processDuplicateReactions(&ROHO2Isom);
		dupReac += processDuplicateReactions(&oleHadd);
		dupReac += processDuplicateReactions(&allyHO2add_rec);
		dupReac += processDuplicateReactions(&allyHO2add_dec);
		dupReac += processDuplicateReactions(&allyCH3add);
		dupReac += processDuplicateReactions(&HO2addOle);

		dupReac += processDuplicateReactions(&OaddOle);
		dupReac += processDuplicateReactions(&ROOHdecom);
		dupReac += processDuplicateReactions(&ROdecom);
		//dupReac += processDuplicateReactions(&Higheroleconsum);
		dupReac += processDuplicateReactions(&O2AddVinR);
		dupReac += processDuplicateReactions(&O2absAllyR);

		dupReac += processDuplicateReactions(&ROHO2Waddington);
		dupReac += processDuplicateReactions(&Waddecom);
		dupReac += processDuplicateReactions(&betaROHO2Elim);
		dupReac += processDuplicateReactions(&betaQOHOOHscis);
		dupReac += processDuplicateReactions(&QOHOOHToEthers);
		dupReac += processDuplicateReactions(&QOHOOHdecom);
		dupReac += processDuplicateReactions(&O2AddQOHOOH);
		dupReac += processDuplicateReactions(&O2QOHOOHElim);
		dupReac += processDuplicateReactions(&O2QOHOOH2KHP);
		dupReac += processDuplicateReactions(&OHKHPdecom);
		dupReac += processDuplicateReactions(&O2QOHOOHIsom);
		dupReac += processDuplicateReactions(&POHOOHdecom);
		dupReac += processDuplicateReactions(&POHOOHTocE);

		dupReac += processDuplicateReactions(&oleOHdecom);
		dupReac += processDuplicateReactions(&oleOHOOHdecom);
		dupReac += processDuplicateReactions(&dieneOOHdecom);
		dupReac += processDuplicateReactions(&cEthOOHDecReac);
		dupReac += processDuplicateReactions(&cEthDecReac);
		dupReac += processDuplicateReactions(&OHcEthOOHdec);
		dupReac += processDuplicateReactions(&OHcEthDec);
		dupReac += processDuplicateReactions(&alphaqohooh);




		std::cout << " - Duplicate reactions merged: " << dupReac <<
			" reactions merged" << std::endl;


		if (reversibeDetailed)
		{
			for (auto& reac : initReac)
				reac.setReversible(true);
			for (auto& reac : hAbsReac)
				reac.setReversible(true);
			for (auto& reac : RIsomReac)
				reac.setReversible(true);
			for (auto& reac : O2addToRReac)
				reac.setReversible(true);
			for (auto& reac : O2addToQOOHReac)
				reac.setReversible(true);
			for (auto& reac : ROOIsomReac)
				reac.setReversible(true);
			for (auto& reac : OOQOOHIsomReac)
				reac.setReversible(true);
			for (auto& reac : oleFromROOReac)
				reac.setReversible(true);
			for (auto& reac : RBetaDecReac)
				reac.setReversible(true);
			//for (auto& reac : oleFromRPlusO2reac)
			//	reac.setReversible(true);
			//for (auto& reac : betaQOOHDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : gammaQOOHDecReac)
			//	reac.setReversible(true);
			for (auto& reac : QOOHDecReac)
				reac.setReversible(true);
			for (auto& reac : QOOHToCEthReac)
				reac.setReversible(true);
			//for (auto& reac : POOH2ToCEthOOHReac)
			//	reac.setReversible(true);
			//for (auto& reac : betaPOOH2DecReac)
			//	reac.setReversible(true);
			//for (auto& reac : gammaPOOH2DecReac)
			//	reac.setReversible(true);
			//for (auto& reac : deltaPOOH2DecReac)
			//	reac.setReversible(true);
			//for (auto& reac : cEthOOHDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : OOQOOHToOLEOOHReac)
			//	reac.setReversible(true);
			//for (auto& reac : OLEOOHDecReac)
			//	reac.setReversible(true);
			for (auto& reac : KHPFormReac)
				reac.setReversible(true);
			for (auto& reac : KHPDecReac)
				reac.setReversible(true);
			//for (auto& reac : cEthDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : allRadFormRac)
			//	reac.setReversible(true);
			//for (auto& reac : alkROFormReac)
			//	reac.setReversible(true);
			//for (auto& reac : alkRODecReac)
			//	reac.setReversible(true);
			for (auto& reac : oleOHadd)
				reac.setReversible(true);
			for (auto& reac : O2AddROH)
				reac.setReversible(true);
			for (auto& reac : ROHisomReac)
				reac.setReversible(true);
			for (auto& reac : ROHO2Isom)
				reac.setReversible(true);

			for (auto& reac : oleHadd)
				reac.setReversible(true);
			for (auto& reac : allyHO2add_rec)
				reac.setReversible(true);
			for (auto& reac : allyHO2add_dec)
				reac.setReversible(true);
			for (auto& reac : allyCH3add)
				reac.setReversible(true);
			for (auto& reac : HO2addOle)
				reac.setReversible(true);

			for (auto& reac : OaddOle)
				reac.setReversible(true);
			for (auto& reac : ROOHdecom)
				reac.setReversible(true);
			for (auto& reac : ROdecom)
				reac.setReversible(true);
			//for (auto& reac : Higheroleconsum)
			//	reac.setReversible(true);
			for (auto& reac : O2AddVinR)
				reac.setReversible(true);
			for (auto& reac : O2absAllyR)
				reac.setReversible(true);

			for (auto& reac : ROHO2Waddington)
				reac.setReversible(true);
			for (auto& reac : Waddecom)
				reac.setReversible(true);
			for (auto& reac : betaROHO2Elim)
				reac.setReversible(true);
			for (auto& reac : betaQOHOOHscis)
				reac.setReversible(true);
			for (auto& reac : QOHOOHToEthers)
				reac.setReversible(true);
			for (auto& reac : QOHOOHdecom)
				reac.setReversible(true);
			for (auto& reac : O2AddQOHOOH)
				reac.setReversible(true);
			for (auto& reac : O2QOHOOHElim)
				reac.setReversible(true);
			for (auto& reac : O2QOHOOH2KHP)
				reac.setReversible(true);
			for (auto& reac : OHKHPdecom)
				reac.setReversible(true);
			for (auto& reac : O2QOHOOHIsom)
				reac.setReversible(true);
			for (auto& reac : POHOOHdecom)
				reac.setReversible(true);
			for (auto& reac : POHOOHTocE)
				reac.setReversible(true);
			for (auto& reac : OOQOOHToOLEOOHReac)
				reac.setReversible(true);
			for (auto& reac : POOH2DecReac)
				reac.setReversible(true);
			for (auto& reac : POOH2ToCEthOOHReac)
				reac.setReversible(true);

			for (auto& reac : oleOHdecom)
				reac.setReversible(true);
			for (auto& reac : oleOHOOHdecom)
				reac.setReversible(true);
			for (auto& reac : dieneOOHdecom)
				reac.setReversible(true);
			for (auto& reac : cEthOOHDecReac)
				reac.setReversible(true);
			for (auto& reac : cEthDecReac)
				reac.setReversible(true);
			for (auto& reac : OHcEthOOHdec)
				reac.setReversible(true);
			for (auto& reac : OHcEthDec)
				reac.setReversible(true);
			for (auto& reac : alphaqohooh)
				reac.setReversible(true);
			
		}

		std::vector<Reaction> allReactions;
		UTL::concatenate(&allReactions, &initReac);
		UTL::concatenate(&allReactions, &hAbsReac);
		UTL::concatenate(&allReactions, &O2addToRReac);
		//if (reversibeDetailed == false)
		//	UTL::concatenate(&allReactions, &O2ElimFromROOReac);
		UTL::concatenate(&allReactions, &O2addToQOOHReac);
		//if (reversibeDetailed == false)
		//	UTL::concatenate(&allReactions, &O2ElimFromOOQOOHReac);
		UTL::concatenate(&allReactions, &RIsomReac);
		UTL::concatenate(&allReactions, &ROOIsomReac);
		//if (reversibeDetailed == false)
		//	UTL::concatenate(&allReactions, &QOOHIsomReac);
		UTL::concatenate(&allReactions, &OOQOOHIsomReac);
		//if (reversibeDetailed == false)
		//	UTL::concatenate(&allReactions, &POOH2IsomReac);
		UTL::concatenate(&allReactions, &oleFromROOReac);
		UTL::concatenate(&allReactions, &RBetaDecReac);
		//UTL::concatenate(&allReactions, &oleFromRPlusO2reac);
		UTL::concatenate(&allReactions, &QOOHDecReac);
		UTL::concatenate(&allReactions, &QOOHToCEthReac);
		UTL::concatenate(&allReactions, &POOH2ToCEthOOHReac);
		UTL::concatenate(&allReactions, &POOH2DecReac);
		//UTL::concatenate(&allReactions, &cEthOOHDecReac);
		UTL::concatenate(&allReactions, &OOQOOHToOLEOOHReac);
		//UTL::concatenate(&allReactions, &OLEOOHDecReac);
		UTL::concatenate(&allReactions, &KHPFormReac);
		UTL::concatenate(&allReactions, &KHPDecReac);
		//UTL::concatenate(&allReactions, &cEthDecReac);
		//UTL::concatenate(&allReactions, &allRadFormRac);
		//UTL::concatenate(&allReactions, &alkROFormReac);
		//UTL::concatenate(&allReactions, &alkRODecReac);
		// JIAXIN
		UTL::concatenate(&allReactions, &oleOHadd);
		UTL::concatenate(&allReactions, &O2AddROH);
		UTL::concatenate(&allReactions, &ROHisomReac);
		UTL::concatenate(&allReactions, &ROHO2Isom);

		UTL::concatenate(&allReactions, &oleHadd);
		UTL::concatenate(&allReactions, &allyHO2add_rec);
		UTL::concatenate(&allReactions, &allyHO2add_dec);
		UTL::concatenate(&allReactions, &allyCH3add);
		UTL::concatenate(&allReactions, &HO2addOle);

		UTL::concatenate(&allReactions, &OaddOle);
		UTL::concatenate(&allReactions, &ROOHdecom);
		UTL::concatenate(&allReactions, &ROdecom);
		//UTL::concatenate(&allReactions, &Higheroleconsum);
		//UTL::concatenate(&allReactions, &O2AddVinR);
		UTL::concatenate(&allReactions, &O2absAllyR);

		UTL::concatenate(&allReactions, &ROHO2Waddington);
		UTL::concatenate(&allReactions, &Waddecom);
		UTL::concatenate(&allReactions, &betaROHO2Elim);
		UTL::concatenate(&allReactions, &betaQOHOOHscis);
		UTL::concatenate(&allReactions, &QOHOOHToEthers);
		UTL::concatenate(&allReactions, &QOHOOHdecom);
		UTL::concatenate(&allReactions, &O2AddQOHOOH);
		UTL::concatenate(&allReactions, &O2QOHOOHElim);
		UTL::concatenate(&allReactions, &O2QOHOOH2KHP);
		UTL::concatenate(&allReactions, &OHKHPdecom);
		UTL::concatenate(&allReactions, &O2QOHOOHIsom);
		UTL::concatenate(&allReactions, &POHOOHdecom);
		UTL::concatenate(&allReactions, &POHOOHTocE);
		UTL::concatenate(&allReactions, &oleOHdecom);
		UTL::concatenate(&allReactions, &O2AddVinR);
		UTL::concatenate(&allReactions, &oleOHOOHdecom);
		UTL::concatenate(&allReactions, &dieneOOHdecom);
		UTL::concatenate(&allReactions, &cEthOOHDecReac);
		UTL::concatenate(&allReactions, &cEthDecReac);
		UTL::concatenate(&allReactions, &OHcEthOOHdec);
		UTL::concatenate(&allReactions, &OHcEthDec);
		UTL::concatenate(&allReactions, &alphaqohooh);


		std::cout << " - List of reactions created: " << allReactions.size() <<
			" reactions" << std::endl;


		// #########################################################################
		// ####################### PRINT SUBMECHANISM ##############################
		// #########################################################################
		std::cout << std::endl;
		UTL::printEmbeddedString('*', "Printing detailed mechanism");
		std::string detailedMechanismPath = subMechFold + "\\kinetic_detailed.inp";
		std::ofstream detailedMechanism(detailedMechanismPath);
		detailedMechanism << "! " << chemOut.molToName(HC) << " COMBUSTION MECHANISM" << std::endl;
		detailedMechanism << "! " << HC << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "ELEMENTS" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "C" << std::endl;
		detailedMechanism << "H" << std::endl;
		detailedMechanism << "O" << std::endl;
		detailedMechanism << "N" << std::endl;
		detailedMechanism << "AR" << std::endl;
		detailedMechanism << "HE" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "END" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "SPECIES" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "! Generic species" << std::endl;
		//detailedMechanism << "O2     H2    H2O2    H2O   N2   AR    HE   CO   CO2" 
		detailedMechanism << "O2                 H2                 H2O2               H2O"
			<< std::endl;
		detailedMechanism << "N2                 AR                 HE                 CO"
			<< std::endl;
		detailedMechanism << "CO2                HO2                OH                 H"
			<< std::endl;
		detailedMechanism << "CH3                C2H5               CH3O               CH3OH" 
			<< std::endl;
		detailedMechanism << "CH3O2              CH3O2H" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "!++++++++++++++++++ " << chemOut.molToName(HC)
			<< " SUBMECHANISM ++++++++++++++++++" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "! Fuel" << std::endl;
		detailedMechanism << chemOut.molToName(HC) << std::endl;
		detailedMechanism << "" << std::endl;
		printSpeciesInFile(&detailedMechanism, Rs, "oleR Alkenyl radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, ROOs, "oleROO alkenyl peroxy radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, QOOHs, "oleQOOH hydroperoxyalkenyl radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, OOQOOHs, "oleO2QOOH Hydroperoxyalkenylperoxy radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, POOH2s, "oleP(OOH)2 dihydroperoxy alkenyl radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, KHPs, "oleKHPs", &chemOut);

		olecEths.insert(olecEths.end(), olecEthOOHs.begin(), olecEthOOHs.end());
		olecEths.erase(std::unique(olecEths.begin(), olecEths.end()), olecEths.end());
		printSpeciesInFile(&detailedMechanism, olecEths, "olecycEths/olecycEthOOHs from oleQOOH/olePOOH elimination", &chemOut);
		printSpeciesInFile(&detailedMechanism, ROHs, "ROHs alcoholic (hydroxyl) alkyl radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, ROHOOs, "ROHOO alcoholic peroxy radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, QOHOOHs, "QOHOOHs alcoholic hydroperoxylalkyl radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, O2QOHOOHs, "O2QOHOOHs alcoholic hydroperoxyl-alkylperoxy radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, POHOOHs, "POHOOHs alcoholic dihydroperoxylalkyl radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, OHKHPs, "OHKHPs", &chemOut);
		ROOHs.insert(ROOHs.end(), alkROs.begin(), alkROs.end());
		ROOHs.erase(std::unique(ROOHs.begin(), ROOHs.end()), ROOHs.end());
		printSpeciesInFile(&detailedMechanism, ROOHs, "RAOOH allylic ROOH and subsequents", &chemOut);
		printSpeciesInFile(&detailedMechanism, Higheroles, "RACH3 Higher oles", &chemOut);
		printSpeciesInFile(&detailedMechanism, OROOHs, "Waddington mechanism", &chemOut); 

		oleOHs1.insert(oleOHs1.end(), OHoleOOHs.begin(), OHoleOOHs.end());
		//std::sort(oleOHs1.begin(), oleOHs1.end());
		oleOHs1.erase(std::unique(oleOHs1.begin(), oleOHs1.end()), oleOHs1.end());
		printSpeciesInFile(&detailedMechanism, oleOHs1, "oleOHs/oleOHOOHs from ROHOO/O2QOHOOH elimination", &chemOut);
		OHcEths.insert(OHcEths.end(), OHcEthOOHs.begin(), OHcEthOOHs.end());
		OHcEths.erase(std::unique(OHcEths.begin(), OHcEths.end()), OHcEths.end());
		printSpeciesInFile(&detailedMechanism, OHcEths, "OHcycEths/OHcycEthOOHs from QOHOOH/POHOOH elimination", &chemOut);
		dienes.insert(dienes.end(), dienes1.begin(), dienes1.end());
		dienes.erase(std::unique(dienes.begin(), dienes.end()), dienes.end());
		printSpeciesInFile(&detailedMechanism, dienes, "Dienes", &chemOut);
		printSpeciesInFile(&detailedMechanism, dieneOOHs, "DieneOOHs", &chemOut);


		//printSpeciesInFile(&detailedMechanism, cEths, "cEth cyclic ethers", &chemOut);
		//printSpeciesInFile(&detailedMechanism, cEthOOHs, "CEth-OOH", &chemOut);
		//printSpeciesInFile(&detailedMechanism, AllRs, "AllR allylic radicals", &chemOut);
		//printSpeciesInFile(&detailedMechanism, AlkROs, "AlkRO alkenyl RO", &chemOut);
		detailedMechanism << "" << std::endl;

		//std::vector<Molecola> fuel_dec(0);
		std::vector<Molecola> R_dec(0);
		std::vector<Molecola> ROO_dec(0);
		std::vector<Molecola> QOOH_dec(0);
		std::vector<Molecola> OOQOOH_dec(0);
		std::vector<Molecola> OLE_dec(0);
		std::vector<Molecola> CO_dec(0);
		std::vector<Molecola> cEth_dec(0);
		std::vector<Molecola> RO_dec(0);
		std::vector<Molecola> KHP_dec(0);
		std::vector<Molecola> POOH2_dec(0);
		std::vector<Molecola> oleOOH_dec(0);
		std::vector<Molecola> cEthOOH_dec(0);
		std::vector<Molecola> oleR_dec(0);
		std::vector<Molecola> oleCO_dec(0);
		std::vector<Molecola> cEthR_dec(0);
		std::vector<Molecola> cEthCO_dec(0);
		std::vector<Molecola> alkRO_dec(0);
		//std::vector<Molecola> oleRO2_elim(0);
		//std::vector<Molecola> oleO2QOOH_elim(0);
		std::vector<Molecola> small_dienes(0);
		std::vector<Molecola> oleKHP(0);
		//std::vector<Molecola> oleKHP_form(0);
		std::vector<Molecola> oleOH(0);
		std::vector<Molecola> oleketo_radical(0);
		std::vector<Molecola> OHoleOOH(0);
		std::vector<Molecola> OHoleKeto(0);
		std::vector<Molecola> oleROH(0);
		std::vector<Molecola> dienesOH(0);
		std::vector<Molecola> ROH(0);
		std::vector<Molecola> dienesR(0);
		std::vector<Molecola> diketo(0);
		std::vector<Molecola> dienesCO(0);
		std::vector<Molecola> OHRO(0);
		std::vector<Molecola> QOHOOH(0);
		std::vector<Molecola> OHketo(0);
		std::vector<Molecola> OHKHP(0);
		std::vector<Molecola> OHoleketo(0);

		std::vector<Molecola> special_dec(0);


		for (auto& mol : allSpecies)
		{
			if (mol == CH3 || mol == CH4 || mol == C2H5 || mol == C2H6)
				continue;
			//if ((mol.parentFuel() == HC || mol.parentFuel_complete() == HC.parentFuel_complete()) || mol.size() == HC.size()) // JIAXIN
			if (mol.size() == HC.size())
				continue;
			species typ = mol.kindOfSPecies();
			switch (typ)
			{
				//case fuel_:
				//	fuel_dec.push_back(mol);
				//	break;
			case R_:
				R_dec.push_back(mol);
				break;
			case ROO_:
				ROO_dec.push_back(mol);
				break;
			case QOOH_:
				QOOH_dec.push_back(mol);
				break;
			case OOQOOH_:
				OOQOOH_dec.push_back(mol);
				break;
			case OLE_:
				OLE_dec.push_back(mol);
				break;
			case CO_:
				CO_dec.push_back(mol);
				break;
			case cEth_:
				cEth_dec.push_back(mol);
				break;
			case RO_:
				RO_dec.push_back(mol);
				break;
			case KHP_:
				KHP_dec.push_back(mol);
				break;
			case POOH2_:
				POOH2_dec.push_back(mol);
				break;
			case oleOOH_:
				oleOOH_dec.push_back(mol);
				break;
			case cEthOOH_:
				cEthOOH_dec.push_back(mol);
				break;
			case oleR_:
				oleR_dec.push_back(mol);
				break;
			case oleCO_:
				oleCO_dec.push_back(mol);
				break;
			case cEthR_:
				cEthR_dec.push_back(mol);
				break;
			case cEthCO_:
				cEthCO_dec.push_back(mol);
				break;
			case alkRO_:
				alkRO_dec.push_back(mol);
				break;
			// JIAXIN
			//case dienes_:
			//	oleRO2_elim.push_back(mol);
			//	break;
			//case dieneOOH_:
			//	oleO2QOOH_elim.push_back(mol);
			//	break;
			case oleKHP_:
				oleKHP.push_back(mol);
				break;
			case dienes_:
				small_dienes.push_back(mol);
				break;
			case oleOH_:
				oleOH.push_back(mol);
				break;
			case oleRO_:
				oleketo_radical.push_back(mol);
				break;
			case OHoleOOH_:
				OHoleOOH.push_back(mol);
				break;
			case OHoleketo_:
				OHoleKeto.push_back(mol);
				break;
			case oleROH_:
				oleROH.push_back(mol);
				break;
			case dienesOH_:
				dienesOH.push_back(mol);
				break;
			case ROH_:
				ROH.push_back(mol);
				break;
			case dienesR_:
				dienesR.push_back(mol);
				break;
			case diketo_:
				diketo.push_back(mol);
				break;
			case dienesCO_:
				dienesCO.push_back(mol);
				break;
			case OHRO_:
				OHRO.push_back(mol);
				break;
			case QOHOOH_:
				QOHOOH.push_back(mol);
				break;
			case COHCHO_:
				OHketo.push_back(mol);
				break;
			case OHKHP_:
				OHKHP.push_back(mol);
				break;
				
			case special_:
				break;
				//oleROO_, oleQOOH_, oleO2QOOH_, olePOOH2_, dienes_, olecEth_, , dieneOOH_,
				//	ROH_, , QOHOOH_, oleROOH_, oleOH_, OHcEth_, O2QOHOOH_, OHoleOOH_, OHKHP_,
				//	POHOOH_, OHRO_, OHcEthOOH_, OROOHs_, COHCHO_
			default:
				UTL::warning("The following decomposition species should not be present");
				std::cout << mol << std::endl;
				break;
			}
		}
		detailedMechanism << "!++++++++++++++++++ DECOMPOSITION PRODUCTS ++++++++++++++++++"
			<< std::endl;
		printSpeciesInFile(&detailedMechanism, R_dec, "R", &chemOut);
		printSpeciesInFile(&detailedMechanism, ROO_dec, "ROO", &chemOut);
		printSpeciesInFile(&detailedMechanism, QOOH_dec, "QOOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, OOQOOH_dec, "OOQOOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, OLE_dec, "OLE", &chemOut);
		printSpeciesInFile(&detailedMechanism, CO_dec, "CO aldehydes/ketones", &chemOut);
		printSpeciesInFile(&detailedMechanism, cEth_dec, "cEth", &chemOut);
		printSpeciesInFile(&detailedMechanism, RO_dec, "RO", &chemOut);
		printSpeciesInFile(&detailedMechanism, KHP_dec, "KHP", &chemOut);
		printSpeciesInFile(&detailedMechanism, POOH2_dec, "POOH2", &chemOut);
		printSpeciesInFile(&detailedMechanism, oleOOH_dec, "OLE-OOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, cEthOOH_dec, "cEth-OOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, oleR_dec, "OLE-R", &chemOut);
		printSpeciesInFile(&detailedMechanism, oleCO_dec, "OLE-CO", &chemOut);
		printSpeciesInFile(&detailedMechanism, cEthR_dec, "cEth-R", &chemOut);
		printSpeciesInFile(&detailedMechanism, cEthCO_dec, "cEth-CO", &chemOut);
		printSpeciesInFile(&detailedMechanism, alkRO_dec, "AlkRO", &chemOut);
		printSpeciesInFile(&detailedMechanism, small_dienes, "smaller dienes", &chemOut);
		printSpeciesInFile(&detailedMechanism, oleKHP, "smaller oleKHPs", &chemOut);
		printSpeciesInFile(&detailedMechanism, oleOH, "oleOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, oleketo_radical, "oleketo_radical", &chemOut);
		printSpeciesInFile(&detailedMechanism, OHoleOOH, "OHoleOOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, OHoleKeto, "OHoleKeto", &chemOut);
		printSpeciesInFile(&detailedMechanism, oleROH, "oleROH", &chemOut);
		printSpeciesInFile(&detailedMechanism, dienesOH, "dienesOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, ROH, "ROH", &chemOut);
		printSpeciesInFile(&detailedMechanism, dienesR, "dienesR", &chemOut);
		printSpeciesInFile(&detailedMechanism, diketo, "diketo", &chemOut);
		printSpeciesInFile(&detailedMechanism, dienesCO, "dienesCO", &chemOut);
		printSpeciesInFile(&detailedMechanism, OHRO, "OHRO", &chemOut);
		printSpeciesInFile(&detailedMechanism, QOHOOH, "QOHOOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, OHketo, "OHketo", &chemOut);
		printSpeciesInFile(&detailedMechanism, OHKHP, "OHKHP", &chemOut);


		detailedMechanism << "" << std::endl;
		detailedMechanism << "END" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "REACTIONS" << std::endl;
		detailedMechanism << std::endl << "! Class01: Unimolecular decompositions:" << std::endl;
		printReactions(&detailedMechanism, initReac, &chemOut);
		detailedMechanism << std::endl << "! Class02: H abstraction:" << std::endl;
		printReactions(&detailedMechanism, hAbsReac, &chemOut);
		detailedMechanism << std::endl << "! Class03: Alkenyl radical isomerization:" << std::endl;
		printReactions(&detailedMechanism, RIsomReac, &chemOut);
		detailedMechanism << std::endl << "! Class04: Alkenyl radical beta scissions:" << std::endl;
		printReactions(&detailedMechanism, RBetaDecReac, &chemOut);
		detailedMechanism << std::endl << "! Class05: H addition to olefins" << std::endl;
		printReactions(&detailedMechanism, oleHadd, &chemOut);
		detailedMechanism << std::endl << "! Class06: O addition to olefins" << std::endl;
		printReactions(&detailedMechanism, OaddOle, &chemOut);
		detailedMechanism << std::endl << "! Class07: O2 addition to vinylic radicals" << std::endl;
		printReactions(&detailedMechanism, O2AddVinR, &chemOut);
		detailedMechanism << std::endl << "! Class08: HO2 addition to allyRs (RA+HO2=RAOOH)" << std::endl;
		printReactions(&detailedMechanism, allyHO2add_rec, &chemOut);
		detailedMechanism << std::endl << "! Class09: HO2 addition to allyRs (RA+HO2=alkROs+OH)" << std::endl;
		printReactions(&detailedMechanism, allyHO2add_dec, &chemOut);
		detailedMechanism << std::endl << "! Class10: RAOOH=RAO+OH" << std::endl;
		printReactions(&detailedMechanism, ROOHdecom, &chemOut);
		detailedMechanism << std::endl << "! Class11: RAO decomposition" << std::endl;
		printReactions(&detailedMechanism, ROdecom, &chemOut);
		detailedMechanism << std::endl << "! Class12: CH3 addition to allyRs (RACH3)" << std::endl;
		printReactions(&detailedMechanism, allyCH3add, &chemOut);
		//detailedMechanism << std::endl << "! ClassXX: Produced higherOle consumption" << std::endl;
		//printReactions(&detailedMechanism, Higheroleconsum, &chemOut);
		detailedMechanism << std::endl << "! Class13: HO2 addition to olefins" << std::endl;
		printReactions(&detailedMechanism, HO2addOle, &chemOut);
		detailedMechanism << std::endl << "! Class14: allyRs HAA by O2 reactions" << std::endl;
		printReactions(&detailedMechanism, O2absAllyR, &chemOut);
		detailedMechanism << std::endl << "! Class15: OH addition to olefins" << std::endl;
		printReactions(&detailedMechanism, oleOHadd, &chemOut);
		detailedMechanism << std::endl << "! Class16: betaROH isomerization" << std::endl;
		printReactions(&detailedMechanism, ROHisomReac, &chemOut);
		detailedMechanism << std::endl << "! Class17: O2 addition to ROHs" << std::endl;
		printReactions(&detailedMechanism, O2AddROH, &chemOut);
		detailedMechanism << std::endl << "! Class18: betaROHO2 Waddington mechanism" << std::endl;
		printReactions(&detailedMechanism, ROHO2Waddington, &chemOut);
		detailedMechanism << std::endl << "! Class19: Waddington product OROOH decompositions" << std::endl;
		printReactions(&detailedMechanism, Waddecom, &chemOut);
		detailedMechanism << std::endl << "! Class20: ROHO2 elimination reactions added" << std::endl;
		printReactions(&detailedMechanism, betaROHO2Elim, &chemOut);
		detailedMechanism << std::endl << "! Class21: ROHOO isomerization" << std::endl;
		printReactions(&detailedMechanism, ROHO2Isom, &chemOut);
		

		detailedMechanism << std::endl << "! Class22: beta-QOHOOH scissions" << std::endl;
		printReactions(&detailedMechanism, betaQOHOOHscis, &chemOut);
		detailedMechanism << std::endl << "! Class22_combine: QOHOOH decomposition reactions" << std::endl;
		printReactions(&detailedMechanism, QOHOOHdecom, &chemOut);

		detailedMechanism << std::endl << "! Class23: QOHOOH to cyclic ethers reactions" << std::endl;
		printReactions(&detailedMechanism, QOHOOHToEthers, &chemOut);
		detailedMechanism << std::endl << "! Class24: O2 addition to QOHOOH reactions" << std::endl;
		printReactions(&detailedMechanism, O2AddQOHOOH, &chemOut);
		detailedMechanism << std::endl << "! Class25: O2QOHOOH HO2 elimination reactions" << std::endl;
		printReactions(&detailedMechanism, O2QOHOOHElim, &chemOut);
		detailedMechanism << std::endl << "! Class26: OH-KHP formation from O2QOHOOH reactions" << std::endl;
		printReactions(&detailedMechanism, O2QOHOOH2KHP, &chemOut);
		detailedMechanism << std::endl << "! Class27: OH-KHP decomposition reactions reactions" << std::endl;
		printReactions(&detailedMechanism, OHKHPdecom, &chemOut);
		detailedMechanism << std::endl << "! Class28: O2QOHOOH isomerization reactions" << std::endl;
		printReactions(&detailedMechanism, O2QOHOOHIsom, &chemOut);
		detailedMechanism << std::endl << "! Class29: POHOOH decomposition reactions" << std::endl;
		printReactions(&detailedMechanism, POHOOHdecom, &chemOut);
		detailedMechanism << std::endl << "! Class30: POHOOH to ether-OHOOHs reactions" << std::endl;
		printReactions(&detailedMechanism, POHOOHTocE, &chemOut);
		detailedMechanism << std::endl << "! Class31: O2 addition to oleR:"<< std::endl;
		printReactions(&detailedMechanism, O2addToRReac, &chemOut);
		detailedMechanism << std::endl << "! Class32: oleROO HO2 elimination:" << std::endl;
		printReactions(&detailedMechanism, oleFromROOReac, &chemOut);
		detailedMechanism << std::endl << "! Class33: oleROO isomerization:" << std::endl;
		printReactions(&detailedMechanism, ROOIsomReac, &chemOut);
		detailedMechanism << std::endl << "! Class34: oleQOOH decomposition:" << std::endl;
		printReactions(&detailedMechanism, QOOHDecReac, &chemOut);
		detailedMechanism << std::endl << "! Class35: oleQOOH to cyclic ethers:" << std::endl;
		printReactions(&detailedMechanism, QOOHToCEthReac, &chemOut);
		detailedMechanism << std::endl << "! Class36: OleOHs (roho2 elim) consumption:" << std::endl;
		printReactions(&detailedMechanism, oleOHdecom, &chemOut);
		detailedMechanism << std::endl << "! Class37: O2 addition to oleQOOH:" << std::endl;
		printReactions(&detailedMechanism, O2addToQOOHReac, &chemOut);
		detailedMechanism << std::endl << "! Class38: oleKHP formation:" << std::endl;
		printReactions(&detailedMechanism, KHPFormReac, &chemOut);
		detailedMechanism << std::endl << "! Class39: oleKHP decomposition:" << std::endl;
		printReactions(&detailedMechanism, KHPDecReac, &chemOut);
		detailedMechanism << std::endl << "! Class40: OHOleOOHs (o2qohooh elim) consumption:" << std::endl;
		printReactions(&detailedMechanism, oleOHOOHdecom, &chemOut);
		detailedMechanism << std::endl << "! Class41: oleO2QOOH isomerization:" << std::endl;
		printReactions(&detailedMechanism, OOQOOHIsomReac, &chemOut);
		detailedMechanism << std::endl << "! Class42: oleP(OOH)2 decomposition:" << std::endl;
		printReactions(&detailedMechanism, POOH2DecReac, &chemOut);
		detailedMechanism << std::endl << "! Class43: P(OOH)2 to cyclic ethers-OOH:" << std::endl;
		printReactions(&detailedMechanism, POOH2ToCEthOOHReac, &chemOut);
		detailedMechanism << std::endl << "! Class44: oleO2QOOH elimination to dienes-OOH:" << std::endl;
		printReactions(&detailedMechanism, OOQOOHToOLEOOHReac, &chemOut);
		detailedMechanism << std::endl << "! Class45: olecEthOOHs (from oleP(OOH)2) consumption:" << std::endl;
		printReactions(&detailedMechanism, cEthOOHDecReac, &chemOut);
		detailedMechanism << std::endl << "! Class46: olecEths (from oleQOOH) consumption:" << std::endl;
		printReactions(&detailedMechanism, cEthDecReac, &chemOut);
		detailedMechanism << std::endl << "! Class47: OHcEthOOHs (from POHOOH2) consumption:" << std::endl;
		printReactions(&detailedMechanism, OHcEthOOHdec, &chemOut);
		detailedMechanism << std::endl << "! Class48: OHcEthDec (from QOHOOH) consumption:" << std::endl;
		printReactions(&detailedMechanism, OHcEthDec, &chemOut);
		detailedMechanism << std::endl << "! Class49: dieneOOHdecom (oleO2QOOH elim) consumption:" << std::endl;
		printReactions(&detailedMechanism, dieneOOHdecom, &chemOut);
		detailedMechanism << std::endl << "! Class50: alpha-QOHOOH To KHPs:" << std::endl;
		printReactions(&detailedMechanism, alphaqohooh, &chemOut);

		//detailedMechanism << std::endl << "! Cyclic ethers-OOH decomposition:" << std::endl;
		//printReactions(&detailedMechanism, cEthOOHDecReac, &chemOut);



		//if (reversibeDetailed == false)
		//{
		//detailedMechanism << std::endl << "! QOOH isomerization:" << std::endl;
		//printReactions(&detailedMechanism, QOOHIsomReac, &chemOut);
		//}

		//if (reversibeDetailed == false)
		//{
		//detailedMechanism << std::endl << "! Oxygen elimination from ROO:" << std::endl;
		//printReactions(&detailedMechanism, O2ElimFromROOReac, &chemOut);
		//}

		//if (reversibeDetailed == false)
		//{
		//detailedMechanism << std::endl << "! Oxygen elimination from OOQOOH:"
		//	<< std::endl;
		//printReactions(&detailedMechanism, O2ElimFromOOQOOHReac, &chemOut);
		//}


		//if (reversibeDetailed == false)
		//{
		//detailedMechanism << std::endl << "! P(OOH)2 isomerization:" << std::endl;
		//printReactions(&detailedMechanism, POOH2IsomReac, &chemOut);
		//}
		//detailedMechanism << std::endl << "! Olefins formation from R + O2:"
		//	<< std::endl;
		//printReactions(&detailedMechanism, oleFromRPlusO2reac, &chemOut);



		//detailedMechanism << std::endl << "! Olefins-OOH decomposition:" << std::endl;
		//printReactions(&detailedMechanism, OLEOOHDecReac, &chemOut);

		//detailedMechanism << std::endl << "! Cyclic ethers decomposition:" << std::endl;
		//printReactions(&detailedMechanism, cEthDecReac, &chemOut);
		//detailedMechanism << std::endl << "! Allylic radicals formation:" << std::endl;
		//printReactions(&detailedMechanism, allRadFormRac, &chemOut);
		//detailedMechanism << std::endl << "! Alkenyl RO formation:" << std::endl;
		//printReactions(&detailedMechanism, alkROFormReac, &chemOut);
		//detailedMechanism << std::endl << "! Alkenyl RO decomposition:" << std::endl;
		//printReactions(&detailedMechanism, alkRODecReac, &chemOut);
		detailedMechanism << "" << std::endl;
		detailedMechanism << "END" << std::endl;
		detailedMechanism.close();

		std::cout << " - Kinetic mechanism printed." << std::endl;

		// #########################################################################
		// ####################### PRINT THERMO FILE ###############################
		// #########################################################################

		std::string detailedThermoPath = subMechFold + "\\thermo_detailed.inp";
		std::ofstream detailedThermo(detailedThermoPath);
		detailedThermo << "! THERMODYNAMIC DATA FOR " << chemOut.molToName(HC)
			<< " SUBMECHANISM" << std::endl;
		detailedThermo << "! " << HC << std::endl;
		detailedThermo << std::endl;
		detailedThermo << "THERMO" << std::endl;
		detailedThermo << "300.000  1000.000  5000.000" << std::endl;
		detailedThermo << "AR                G         1  0    0      0G   200.000  6000.00  1000.00      1" << std::endl;
		detailedThermo << "+2.50000000E+00+0.00000000E+00+0.00000000E+00+0.00000000E+00+0.00000000E+00    2" << std::endl;
		detailedThermo << "-7.45375000E+02+4.37967491E+00+2.50000000E+00+0.00000000E+00+0.00000000E+00    3" << std::endl;
		detailedThermo << "+0.00000000E+00+0.00000000E+00-7.45375000E+02+4.37967491E+00+0.00000000E+00    4" << std::endl;
		detailedThermo << "N2                G         2    0    0    0G   200.000  6000.00  1000.00      1" << std::endl;
		detailedThermo << "+2.95257637E+00+1.39690040E-03-4.92631603E-07+7.86010195E-11-4.60755204E-15    2" << std::endl;
		detailedThermo << "-9.23948688E+02+5.87188762E+00+3.53100528E+00-1.23660988E-04-5.02999433E-07    3" << std::endl;
		detailedThermo << "+2.43530612E-09-1.40881235E-12-1.04697628E+03+2.96747038E+00+0.00000000E+00    4" << std::endl;
		detailedThermo << "HE                G        1    0    0    0 G   200.000  6000.00  1000.00      1" << std::endl;
		detailedThermo << "+2.50000000E+00+0.00000000E+00+0.00000000E+00+0.00000000E+00+0.00000000E+00    2" << std::endl;
		detailedThermo << "-7.45375000E+02+9.28723974E-01+2.50000000E+00+0.00000000E+00+0.00000000E+00    3" << std::endl;
		detailedThermo << "+0.00000000E+00+0.00000000E+00-7.45375000E+02+9.28723974E-01+0.00000000E+00    4" << std::endl;
		for (auto& spec : allSpecies)
		{
			if (spec == N2)
				continue;
			detailedThermo << thermOut.NASAOutput(&spec);
		}
		detailedThermo << "END";
		detailedThermo.close();
		std::cout << " - Thermodynamic data file printed." << std::endl;

		UTL::concatenateUnique(&totalReactionsList, &allReactions);
		UTL::concatenateUnique(&totalSpeciesList, &allSpecies);

		currentTime = std::chrono::steady_clock::now();
		timeFile << "Detailed generation for " << chemOut.molToName(HC) << " ended at " <<
			std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count() <<
			" s" << std::endl;

		// #########################################################################
		// ########################### LUMP MECHANISM ##############################
		// #########################################################################
		//if (generateLumped)
		//{
		//	UTL::printTitle("GENERATING LUMPED MECHANISM");
		//	// generate mechanisms
		//	std::string lumpedMechFilePath = subMechFold + "\\kinetic_lumped.inp";
		//	std::ofstream lumpedMech(lumpedMechFilePath);
		//	
		//	UTL::createDirectory("temp");
		//	UTL::createDirectory("temp/CKmech");
		//	std::ofstream kinCK("temp/CKmech/kin.txt");
		//	kinCK << "ELEMENTS" << std::endl;
		//	kinCK << "C" << std::endl << "H" << std::endl << "O" << std::endl
		//		<< "N" << std::endl << "AR" << std::endl << "HE" << std::endl;
		//	kinCK << "END" << std::endl << std::endl;
		//	kinCK << "SPECIES" << std::endl;
		//	kinCK << baseMech.speciesText.str() << std::endl;
		//	for (auto& name : allNames)
		//		if (UTL::isPresent(&baseMech.mechSpecies, name) == false)
		//			kinCK << name << std::endl;
		//	kinCK << "END" << std::endl << std::endl;
		//	kinCK << "REACTIONS" << std::endl;
		//	if (lumpWithCoreMech)
		//		kinCK << baseMech.reactionsText.str() << std::endl;
		//	for (auto& reac : allReactions)
		//	{
		//		if (isIncluded(reac, &(baseMech.reacs), &chemOut) == false)
		//			printReaction(&kinCK, reac, &chemOut);
		//	}
		//	kinCK << "END" << std::endl << std::endl;
		//	kinCK.close();
		//
		//	std::ofstream therCK("temp/CKmech/ther.txt");
		//	therCK << "THERMO" << std::endl;
		//	therCK << "300.000  1000.000  5000.000" << std::endl;
		//	therCK << baseMech.thermoText.str() << std::endl;
		//	for (auto& spec : allSpecies)
		//		if (spec.isSpecialMolecule() == 0)
		//			therCK << thermOut.NASAOutput(&spec);
		//	kinCK << "END" << std::endl;
		//	therCK.close();
		//
		//	// run the simulations
		//	//Simulation sim("temp/CKmech/kin.txt", "temp/CKmech/ther.txt", Temps, Press,
		//	//	&allSpecies, &chemOut, numCores);
		//	Simulation* sim;
		//
		//	if (equilibriumLumping)
		//	{
		//		sim = new Simulation(Temps, &allSpecies, &chemOut, &thermOut, numCores);
		//	}
		//	else 
		//	{
		//		sim = new Simulation("temp/CKmech/kin.txt", "temp/CKmech/ther.txt", Temps, Press,
		//			&allSpecies, &chemOut, numCores);
		//		if (useBatch)
		//		{
		//			sim->setBatchParameters(HC, eqRatio, 5000, MAX_SIMULATION_TIME);
		//		}
		//		else
		//		{
		//			sim->setCSTRParameters(HC, eqRatio, tau, 5000, MAX_SIMULATION_TIME);
		//		}
		//	}
		//
		//	sim->solve();
		//	std::string resultsDirPath = subMechFold + "\\simResults";
		//	UTL::createDirectory(resultsDirPath);
		//	sim->printResults(resultsDirPath);
		//	
		//	std::ofstream targetRatesFile(resultsDirPath + "\\targetRates.txt");
		//	// lump reactions
		//
		//	if (initReac.size() > 0)
		//	{
		//		LumpedReaction initReacLump(initReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(initReacLump);
		//		lumpedMech << initReacLump.print();
		//		targetRatesFile << initReacLump.printTargetRates();
		//	}
		//	if (hAbsReacByO2.size() > 0)
		//	{
		//		LumpedReaction hAbsReacByO2Lump(hAbsReacByO2, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(hAbsReacByO2Lump);
		//		lumpedMech << hAbsReacByO2Lump.print();
		//		targetRatesFile << hAbsReacByO2Lump.printTargetRates();
		//	}
		//	if (hAbsReacByOH.size() > 0)
		//	{
		//		LumpedReaction hAbsReacByOHLump(hAbsReacByOH, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(hAbsReacByOHLump);
		//		lumpedMech << hAbsReacByOHLump.print();
		//		targetRatesFile << hAbsReacByOHLump.printTargetRates();
		//	}
		//	if (hAbsReacByH.size() > 0)
		//	{
		//		LumpedReaction hAbsReacByHLump(hAbsReacByH, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(hAbsReacByHLump);
		//		lumpedMech << hAbsReacByHLump.print();
		//		targetRatesFile << hAbsReacByHLump.printTargetRates();
		//	}
		//	if (hAbsReacByO.size() > 0)
		//	{
		//		LumpedReaction hAbsReacByOLump(hAbsReacByO, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(hAbsReacByOLump);
		//		lumpedMech << hAbsReacByOLump.print();
		//		targetRatesFile << hAbsReacByOLump.printTargetRates();
		//	}
		//	if (hAbsReacByHO2.size() > 0)
		//	{
		//		LumpedReaction hAbsReacByHO2Lump(hAbsReacByHO2, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(hAbsReacByHO2Lump);
		//		lumpedMech << hAbsReacByHO2Lump.print();
		//		targetRatesFile << hAbsReacByHO2Lump.printTargetRates();
		//	}
		//	if (hAbsReacByCH3.size() > 0)
		//	{
		//		LumpedReaction hAbsReacByCH3Lump(hAbsReacByCH3, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(hAbsReacByCH3Lump);
		//		lumpedMech << hAbsReacByCH3Lump.print();
		//		targetRatesFile << hAbsReacByCH3Lump.printTargetRates();
		//	}
		//	// JIAXIN
		//	if (hAbsReacByCH3O.size() > 0)
		//	{
		//		LumpedReaction hAbsReacByCH3OLump(hAbsReacByCH3O, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(hAbsReacByCH3OLump);
		//		lumpedMech << hAbsReacByCH3OLump.print();
		//		targetRatesFile << hAbsReacByCH3OLump.printTargetRates();
		//	}
		//	if (hAbsReacByCH3O2.size() > 0)
		//	{
		//		LumpedReaction hAbsReacByCH3O2Lump(hAbsReacByCH3O2, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(hAbsReacByCH3O2Lump);
		//		lumpedMech << hAbsReacByCH3O2Lump.print();
		//		targetRatesFile << hAbsReacByCH3O2Lump.printTargetRates();
		//	}
		//	if (hAbsReacByC2H5.size() > 0)
		//	{
		//		LumpedReaction hAbsReacByC2H5Lump(hAbsReacByC2H5, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(hAbsReacByC2H5Lump);
		//		lumpedMech << hAbsReacByC2H5Lump.print();
		//		targetRatesFile << hAbsReacByC2H5Lump.printTargetRates();
		//	}
		//	
		//	if (O2addToRReac.size() > 0)
		//	{
		//		LumpedReaction O2addToRReacLump(O2addToRReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(O2addToRReacLump);
		//		lumpedMech << O2addToRReacLump.print();
		//		targetRatesFile << O2addToRReacLump.printTargetRates();
		//	}
		//	//if (O2ElimFromROOReac.size() > 0)
		//	//{
		//	//	LumpedReaction O2ElimFromROOReacLump(O2ElimFromROOReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(O2ElimFromROOReacLump);
		//	//	lumpedMech << O2ElimFromROOReacLump.print();
		//	//	targetRatesFile << O2ElimFromROOReacLump.printTargetRates();
		//	//}
		//	if (O2addToQOOHReac.size() > 0)
		//	{
		//		LumpedReaction O2addToQOOHReacLump(O2addToQOOHReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(O2addToQOOHReacLump);
		//		lumpedMech << O2addToQOOHReacLump.print();
		//		targetRatesFile << O2addToQOOHReacLump.printTargetRates();
		//	}
		//	//if (O2ElimFromOOQOOHReac.size() > 0)
		//	//{
		//	//	LumpedReaction O2ElimFromOOQOOHReacLump(O2ElimFromOOQOOHReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(O2ElimFromOOQOOHReacLump);
		//	//	lumpedMech << O2ElimFromOOQOOHReacLump.print();
		//	//	targetRatesFile << O2ElimFromOOQOOHReacLump.printTargetRates();
		//	//}
		//	if (ROOIsomReac.size() > 0)
		//	{
		//		LumpedReaction ROOIsomReacLump(ROOIsomReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(ROOIsomReacLump);
		//		lumpedMech << ROOIsomReacLump.print();
		//		targetRatesFile << ROOIsomReacLump.printTargetRates();
		//	}
		//	//if (QOOHIsomReac.size() > 0)
		//	//{
		//	//	LumpedReaction QOOHIsomReacLump(QOOHIsomReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(QOOHIsomReacLump);
		//	//	lumpedMech << QOOHIsomReacLump.print();
		//	//	targetRatesFile << QOOHIsomReacLump.printTargetRates();
		//	//}
		//	if (OOQOOHIsomReac.size() > 0)
		//	{
		//		LumpedReaction OOQOOHIsomReacLump(OOQOOHIsomReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(OOQOOHIsomReacLump);
		//		lumpedMech << OOQOOHIsomReacLump.print();
		//		targetRatesFile << OOQOOHIsomReacLump.printTargetRates();
		//	}
		//	//if (POOH2IsomReac.size() > 0)
		//	//{
		//	//	LumpedReaction POOH2IsomReacLump(POOH2IsomReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(POOH2IsomReacLump);
		//	//	lumpedMech << POOH2IsomReacLump.print();
		//	//	targetRatesFile << POOH2IsomReacLump.printTargetRates();
		//	//}
		//	if (oleFromROOReac.size() > 0)
		//	{
		//		LumpedReaction oleFromROOReacLump(oleFromROOReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(oleFromROOReacLump);
		//		lumpedMech << oleFromROOReacLump.print();
		//		targetRatesFile << oleFromROOReacLump.printTargetRates();
		//	}
		//	if (RBetaDecReac.size() > 0)
		//	{
		//		LumpedReaction RBetaDecReacLump(RBetaDecReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(RBetaDecReacLump);
		//		lumpedMech << RBetaDecReacLump.print();
		//		targetRatesFile << RBetaDecReacLump.printTargetRates();
		//	}
		//	//if (oleFromRPlusO2reac.size() > 0)
		//	//{
		//	//	LumpedReaction oleFromRPlusO2reacLump(oleFromRPlusO2reac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(oleFromRPlusO2reacLump);
		//	//	lumpedMech << oleFromRPlusO2reacLump.print();
		//	//	targetRatesFile << oleFromRPlusO2reacLump.printTargetRates();
		//	//}
		//	//if (betaQOOHDecReac.size() > 0)
		//	//{
		//	//	LumpedReaction betaQOOHDecReacLump(betaQOOHDecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(betaQOOHDecReacLump);
		//	//	lumpedMech << betaQOOHDecReacLump.print();
		//	//	targetRatesFile << betaQOOHDecReacLump.printTargetRates();
		//	//}
		//	//if (gammaQOOHDecReac.size() > 0)
		//	//{
		//	//	LumpedReaction gammaQOOHDecReacLump(gammaQOOHDecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(gammaQOOHDecReacLump);
		//	//	lumpedMech << gammaQOOHDecReacLump.print();
		//	//	targetRatesFile << gammaQOOHDecReacLump.printTargetRates();
		//	//}
		//	//if (deltaQOOHDecReac.size() > 0)
		//	//{
		//	//	LumpedReaction deltaQOOHDecReacLump(deltaQOOHDecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(deltaQOOHDecReacLump);
		//	//	lumpedMech << deltaQOOHDecReacLump.print();
		//	//	targetRatesFile << deltaQOOHDecReacLump.printTargetRates();
		//	//}
		//	if (QOOHToCEthReac.size() > 0)
		//	{
		//		LumpedReaction QOOHToCEthReacLump(QOOHToCEthReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(QOOHToCEthReacLump);
		//		lumpedMech << QOOHToCEthReacLump.print();
		//		targetRatesFile << QOOHToCEthReacLump.printTargetRates();
		//	}
		//	//if (POOH2ToCEthOOHReac.size() > 0)
		//	//{
		//	//	LumpedReaction POOH2ToCEthOOHReacLump(POOH2ToCEthOOHReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(POOH2ToCEthOOHReacLump);
		//	//	lumpedMech << POOH2ToCEthOOHReacLump.print();
		//	//	targetRatesFile << POOH2ToCEthOOHReacLump.printTargetRates();
		//	//}
		//	//if (betaPOOH2DecReac.size() > 0)
		//	//{
		//	//	LumpedReaction betaPOOH2DecReacLump(betaPOOH2DecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(betaPOOH2DecReacLump);
		//	//	lumpedMech << betaPOOH2DecReacLump.print();
		//	//	targetRatesFile << betaPOOH2DecReacLump.printTargetRates();
		//	//}
		//	//if (gammaPOOH2DecReac.size() > 0)
		//	//{
		//	//	LumpedReaction gammaPOOH2DecReacLump(gammaPOOH2DecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(gammaPOOH2DecReacLump);
		//	//	lumpedMech << gammaPOOH2DecReacLump.print();
		//	//	targetRatesFile << gammaPOOH2DecReacLump.printTargetRates();
		//	//}
		//	//if (deltaPOOH2DecReac.size() > 0)
		//	//{
		//	//	LumpedReaction deltaPOOH2DecReacLump(deltaPOOH2DecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(deltaPOOH2DecReacLump);
		//	//	lumpedMech << deltaPOOH2DecReacLump.print();
		//	//	targetRatesFile << deltaPOOH2DecReacLump.printTargetRates();
		//	//}
		//	//if (cEthOOHDecReac.size() > 0)
		//	//{
		//	//	LumpedReaction cEthOOHDecReacLump(cEthOOHDecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(cEthOOHDecReacLump);
		//	//	lumpedMech << cEthOOHDecReacLump.print();
		//	//	targetRatesFile << cEthOOHDecReacLump.printTargetRates();
		//	//}
		//	if (OOQOOHToOLEOOHReac.size() > 0)
		//	{
		//		LumpedReaction OOQOOHToOLEOOHReacLump(OOQOOHToOLEOOHReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(OOQOOHToOLEOOHReacLump);
		//		lumpedMech << OOQOOHToOLEOOHReacLump.print();
		//		targetRatesFile << OOQOOHToOLEOOHReacLump.printTargetRates();
		//	}
		//	//if (OLEOOHDecReac.size() > 0)
		//	//{
		//	//	LumpedReaction OLEOOHDecReacLump(OLEOOHDecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(OLEOOHDecReacLump);
		//	//	lumpedMech << OLEOOHDecReacLump.print();
		//	//	targetRatesFile << OLEOOHDecReacLump.printTargetRates();
		//	//}
		//	if (KHPFormReac.size() > 0)
		//	{
		//		LumpedReaction KHPFormReacLump(KHPFormReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(KHPFormReacLump);
		//		lumpedMech << KHPFormReacLump.print();
		//		targetRatesFile << KHPFormReacLump.printTargetRates();
		//	}
		//	//if (KHPDecReac.size() > 0)
		//	//{
		//	//	LumpedReaction KHPDecReacLump(KHPDecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(KHPDecReacLump);
		//	//	lumpedMech << KHPDecReacLump.print();
		//	//	targetRatesFile << KHPDecReacLump.printTargetRates();
		//	//}
		//	//if (cEthDecReac.size() > 0)
		//	//{
		//	//	LumpedReaction cEthDecReacLump(cEthDecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(cEthDecReacLump);
		//	//	lumpedMech << cEthDecReacLump.print();
		//	//	targetRatesFile << cEthDecReacLump.printTargetRates();
		//	//}
		//	//if (allRadFormRac.size() > 0)
		//	//{
		//	//	LumpedReaction allRadFormRacLump(allRadFormRac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(allRadFormRacLump);
		//	//	lumpedMech << allRadFormRacLump.print();
		//	//	targetRatesFile << allRadFormRacLump.printTargetRates();
		//	//}
		//	//if (alkROFormReac.size() > 0)
		//	//{
		//	//	LumpedReaction alkROFormReacLump(alkROFormReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(alkROFormReacLump);
		//	//	lumpedMech << alkROFormReacLump.print();
		//	//	targetRatesFile << alkROFormReacLump.printTargetRates();
		//	//}
		//	//if (alkRODecReac.size() > 0)
		//	//{
		//	//	LumpedReaction alkRODecReacLump(alkRODecReac, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(alkRODecReacLump);
		//	//	lumpedMech << alkRODecReacLump.print();
		//	//	targetRatesFile << alkRODecReacLump.printTargetRates();
		//	//}
		//	if (oleOHadd.size() > 0)
		//	{
		//		LumpedReaction oleOHaddLump(oleOHadd, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(oleOHaddLump);
		//		lumpedMech << oleOHaddLump.print();
		//		targetRatesFile << oleOHaddLump.printTargetRates();
		//	}
		//	if (O2AddROH.size() > 0)
		//	{
		//		LumpedReaction O2AddROHLump(O2AddROH, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(O2AddROHLump);
		//		lumpedMech << O2AddROHLump.print();
		//		targetRatesFile << O2AddROHLump.printTargetRates();
		//	}
		//	if (ROHisomReac.size() > 0)
		//	{
		//		LumpedReaction ROHisomReacLump(ROHisomReac, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(ROHisomReacLump);
		//		lumpedMech << ROHisomReacLump.print();
		//		targetRatesFile << ROHisomReacLump.printTargetRates();
		//	}
		//	if (ROHO2Isom.size() > 0)
		//	{
		//		LumpedReaction ROHO2IsomLump(ROHO2Isom, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(ROHO2IsomLump);
		//		lumpedMech << ROHO2IsomLump.print();
		//		targetRatesFile << ROHO2IsomLump.printTargetRates();
		//	}
		//
		//
		//	if (oleHadd.size() > 0)
		//	{
		//		LumpedReaction oleHaddLump(oleHadd, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(oleHaddLump);
		//		lumpedMech << oleHaddLump.print();
		//		targetRatesFile << oleHaddLump.printTargetRates();
		//	}
		//	if (allyHO2add_rec.size() > 0)
		//	{
		//		LumpedReaction allyHO2add_recLump(allyHO2add_rec, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(allyHO2add_recLump);
		//		lumpedMech << allyHO2add_recLump.print();
		//		targetRatesFile << allyHO2add_recLump.printTargetRates();
		//	}
		//	if (allyHO2add_dec.size() > 0)
		//	{
		//		LumpedReaction allyHO2add_decLump(allyHO2add_dec, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(allyHO2add_decLump);
		//		lumpedMech << allyHO2add_decLump.print();
		//		targetRatesFile << allyHO2add_decLump.printTargetRates();
		//	}
		//	if (allyCH3add.size() > 0)
		//	{
		//		LumpedReaction allyCH3addLump(allyCH3add, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(allyCH3addLump);
		//		lumpedMech << allyCH3addLump.print();
		//		targetRatesFile << allyCH3addLump.printTargetRates();
		//	}
		//	if (HO2addOle.size() > 0)
		//	{
		//		LumpedReaction HO2addOleLump(HO2addOle, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(HO2addOleLump);
		//		lumpedMech << HO2addOleLump.print();
		//		targetRatesFile << HO2addOleLump.printTargetRates();
		//	}
		//
		//	if (OaddOle.size() > 0)
		//	{
		//		LumpedReaction OaddOleLump(OaddOle, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(OaddOleLump);
		//		lumpedMech << OaddOleLump.print();
		//		targetRatesFile << OaddOleLump.printTargetRates();
		//	}
		//	if (ROOHdecom.size() > 0)
		//	{
		//		LumpedReaction ROOHdecomLump(ROOHdecom, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(ROOHdecomLump);
		//		lumpedMech << ROOHdecomLump.print();
		//		targetRatesFile << ROOHdecomLump.printTargetRates();
		//	}
		//	if (ROdecom.size() > 0)
		//	{
		//		LumpedReaction ROdecomLump(ROdecom, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(ROdecomLump);
		//		lumpedMech << ROdecomLump.print();
		//		targetRatesFile << ROdecomLump.printTargetRates();
		//	}
		//	//if (Higheroleconsum.size() > 0)
		//	//{
		//	//	LumpedReaction HigheroleconsumLump(Higheroleconsum, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(HigheroleconsumLump);
		//	//	lumpedMech << HigheroleconsumLump.print();
		//	//	targetRatesFile << HigheroleconsumLump.printTargetRates();
		//	//}
		//	//if (O2AddVinR.size() > 0)
		//	//{
		//	//	LumpedReaction O2AddVinRLump(O2AddVinR, sim, &thermOut, &chemOut);
		//	//	totalLumpedReactionsList.push_back(O2AddVinRLump);
		//	//	lumpedMech << O2AddVinRLump.print();
		//	//	targetRatesFile << O2AddVinRLump.printTargetRates();
		//	//}
		//	if (O2absAllyR.size() > 0)
		//	{
		//		LumpedReaction O2absAllyRLump(O2absAllyR, sim, &thermOut, &chemOut);
		//		totalLumpedReactionsList.push_back(O2absAllyRLump);
		//		lumpedMech << O2absAllyRLump.print();
		//		targetRatesFile << O2absAllyRLump.printTargetRates();
		//	}
		//
		//	targetRatesFile.close();
		//	
		//	lumpedMech.close();
		//	currentTime = std::chrono::steady_clock::now();
		//	timeFile << "Lumped generation for " << chemOut.molToName(HC) << " ended at " <<
		//		std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count() <<
		//		" s" << std::endl;
		//}
	}

	// #########################################################################
	// ######################## FILL MISSING REACTIONS #########################
	// #########################################################################
	//std::vector<bool> isSpeciesNonDecomposing(totalSpeciesList.size(), false);

	// total species list clean.
	std::vector<std::string> totalSpecInChIList(totalSpeciesList.size());

	for (int i = 0; i < totalSpeciesList.size(); i++)
	{
		if (!totalSpeciesList[i].isSpecialMolecule())
		{
			totalSpecInChIList[i] = totalSpeciesList[i].inchiName();
		}
		else
		{
			totalSpecInChIList[i] = std::to_string(i);
		}

		//std::cout << "totalSpecInChIList[i]:" << totalSpecInChIList[i] << std::endl;
		
	}

	//totalSpecInChIList.push_back("InChI=1S/C2H2/c1-2/h1-2H");
	std::unordered_map<std::string, size_t> firstpos;
	for (size_t i = 0; i < totalSpecInChIList.size(); ++i) {
		const std::string& species = totalSpecInChIList[i];

		if (firstpos.find(species) != firstpos.end())
		{
			// index of first appear if possible;
			size_t firstIndex = firstpos[species];
			// allNames[firstIndex] to allNames[i]
			totalSpeciesList[i] = totalSpeciesList[firstIndex];
		}
		else {
			// index of first appear
			firstpos[species] = i;
		}
	}


	for (int i = 0; i < totalReactionsList.size(); i++)
	{
		//if (!totalSpeciesList[i].isSpecialMolecule())
		//{
		//totalReactionsList[i] = totalReactionsList[i].reac_clean();
		std::vector<Molecola> reacs = totalReactionsList[i].Reactants();
		//std::vector<Molecola&> reac = reacs;
		for (int j = 0; j < reacs.size(); j++)
			if (!reacs[j].isSpecialMolecule())
			{
				std::string inchi = reacs[j].inchiName();
				int posr = find_pos(totalSpecInChIList, inchi);
				//std::cout << "posr=" << totalSpeciesList[posr] << std::endl;
				totalReactionsList[i].Reactants()[j] = totalSpeciesList[posr];
			}
		std::vector<Molecola> prods = totalReactionsList[i].Products();
		for (int z = 0; z < prods.size(); z++)
			if (!prods[z].isSpecialMolecule())
			{
				std::string inchi1 = prods[z].inchiName();
				
				std::cout << inchi1 << std::endl;
				int posp = find_pos(totalSpecInChIList, inchi1);
				//if (inchi1 == "InChI=1S/C6H9O/c1-2-3-4-5-6-7/h3-6H,2H2,1H3")
				//{
					//std::cout << "posp=" << posp << std::endl;
					//std::cout << "speciesname=" << totalSpeciesList[posp] << std::endl;
					////std::cout << "pos1=" << totalReactionsList[i] << std::endl;
					//std::cout << "prodz=" << totalReactionsList[i].Products()[z] << std::endl;
				if (posp != -1)
					totalReactionsList[i].Products()[z] = totalSpeciesList[posp];
					//std::cout << "prodz1=" << totalReactionsList[i].Products()[z] << std::endl;
					//std::cout << " " << std::endl;
					//std::cout << "pos2=" << totalReactionsList[i] << std::endl;
				//}
			}
	}
	// total reactions list clean. 
	//totalReactionsList
	//for (int i = 0; i < totalReactionsList.size(); i++)
	//{
	//	//if (!totalSpeciesList[i].isSpecialMolecule())
	//	//{
	//	//totalReactionsList[i] = totalReactionsList[i].reac_clean();
	//	std::vector<Molecola*> reacs = totalReactionsList[i].reactantList();
	//	//std::vector<Molecola&> reac = reacs;
	//	for (int j = 0; j < reacs.size(); j++)
	//		if (!reacs[j]->isSpecialMolecule())
	//		{
	//			std::string inchi = reacs[j]->inchiName();
	//			int posr = find_pos(totalSpecInChIList, inchi);
	//			//std::cout << "posr=" << totalSpeciesList[posr] << std::endl;
	//			totalReactionsList[i].reactantList()[j] = &totalSpeciesList[posr];
	//		}
	//	std::vector<Molecola*> prods = totalReactionsList[i].productList();
	//	for (int z = 0; z < prods.size(); z++)
	//		if (!prods[z]->isSpecialMolecule())
	//		{
	//			std::string inchi1 = prods[z]->inchiName();
	//			int posp = find_pos(totalSpecInChIList, inchi1);
	//			if (inchi1 == "InChI=1S/C6H9O/c1-2-3-4-5-6-7/h3-6H,2H2,1H3")
	//			{
	//				std::cout << "posp=" << posp << std::endl;
	//				std::cout << "speciesname=" << totalSpeciesList[posp] << std::endl;
	//				std::cout << "pos1=" << totalReactionsList[i].productList()[z] << std::endl;
	//				totalReactionsList[i].productList().at(z) = &totalSpeciesList[posp];
	//				std::cout << "pos2=" << totalReactionsList[i] << std::endl;
	//			}
	//		}
	//
	//	//}
	//	//else
	//	//{
	//	//	totalSpecInChIList[i] = std::to_string(i);
	//	//}
	//}



	std::cout << std::endl;
	UTL::printEmbeddedString('#', "PROCESSING COMPLETE MECHANISM");

	int numOfNonConsumingSpecies = 0;
	int initialReacNum = totalReactionsList.size();
	int initialSpeciesNum = totalSpeciesList.size();
	std::vector<Reaction> newlyAddedReactions;
	// JIAXIN
	//for (int i = 0; i < totalSpeciesList.size(); i++)
	//{
	//	Molecola spec = totalSpeciesList[i];
	//	//std::cout << "spec: " << spec << std::endl;
	//	if (!isThereDecompositionPath(spec, &(baseMech.reacs), &totalReactionsList,
	//		&chemOut))
	//	{
	//		numOfNonConsumingSpecies++;
	//		std::vector<Reaction> decReac;
	//		std::vector<Reaction> tempReac;
	//		std::vector<Molecola> reactantVector = { spec };
	//		species typ = spec.kindOfSPecies();
	//		switch (typ)
	//		{
	//		case R_:
	//			decReac = RBetaDecompositioReactions(reactantVector, &k);
	//			tempReac = olefinsFromRPlusO2Reactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			tempReac = O2AdditionToRReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			break;
	//		case ROO_:
	//			decReac = olefinsFromROOReactions(reactantVector, &k);
	//			tempReac = ROOIsomerizationReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			tempReac = O2EliminationFromROOReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			break;
	//		case QOOH_:
	//			decReac = QOOHDecompositionReaction(reactantVector, &k);
	//			tempReac = QOOHDecompositionReaction(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			tempReac = QOOHToCEthReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			tempReac = O2AdditionToQOOHReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			tempReac = QOOHIsomerizationReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			break;
	//		case OOQOOH_:
	//			decReac = KHPFormationReactions(reactantVector, &k);
	//			tempReac = OOQOOHToOLEOOHReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			tempReac = O2EliminationFromOOQOOHReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			tempReac = OOQOOHIsomerizationReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			break;
	//		case OLE_:
	//			decReac = allylicRadicalsFormationReactions(reactantVector, &k);
	//			break;
	//		case CO_:
	//			if (spec.isAldehyde())
	//				decReac = aldehydesDecompositionReactions(reactantVector, &k);
	//			else if (spec.isKetone())
	//			{
	//				decReac = ketonesDecompositionReactions(reactantVector, &k);
	//				//std::cout << "Keto found: " << std::endl;
	//				//for (auto& reac : decReac)
	//				//	std::cout << "   - " << reac << std::endl;
	//			}
	//			break;
	//		case cEth_:
	//			decReac = cyclicEthersDecompositionReactions(reactantVector, &k);
	//			break;
	//		case RO_:
	//			UTL::warning("RO found");
	//			break;
	//		case KHP_:
	//			decReac = KHPDecompositionReactions(reactantVector, &k);
	//			break;
	//		case POOH2_:
	//			decReac = POOH2DecompositionReaction(reactantVector, &k);
	//			tempReac = POOH2DecompositionReaction(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			tempReac = POOH2ToCEthOOHReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			tempReac = POOH2IsomerizationReactions(reactantVector, &k);
	//			UTL::concatenate(&decReac, &tempReac);
	//			break;
	//		case oleOOH_:
	//			decReac = OLEOOHDecompositionReactions(reactantVector, &k);
	//			break;
	//		case cEthOOH_:
	//			decReac = cyclicEtherOOHDecompositionReactions(reactantVector, &k);
	//			break;
	//		case oleR_:
	//			decReac = alkenylROFormationReactions(reactantVector, &k);
	//			break;
	//		case oleCO_:
	//			if (spec.isOleAldehyde()) 
	//				decReac = aldehydeOlefinsDecompositionReactions(reactantVector, &k);
	//			else if (spec.isOleKetone())
	//			{
	//				decReac = ketonesOlefinsDecompositionReactions(reactantVector, &k);
	//			}
	//			break;
	//		case cEthR_:
	//			UTL::warning("cEthR found");
	//			break;
	//		case cEthCO_:
	//			UTL::warning("cEthCO found");
	//			break;
	//		case alkRO_:
	//			decReac = alkenylRODecompositionReactions(reactantVector, &k);
	//			break;
	//		case special_:
	//			UTL::warning("Special found");
	//			break;
	//		default:
	//			UTL::error("unexpected species found");
	//			break;
	//		}
	//		processDuplicateReactions(&decReac);
	//		//for (auto& reac : decReac)
	//		//	reac.setReversible(true);
	//		UTL::concatenate(&totalReactionsList, &decReac);
	//		addNewSpecies(&totalSpeciesList, &decReac);
	//		UTL::concatenate(&newlyAddedReactions, &decReac);
	//	}
	//}
	std::cout << " - Fixed " << numOfNonConsumingSpecies
		<< " species without consumption pathways." << std::endl;
	std::cout << "      " << (totalReactionsList.size() - initialReacNum)
		<< " reactions added" << std::endl;
	//std::cout << totalReactionsList.size() << std::endl;
	//std::cout << initialReacNum << std::endl;
	std::cout << "      " << (totalSpeciesList.size() - initialSpeciesNum)
		<< " species added" << std::endl;
	//std::cout << totalSpeciesList.size() << std::endl;
	//std::cout << initialSpeciesNum << std::endl;

	// #########################################################################
	// ################# PRINT COMPLETE DETAILED MECHANISM #####################
	// #########################################################################
	std::ofstream fullDetailedMech(outFolderPath + "\\completeDetailedMech.inp");
	fullDetailedMech << "ELEMENTS" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "C" << std::endl;
	fullDetailedMech << "H" << std::endl;
	fullDetailedMech << "N" << std::endl;
	fullDetailedMech << "O" << std::endl;
	fullDetailedMech << "AR" << std::endl;
	fullDetailedMech << "HE" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "END" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "SPECIES" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "" << std::endl;
	std::vector<std::string> baseMechSpeciesNames = { "O2" , "H2", "H2O2",
		"H2O", "N2", "AR", "HE", "CO", "CO2", "HO2", "OH", "H" };
	UTL::concatenateUnique(&baseMechSpeciesNames, &(baseMech.mechSpecies));
	printSpeciesInFile(&fullDetailedMech, baseMechSpeciesNames,
		"Base mechanism");
	//detect all the parent fuels present in the species
	std::vector<Molecola> allParentFules = fuels;
	//for (auto& spec : totalSpeciesList)
	//{
	//	if (spec.kindOfSPecies() == dienes_ || spec.kindOfSPecies() == dieneOOH_ || spec.kindOfSPecies() == dienesOH_
	//		|| spec.kindOfSPecies() == dienesR_ || spec.kindOfSPecies() == dienesCO_ || spec.size() < 5 || !spec.isLinear())
	//		continue;
	//	if (spec.isSpecialMolecule() == 0 && spec.kindOfSPecies() == OLE_)
	//		UTL::addUnique(&allParentFules, spec.parentFuelOH());
	//}
	std::vector<std::string> classesNames = {
	"R","ROO","QOOH","OOQOOH","OLE",
	"CO","cEth","RO","KHP","POOH2",
	"oleOOH","cEthOOH","oleR","oleCO","cEthR",
	"cEthCO","alkRO","oleKHP","dienes","oleOH",
	"oleRO","OHoleOOH","OHoleketo","oleROH","dienesOH",
	"ROH","dienesR","diketo","dienesCO","OHRO",
	"QOHOOH","COHCHO","OHKHP","ROHOO", "OROOHs",
	"OHcEth", "O2QOHOOH", "POHOOH", "OHcEthOOH", "olecEth",
	"olePOOH2", "dieneOOH", "olecEthOOH", "oleROO", "oleQOOH",
	"oleO2QOOH", "olecEthR", "oleROOH", "oleQO", "OHketo", "OHdiketo"};

	std::vector<species> classes =
	{
	R_, ROO_, QOOH_, OOQOOH_, OLE_,
	CO_, cEth_, RO_, KHP_, POOH2_,
	oleOOH_, cEthOOH_, oleR_, oleCO_, cEthR_, 
	cEthCO_, alkRO_, oleKHP_, dienes_, oleOH_,
	oleRO_, OHoleOOH_, OHoleketo_, oleROH_, dienesOH_,
	ROH_, dienesR_,diketo_,dienesCO_,OHRO_,
	QOHOOH_, COHCHO_, OHKHP_, ROHOO_, OROOHs_,
	OHcEth_, O2QOHOOH_, POHOOH_, OHcEthOOH_, olecEth_, 
	olePOOH2_, dieneOOH_, olecEthOOH_, oleROO_, oleQOOH_, 
	oleO2QOOH_, olecEthR_, oleROOH_, oleQO_, OHketo_, OHdiketo_
	};

	//for (auto& parFuel : allParentFules)
	//{
		//std::string fuelName = chemOut.molToName(parFuel);
		//fullDetailedMech << "!";
		//for (int i = 0; i < (76 - fuelName.size()) / 2; i++)
		//	fullDetailedMech << "#";
		//fullDetailedMech << " " << fuelName << " ";
		//for (int i = 0; i < (76 - fuelName.size()) / 2; i++)
		//	fullDetailedMech << "#";
		//fullDetailedMech << std::endl;
		//if (UTL::isPresent(&totalSpeciesList, parFuel)
		//	&& !UTL::isPresent(&baseMechSpeciesNames, chemOut.molToName(parFuel)))
		//	printSpeciesInFile(&fullDetailedMech, std::vector<Molecola> {parFuel},
		//		"fuel", & chemOut);
		for (int i = 0; i < classes.size(); i++)
		{
			std::vector<std::string> tempList;
			for (auto& spec : totalSpeciesList)
			{
				if (spec.isSpecialMolecule() != 0)
					continue;
				if (spec.kindOfSPecies() != classes[i])
					continue;
				// JIAXIN
				//if (!(spec.parentFuel() == parFuel))
				//	continue;
				//if (spec.size() != parFuel.size())
				//	continue;
				std::string specName = chemOut.molToName(spec);
				if (specName == "C2H2")
					std::cout << "hereeee" << std::endl;
				if (!UTL::isPresent(&baseMechSpeciesNames, specName))
					tempList.push_back(specName);
			}
			if (tempList.size() > 0)
			{
				printSpeciesInFile(&fullDetailedMech, tempList,
					classesNames[i]);
			}
		}
	//}

	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "END" << std::endl;
	fullDetailedMech << "" << std::endl;

	// print reactions
	fullDetailedMech << "REACTIONS" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "! BASE MECHANISM START" << std::endl;
	fullDetailedMech << baseMech.reactionsText.str() << std::endl;
	fullDetailedMech << "! BASE MECHANISM END" << std::endl;
	fullDetailedMech << "" << std::endl;

	std::vector<std::string> reactionsLabels = {
		"Class01: Unimolecular decompositions",
        "Class02: H abstraction",
        "Class03: Alkenyl radical isomerization",
        "Class04: Alkenyl radical beta-scissions",
        "Class05: H addition to olefins",
        "Class06: O addition to olefins",
        "Class07: O2 addition to vinylic radicals",
        "Class08: HO2 addition to allyRs (RA+HO2=RAOOH)",
        "Class09: HO2 addition to allyRs (RA+HO2=alkROs+OH)",
        "Class10: RAOOH=RAO+OH",
        "Class11: RAO decomposition",
        "Class12: CH3 addition to allyRs (RACH3)",
        "Class13: HO2 addition to olefins",
        "Class14: allyRs HAA by O2 reactions",
        "Class15: OH addition to olefins",
        "Class16: betaROH isomerization",
        "Class17: O2 addition to ROHs",
        "Class18: betaROHO2 Waddington mechanism",
        "Class19: Waddington product OROOH decompositions",
        "Class20: ROHO2 elimination reactions added",
        "Class21: ROHOO isomerization",
        "Class22: beta-QOHOOH scissions",
        "Class22_combine: QOHOOH decomposition reactions",
        "Class23: QOHOOH to cyclic ethers reactions",
        "Class24: O2 addition to QOHOOH reactions",
        "Class25: O2QOHOOH HO2 elimination reactions",
        "Class26: OH-KHP formation from O2QOHOOH reactions",
        "Class27: OH-KHP decomposition reactions",
        "Class28: O2QOHOOH isomerization reactions",
        "Class29: POHOOH decomposition reactions",
        "Class30: POHOOH to ether-OHOOHs reactions",
        "Class31: O2 addition to oleR",
        "Class32: oleROO HO2 elimination",
        "Class33: oleROO isomerization",
        "Class34: oleQOOH decomposition",
        "Class35: oleQOOH to cyclic ethers",
        "Class36: OleOHs (roho2elim) consumption",
        "Class37: O2 addition to oleQOOH",
        "Class38: oleKHP formation",
        "Class39: oleKHP decomposition",
        "Class40: OHOleOOHs (o2qohoohelim) consumption",
        "Class41: oleO2QOOH isomerization",
        "Class42: oleP(OOH)2 decomposition",
        "Class43: P(OOH)2 to cyclic ethers-OOH",
        "Class44: oleO2QOOH elimination to dienes-OOH",
        "Class45: olecEthOOHs (fromoleP(OOH)2) consumption",
        "Class46: olecEths (fromoleQOOH) consumption",
        "Class47: OHcEthOOHs (fromPOHOOH2) consumption",
        "Class48: OHcEthDec (fromQOHOOH) consumption",
        "Class49: dieneOOHdecom (oleO2QOOHelim) consumption",
		"Class50: alpha-QOHOOH To KHPs",
		//"Initiation reactions",
		//"H abstraction",
		//"O2 addition to R",
		//"O2 elimination from ROO",
		//"O2 addition to QOOH",
		//"O2 elimination from OOQOOH",
		//"Isomerization R",
		//"Isomerization ROO",
		//"Isomerization QOOH",
		//"Isomerization OOQOOH",
		//"Isomerization P(OOH)2",
		//"ROO to olefins",
		//"Radicals beta decomposition",
		//"Olefins from R + O2",
		//"Decomposition beta-QOOH",
		//"Decomposition gamma-QOOH",
		//"Decomposition delta-QOOH",
		//"Cyclic ethers from QOOH",
		//"Ethers-OOH from P(OOH)2",
		//"Beta-P(OOH)2 decomposition",
		//"Gamma-P(OOH)2 decomposition",
		//"Delta P(OOH)2 decomposition",
		//"Ether-OOH decomposition",
		//"Olefins-OOH from OOQOOH",
		//"Ole-OOH decomposition",
		//"OOQOOH conversion to ketohydroperoxide",
		//"Ketohydroperoxides decomposition",
		//"Cyclic ethers decomposition",
		//"Allylic radicals formation",
		//"Alkenyl RO formation",
		//"Alkenyl RO decomposition",
		//"Aldehydes decomposition",
		//"Olefin aldehydes decomposition",
		//"Ketones decomposition",
		//"Olefins ketones decomposition"
	};
	//for (auto& parFuel : allParentFules)
	//{
	//	std::string fuelName = chemOut.molToName(parFuel);
	//	fullDetailedMech << "!";
	//	for (int i = 0; i < (76 - fuelName.size()) / 2; i++)
	//		fullDetailedMech << "#";
	//	fullDetailedMech << " " << fuelName << " ";
	//	for (int i = 0; i < (76 - fuelName.size()) / 2; i++)
	//		fullDetailedMech << "#";
	//	fullDetailedMech << std::endl;
		for (auto& label : reactionsLabels)
		{
			std::vector<Reaction> tempList;
			for (auto& reac : totalReactionsList)
			{
				if (reac.reactionLabel() == label)
					//if (reac.parentFuel_complete() == parFuel.parentFuel_complete()) // JIAXIN
					//if (reac.reac_size() == parFuel.size()) // JIAXIN
						//if (isIncluded(reac, &(baseMech.reacs), &chemOut) == false)
							tempList.push_back(reac);
			}
			if (tempList.size() > 0)
			{
				fullDetailedMech << "! " << label << std::endl;
				printReactions(&fullDetailedMech, tempList, &chemOut);
				fullDetailedMech << std::endl;
			}
		}
	//}
	//for (auto& reac : totalReactionsList)
	//	if(isIncluded(reac, &(baseMech.reacs), &chemOut) == false)
	//		printReaction(&fullDetailedMech, reac, &chemOut);
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "END" << std::endl;
	fullDetailedMech.close();

	// print thermo file
	std::ofstream fullDetailedThermo(outFolderPath + "\\completeDetailedThermo.inp");
	fullDetailedThermo << "THERMO" << std::endl;
	fullDetailedThermo << "300.000  1000.000  5000.000" << std::endl;
	fullDetailedThermo << "! BASE MECHANISM START" << std::endl;
	fullDetailedThermo << baseMech.thermoText.str() << std::endl;
	fullDetailedThermo << "! BASE MECHANISM END" << std::endl;
	for (auto& spec : totalSpeciesList)
	{
		if (spec == N2)
			continue;
		fullDetailedThermo << thermOut.NASAOutput(&spec);
	}
	fullDetailedThermo << "END";
	fullDetailedThermo.close();

	// PRINT SPECIES LIST
	std::ofstream speciesListFile(outFolderPath + "\\speciesList.csv");
	speciesListFile << "Mech. name,structure,InChI,class,LumpedSpecies" << std::endl;
	std::vector<std::string> totalSpeciesNamesList(totalSpeciesList.size(), "");
	std::vector<std::string> totalSpeciesInChIsList(totalSpeciesList.size(), "");
	for (int i = 0; i < totalSpeciesList.size(); i++)
	{
		//if (totalSpeciesList[i].size() > 4)
		// JIAXIN
		if (totalSpeciesList[i].size() > 0)
		{
			totalSpeciesNamesList[i] = chemOut.molToName(totalSpeciesList[i]);
			totalSpeciesInChIsList[i] = totalSpeciesList[i].inchiName();
			speciesListFile << "\"" << totalSpeciesNamesList[i] << "\","
				<< totalSpeciesList[i] << ",\"" <<
				totalSpeciesInChIsList[i] << "\"," <<
				speciesToText(totalSpeciesList[i].kindOfSPecies()) << std::endl;
				// JIAXIN //speciesToText(totalSpeciesList[i].kindOfSPecies()) << "," << chemOut.molToNameLump(totalSpeciesList[i]) << std::endl;
		}
	}

	for (int i = 0; i < totalSpeciesNamesList.size(); i++)
	{
		if (totalSpeciesNamesList[i] == "")
			continue;
		for (int j = i + 1; j < totalSpeciesNamesList.size(); j++)
		{
			if (totalSpeciesNamesList[j] == "")
				continue;
			if (totalSpeciesNamesList[i] == totalSpeciesNamesList[j])
			{
				std::stringstream warningMessage;
				warningMessage << "Two different species found having the same naming: "
					<< std::endl << "    - " << totalSpeciesNamesList[i] << "  " <<
					totalSpeciesList[i] << "   " << totalSpeciesInChIsList[i]
					<< std::endl << "    - " << totalSpeciesNamesList[j] << "  " <<
					totalSpeciesList[j] << "   " << totalSpeciesInChIsList[j] << std::endl;
				warningMessage << "    Probably there is a conflict in the glossary file" << std::endl;
				warningMessage << "    Either remove definition for " << totalSpeciesNamesList[i]
					<< " or add definition for missing InChI." << std::endl << std::endl;
				UTL::warning(warningMessage.str());
			}
		}
	}
	speciesListFile.close();

	chemOut.printLongNameSpeciesMessage();

	// #########################################################################
	// ####################### COMPLETE LUMPED MECHANISM #######################
	// #########################################################################
	if (generateLumped)
	{
		//newlyAddedReactions
		for (auto& parFuel : allParentFules)
		{
			std::vector<Reaction> pertinentReaction;
			for (auto& R : newlyAddedReactions)
			{
				std::vector<Molecola*> reacs = R.reactantList();
				int maxSize = 0;
				Molecola relReac;
				for (auto& reac : reacs)
				{
					if (reac->isSpecialMolecule() == 0 && reac->size() > maxSize)
					{
						maxSize = reac->size();
						relReac = *reac;
					}
				}
				if (relReac.parentFuel() == parFuel)
					pertinentReaction.push_back(R);
			}
			std::vector<std::vector<Reaction>> pertReacByType(reactionsLabels.size());
			for (int i = 0; i < reactionsLabels.size(); i++)
			{
				for (auto& R : pertinentReaction)
					if (R.reactionLabel() == reactionsLabels[i])
						pertReacByType[i].push_back(R);
			}
			for (auto& reacVec : pertReacByType)
			{
				if (reacVec.size() == 0)
					continue;
				LumpedReaction lump(reacVec, &thermOut, Temps, &chemOut);
				totalLumpedReactionsList.push_back(lump);
			}
		}

		// generate species list
		std::vector<std::string> allLumpedNames;
		std::vector<Molecola>    allLumpedMols;
		for (auto& spec : totalSpeciesList)
		{
			if (UTL::addUnique(&allLumpedNames, chemOut.molToNameLump(spec)) == -1)
				allLumpedMols.push_back(spec);
		}

		std::vector<std::string> lumpedSpeciesForMech = { "O2" , "H2", "H2O2",
			"H2O", "N2", "AR", "HE", "CO", "CO2", "HO2", "OH", "H" };
		UTL::concatenateUnique(&lumpedSpeciesForMech, &(baseMech.mechSpecies));
		UTL::concatenateUnique(&lumpedSpeciesForMech, &allLumpedNames);

		std::ofstream fullLumpedMech(outFolderPath + "\\completeLumpedMech.inp");
		fullLumpedMech << "ELEMENTS" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "C" << std::endl;
		fullLumpedMech << "H" << std::endl;
		fullLumpedMech << "N" << std::endl;
		fullLumpedMech << "O" << std::endl;
		fullLumpedMech << "AR" << std::endl;
		fullLumpedMech << "HE" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "END" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "SPECIES" << std::endl;
		fullLumpedMech << "" << std::endl;
		printSpeciesInFile(&fullLumpedMech, lumpedSpeciesForMech, "all");
		fullLumpedMech << "END" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "REACTIONS" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "! BASE MECHANISM START" << std::endl;
		fullLumpedMech << baseMech.reactionsText.str() << std::endl;
		fullLumpedMech << "! BASE MECHANISM END" << std::endl;
		fullLumpedMech << "" << std::endl;
		for (auto& R : totalLumpedReactionsList)
			fullLumpedMech << R.print() << std::endl;
		fullLumpedMech << "END" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech.close();

		// print thermo file
		std::ofstream fullLumpedThermo(outFolderPath + "\\completeLumpedThermo.inp");
		fullLumpedThermo << "THERMO" << std::endl;
		fullLumpedThermo << "300.000  1000.000  5000.000" << std::endl;
		fullLumpedThermo << "! BASE MECHANISM START" << std::endl;
		fullLumpedThermo << baseMech.thermoText.str() << std::endl;
		fullLumpedThermo << "! BASE MECHANISM END" << std::endl;
		std::vector<std::string> addedSpecies;
		for (auto& reac : totalLumpedReactionsList)
			if (UTL::addUnique(&addedSpecies, reac.relevantReactantName) == -1)
				fullLumpedThermo << reac.relevantReactantThermo;
		for (int i = 0; i < allLumpedNames.size(); i++)
			if (UTL::isPresent(&addedSpecies, allLumpedNames[i]) == false)
				fullLumpedThermo << thermOut.NASAOutputLumped(allLumpedNames[i],
					{ allLumpedMols[i] }, { 300, 5000 }, { {1},{1} });

		fullLumpedThermo << "END" << std::endl;
		fullLumpedThermo << "" << std::endl;
		fullLumpedThermo.close();
	}

	currentTime = std::chrono::steady_clock::now();
	timeFile << "Mechanism merging ended at " <<
		std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count() <<
		" s" << std::endl;
	timeFile.close();
	return 0;
}






//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION DEFINITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void generateR(Molecola HC, std::vector<Molecola>* Rvec, std::vector<int>* cToRMap)
{
	std::vector<Molecola> tempRvec;
	std::vector<int> tempMap(HC.size() + 1);
	for (int i = 1; i <= HC.size(); i++)
	{
		// Generate all the radicals for the 1 groups (carbon)       
		// exept the quaternary carbon                               
		if (HC.tipo(i) != 1 && HC.tipo(i) != 10) continue; // JIAXIN
		if (HC.tipoC(i) == Cq) continue;
		Molecola newrad = HC;
		if (HC.tipo(i) == 1) newrad.tipo(i, 2);
		if (HC.tipo(i) == 10) newrad.tipo(i, 11);

		int isomero = 0;
		for (int j = 0; j < tempRvec.size(); j++) // search if it already exist
		{
			if (newrad == tempRvec[j])
			{
				tempRvec[j].isomeri++;
				tempMap[i] = j; // the radical in the i position is equal to the one in 
				//j position 
				isomero = 1;
			}
		}
		if (isomero == 0)  // if it does not already exist...     
		{
			tempRvec.push_back(newrad);
			tempMap[i] = tempRvec.size() - 1;
		}
	}

	*Rvec = tempRvec;
	*cToRMap = tempMap;
}


void generateROO(Molecola HC, std::vector<Molecola>* ROOvec, std::vector<int>* cToROOMap)
{
	//generateR(HC, ROOvec, cToROOMap);
	std::vector<Molecola> tempRvec;
	std::vector<int> tempMap(HC.size() + 1);
	for (int i = 1; i <= HC.size(); i++)
	{
		// Generate all the radicals for the 1 groups (carbon)       
		// exept the quaternary carbon                               
		if (HC.tipo(i) != 1) continue;
		if (HC.tipoC(i) == Cq) continue;
		if (HC.isole(i)) continue; // JIAXIN add conditions to ignore the vinylic sites.
		Molecola newrad = HC;
		newrad.tipo(i, 2);
		int isomero = 0;
		for (int j = 0; j < tempRvec.size(); j++) // search if it already exist
		{
			if (newrad == tempRvec[j])
			{
				tempRvec[j].isomeri++;
				tempMap[i] = j; // the radical in the i position is equal to the one in 
				//j position 
				isomero = 1;
			}
		}
		if (isomero == 0)  // if it does not already exist...     
		{
			tempRvec.push_back(newrad);
			tempMap[i] = tempRvec.size() - 1;
		}
	}

	*ROOvec = tempRvec;
	*cToROOMap = tempMap;

	// For above, generate ROO without vinylic sites.

	for (int i = 0; i < ROOvec->size(); i++)
	{
		(*ROOvec)[i].Crad_to_ROO((*ROOvec)[i].trova(2));
	}
}

void generateQOOH(Molecola HC, std::vector<Molecola>* QOOHvec,
	std::vector<std::vector<int>>* cToQOOHMap)
{
	std::vector<Molecola> tempQOOHvec;
	std::vector<std::vector<int>> tempMap(HC.size() + 1);
	for (auto& vec : tempMap)
		vec.resize(HC.size() + 1);

	std::vector<Molecola> ROOvec;
	std::vector<int> cToROOMap;
	generateROO(HC, &ROOvec, &cToROOMap);

	for (int i = 0; i < ROOvec.size(); i++)  // i iterates trough the roo   
	{
		int isomero = 0;
		int quat = 0;
		int dist1[SIZEMAX + 1]; // distance
		int dist2[SIZEMAX + 1];
		int dist3[SIZEMAX + 1];
		int dist4[SIZEMAX + 1];
		int dist[SIZEMAX + 1]; // JIAXIN add one further position.

		for (int k = 0; k <= (SIZEMAX); k++)
			dist1[k] = dist2[k] = dist3[k] = dist[k] = 0;
		//--------------------------------find the position of the i-th COO*      
		int pos_oo = 1;
		while (ROOvec[i].tipo(pos_oo) != 3) pos_oo++;
		//-----------------find all the postions of all the possible radicals...  
		ROOvec[i].scorri(pos_oo, 1, dist1);
		ROOvec[i].scorri(pos_oo, 2, dist2);
		ROOvec[i].scorri(pos_oo, 3, dist3);
		ROOvec[i].scorri(pos_oo, 4, dist4); // JIAXIN add one further position.
		for (int kk = 1; kk <= SIZEMAX; kk++)
			dist[kk] = dist1[kk] + dist2[kk] + dist3[kk] + dist4[kk];  // JIAXIN add one further position.
		//------------------------------------------- ...and saves them in dist[]   

		for (int pos_rad = 1; pos_rad <= ROOvec[i].size(); pos_rad++) // Jiaxin finds the position of C-OO*. 
			if (dist[pos_rad] == 1)					//  iterates trough this COO*  
			{										//  all the positions of the   
				isomero = 0;						//  found radicals              
				quat = 0;
				Molecola temp = ROOvec[i];
				temp.tipo(pos_oo, 4);
				temp.tipo(pos_rad, 2);
				if (temp.tipoC(pos_rad) == Cq || HC.isole(pos_rad)) quat = 1;  // JIAXIN  // if the carbon is quaternary              
				// the next time it has to be skipped      
				for (int j = 0; j < tempQOOHvec.size(); j++)

					if (tempQOOHvec[j] == temp)
					{
						isomero = 1;
						tempQOOHvec[j].isomeri++;
						tempMap[pos_oo][pos_rad] = j;
					};   // end if e for j
				if (!isomero && !quat) {
					tempQOOHvec.push_back(temp);
					tempMap[pos_oo][pos_rad] = tempQOOHvec.size() - 1;
				
				}
			} // end if and for pos_rad

	}  // end for i 

	*QOOHvec = tempQOOHvec;
	*cToQOOHMap = tempMap;
}

void generateOOQOOH(Molecola HC, std::vector<Molecola>* OOQOOHvec,
	std::vector<std::vector<int>>* cToOOQOOHMap)
{
	generateQOOH(HC, OOQOOHvec, cToOOQOOHMap);
	for (auto& mol : (*OOQOOHvec))
		mol.Crad_to_ROO(mol.trova(2));
}

std::vector<Reaction> initiationReactions(Molecola HC, Kinox* k)
{
	std::vector<Reaction> reactions;
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

	//std::cout << "atomsInBond.size() = " << atomsInBond.size() << std::endl;
	//std::cout << "atomsInBond[1][0] = " << atomsInBond[1][0] << std::endl;
	//std::cout << "atomsInBond[1][1] = " << atomsInBond[1][1] << std::endl;
	//std::cout << "atomsInBond[2][0] = " << atomsInBond[2][0] << std::endl;
	//std::cout << "atomsInBond[2][1] = " << atomsInBond[2][1] << std::endl;
	//std::cout << "atomsInBond[3][0] = " << atomsInBond[3][0] << std::endl;
	//std::cout << "atomsInBond[3][1] = " << atomsInBond[3][1] << std::endl;

	int numberOfReactions = 0;
	for (int i = 0; i < atomsInBond.size(); i++)
	{
		bool alreadyPresent = false;	// tells if the an equivalent reaction 
		// is already present
		int isPossible = HC.spezza(atomsInBond[i][0], atomsInBond[i][1], &temp1, &temp2);
		if (isPossible == 0) continue; // JIAXIN. this line ignores the rest of single bonds; change from break to continue;

		for (int j = 1; j <= HC.size(); j++)
		{
			if ((temp1 == firstRad[j] && temp2 == secondRad[j]) || (temp2 == firstRad[j]
				&& temp1 == secondRad[j])) // check if an equivalent reaction is already saved
			{
				firstRad[j].isomeri++;	// the number of equivalent reaction is stored 
				// in the firsRad.isomeri just for convinience
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

	for (int i = 1; i <= numberOfReactions; i++)
	{
		Radicale firstRadicalType = firstRad[i].tipoR(firstRad[i].trova(2));
		Radicale secondRadicalType = secondRad[i].tipoR(secondRad[i].trova(2));

		int isAllorNot = 0;
		if (firstRad[i].isAllylic(firstRad[i].trova(2)) || secondRad[i].isAllylic(secondRad[i].trova(2))) isAllorNot = 1;
		reactionComment reacomm = k->v_initiation(firstRadicalType, firstRad[i].isomeri, isAllorNot);
		if (isAllorNot == 1)
			reactions.push_back(Reaction(HC, std::vector<Molecola>{ firstRad[i], secondRad[i] },
				new double[3] { k->A, k->n, k->E }, "Class01: Unimolecular decompositions", reacomm));
		else
		reactions.push_back(Reaction(std::vector<Molecola>{ firstRad[i], secondRad[i] }, HC,
			new double[3] { k->A, k->n, k->E }, "Class01: Unimolecular decompositions", reacomm));
	}

	return reactions;
}

std::vector<Reaction> hAbstractionReactions(Molecola HC, Molecola absR, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	if (HC.kindOfSPecies() != OLE_) // JIAXIN
	{
		UTL::error("hAbstractionReactions called with HC that is not an hydrocarbon!");
		return std::vector<Reaction> {};
	}
	std::vector<Molecola> Rs;
	std::vector<int> cToRMap;
	generateR(HC, &Rs, &cToRMap);
	std::string nameHrad = k->nameHAbsRad(absR);
	std::string label = "Class02: H abstraction";
	Molecola HabsProd = absR;
	if (HabsProd.addH() == 0)
	{
		UTL::error("hAbstractionReactions called with absR that cannot host an hydrogen!");
		return std::vector<Reaction> {};
	}
	for (auto& R : Rs)
	{
		int pos_rad = R.trova(2);
		int isAllorNot = 0;
		if(HC.isole(pos_rad) == 1)
			isAllorNot = 2;
		if (HC.isAllylic(pos_rad))
			isAllorNot = 1;

		/* JIAXIN.add more conditions. P_1, S_01, S_11. (P is default P_1)
		1. For a radical, if it's not allylic && is S type && is not ole,
			scorri(R,1,found) , if any tipo(found[1]) == 'primary', return S_01, else return S_11.
		*/
		std::string OH_label = "none";

		if (absR == OH)
		{
			if (R.tipoC(pos_rad) == Cs && !R.isAllylic(pos_rad) && !R.isole(pos_rad))
			{
				int found[SIZEMAX + 1];
				R.scorri(pos_rad, 1, found);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					if (found[j] == 1 && R.tipoC(j) == Cp)
					{
						OH_label = "_01";
						break;
					}
					else OH_label = "_11";
				}
			}
		}


		//JIAXIN
		//std::cout << "Habs index=" << isAllorNot << std::endl;
		reactionComment reacomm = k->v_h_abstraction(absR, HC.tipoC(R.trova(2)), isAllorNot, 
			HC.numAbstractableH(R.trova(2)), R.isomeri, OH_label);
		reactions.push_back(Reaction(std::vector<Molecola>{ HC, absR },
			std::vector<Molecola>{ R, HabsProd }, new double[3] { k->A, k->n, k->E },
			label, reacomm));
	}

	return reactions;
}

std::vector<Reaction> O2AdditionToRReactions(std::vector<Molecola> Rs, Kinox* k) // Class.
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& R : Rs)
	{
		//JIAXIN std::cout << "R.kindOfSPecies   = " << R.kindOfSPecies() << std::endl;
		if (R.kindOfSPecies() != oleR_) // Jiaxin Convert R_ to oleR_.
		{
			UTL::error("O2AdditionToRReactions called on a non radical species");
			continue;
		}

		// JIAXIN
		int isAllorNot = 0;
		//std::cout << "R_trova=" << R.trova(2) << std::endl;
		if (R.isole(R.trova(2)))  // read the first column of molecule/radical, returns radical(2) index.
		{
			//std::cout << "Is olefin? = 1" << std::endl;
			isAllorNot = 2;
		}
		if (R.isAllylic(R.trova(2)))
		{
			isAllorNot = 1;
		}
		// JIAXIN
		//std::cout << "All index=" << isAllorNot << std::endl;
		//std::cout << "All index=" << R.isAllylic(R.trova(2)) << std::endl;
		if (isAllorNot == 2)
			continue;
			

		Molecola ROO = R;
		ROO.Crad_to_ROO(ROO.trova(2));
		Radicale radicalType = R.tipoR(R.trova(2));
		reactionComment reacomm = k->v_o2_add_r(radicalType, isAllorNot); // JIAXIN
		reactions.push_back(Reaction(std::vector<Molecola>{ R, O2 }, ROO, new double[3]
			{ k->A, k->n, k->E }, "Class31: O2 addition to oleR", reacomm));
		//chemout.wrireaDetailed( 3, k.A, k.n, k.E, r[i], roo[i]);
	}
	return reactions;
}

std::vector<Reaction> O2EliminationFromROOReactions(std::vector<Molecola> ROOs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& ROO : ROOs)
	{
		if (ROO.kindOfSPecies() != ROO_)
		{
			UTL::error("O2EliminationFromROOReactions called on a non ROO species");
			continue;
		}
		Molecola R = ROO;
		R.COOrad_to_Crad(R.trova(3));
		Radicale radicalType = R.tipoR(R.trova(2));
		reactionComment reacomm = k->v_o2_rem_roo(radicalType);
		reactions.push_back(Reaction(ROO, std::vector<Molecola>{ R, O2 }, new double[3]
			{ k->A, k->n, k->E }, "O2 elimination from ROO", reacomm));
	}
	return reactions;
}

std::vector<Reaction> O2AdditionToQOOHReactions(std::vector<Molecola> QOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& QOOH : QOOHs)
	{
		if (QOOH.kindOfSPecies() != oleQOOH_) // JIAXIN change from QOOH_ to oleQOOH_
		{
			UTL::error("O2AdditionToQOOHReactions called on a non oleQOOH species");
			continue;
		}






		// JIAXIN
		int isAllorNot = 0;
		//std::cout << "R_trova=" << R.trova(2) << std::endl;
		if (QOOH.isole(QOOH.trova(2)))  // read the first column of molecule/radical, returns radical(2) index.
		{
			//std::cout << "Is olefin? = 1" << std::endl;
			isAllorNot = 2;
		}
		if (QOOH.isAllylic(QOOH.trova(2)))
		{
			isAllorNot = 1;
		}
		// JIAXIN
		//std::cout << "All index=" << isAllorNot << std::endl;
		//std::cout << "All index=" << QOOH.isAllylic(QOOH.trova(2)) << std::endl;
		if (isAllorNot == 2)
			continue;

		int num_C = 0;
		if (QOOH.tipoR(QOOH.trova(2)) == Rs) num_C = QOOH.numberOfC();
		

		//std::cout << "HERE num C=" << num_C << std::endl;

		Molecola OOQOOH = QOOH;
		OOQOOH.Crad_to_ROO(OOQOOH.trova(2));
		Radicale radicalType = QOOH.tipoR(QOOH.trova(2));
		reactionComment reacomm = k->v_o2_add_qooh(radicalType, isAllorNot, num_C); // JIAXIN, the C number dependency.
		reactions.push_back(Reaction(std::vector<Molecola>{ QOOH, O2 }, OOQOOH,
			new double[3] { k->A, k->n, k->E }, "Class37: O2 addition to oleQOOH", reacomm));
	}
	return reactions;
}

std::vector<Reaction> O2EliminationFromOOQOOHReactions(std::vector<Molecola> OOQOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& OOQOOH : OOQOOHs)
	{
		if (OOQOOH.kindOfSPecies() != OOQOOH_)
		{
			UTL::error("O2EliminationFromOOQOOHReactions called on a non OOQOOH species");
			continue;
		}
		Molecola QOOH = OOQOOH;
		QOOH.COOrad_to_Crad(QOOH.trova(3));
		Radicale radicalType = QOOH.tipoR(QOOH.trova(2));
		reactionComment reacomm = k->v_o2_rem_ooqooh(radicalType);
		reactions.push_back(Reaction(OOQOOH, std::vector<Molecola>{ QOOH, O2 },
			new double[3] { k->A, k->n, k->E }, "O2 elimination from OOQOOH", reacomm));
	}
	return reactions;
}

std::vector<Reaction> RIsomerizationReaction(std::vector<Molecola> Rs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& R_reac : Rs)
	{
		if (R_reac.kindOfSPecies() != oleR_) // Jiaxin R_ to oleR_.
 		{
			UTL::error("RIsomerizationReaction called on a non R species");
			continue;
		}

		if (R_reac.isole(R_reac.trova(2)) == 1) continue; // JIAXIN. skip all the viny radicals.

		// return the positions of the two C=C carbons.
		int pos_ole1 = R_reac.posOle(); // position of the first C=C carbon.
		int pos_ole2 = pos_ole1 + 1;

		if (R_reac.inchiName() == "InChI=1S/C6H11/c1-3-5-6-4-2/h3H,1-2,4-6H2")
			std::cout << "here" << std::endl;

		for (int dist = 3; dist < 6; dist++)  // look for 5,6 and 7 atom ring isomerizations
		{
			Anello ring = a5;
			switch (dist)
			{
			case 3:
				ring = a5;
				break;
			case 4:
				ring = a6;
				break;
			case 5:
				ring = a7;
				break;
			default:
				continue;
				break;
			}
			int pos_rad = R_reac.trova(2);
			int trovati[SIZEMAX + 1];
			//-----------------------------isomerization involving 5 atom ring
			R_reac.scorri(pos_rad, dist, trovati);		// flag in the vector "trovati" all the carbons that are at 3 atom distance from the radical
			for (int j = 1; j <= SIZEMAX; j++)
			{
				if (trovati[j] == 1 && R_reac.numAbstractableH(j) != 0)
				{
					Molecola R_prod = R_reac;
					R_prod.Crad_to_C(R_reac.trova(2));
					R_prod.C_to_Crad(j);
					if (R_prod == R_reac)
						continue;
					// JIAXIN
					if (R_prod.isole(j) == 1) continue;
					int isAllorNot = 0;
					if (R_prod.isAllylic(R_prod.trova(2)))
						isAllorNot = 1;
					// compare if C=C (pos_ole1,2) is between radical and H site.
					int pos_radd = pos_rad;
					int j1 = j;
					if (pos_radd > j1) std::swap(pos_radd, j1);
					if ((pos_ole1 >= pos_radd && pos_ole1 <= j1) && (pos_ole2 >= pos_radd && pos_ole2 <= j1))
						continue;

					reactionComment reacomm = k->v_isomerization_r(R_reac.tipoR(pos_rad), R_reac.tipoH(j),
						ring, R_reac.numAbstractableH(j), isAllorNot);
					reactions.push_back(Reaction(R_reac, R_prod, new double[3] { k->A, k->n, k->E },
						"Class03: Alkenyl radical isomerization", reacomm));
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> ROOIsomerizationReactions(std::vector<Molecola> ROOs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& ROO : ROOs)
	{
		//if (ROO.kindOfSPecies() != ROO_)   // JIAXIN convert ROO_ to oleROO_
		if (ROO.kindOfSPecies() != oleROO_)
		{
			UTL::error("ROOIsomerizationReactions called on a non ROO species");
			continue;
		}
		if (ROO.inchiName() == "InChI=1S/C7H13O2/c1-3-5-6-7(4-2)9-8/h4,7H,2-3,5-6H2,1H3")
			std::cout << "Is olefin? = 1" << std::endl;
		int pos_o2 = ROO.trova(3);
		int trovati[SIZEMAX + 1];
		int pos_ole1 = ROO.posOle(); // position of the first C=C carbon.
		int pos_ole2 = pos_ole1 + 1;
		// JIAXIN
		int isOOally = 0;
		//std::cout << "R_trova=" << R.trova(2) << std::endl;
		if (ROO.isAllylic(pos_o2))	isOOally = 1;
		int found1[SIZEMAX + 1];
		int found2[SIZEMAX + 1];
		int found3[SIZEMAX + 1];
		int found4[SIZEMAX + 1];
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
			case 4: // Jiaxin add one more
				ringSize = a8;
				break;
			default:
				break;
			}
			ROO.scorri(pos_o2, dist, trovati);
			
		//if (dist == 1) found1 = found2;
		//if (dist == 2) found2 = trovati;
		//if (dist == 3) found3 = trovati;
		//if (dist == 4) found4 = trovati;

			
			for (int j = 1; j <= SIZEMAX; j++)          // j iterates trough the found H
			{
				int isole = 0;
				if (ROO.isole(j)) isole = 1;
				// Jiaxin modify the judging condition.
				int posO21 = pos_o2;
				int j1 = j;
				if (posO21 > j1) std::swap(posO21, j1);
					if ((pos_ole1 >= posO21 && pos_ole1 <= j1) && (pos_ole2 >= posO21 && pos_ole2 <= j1))
						continue;

				if (trovati[j] == 1 && ROO.numAbstractableH(j) != 0 && isole != 1)
				{
					// JIAXIN
					int isHally = 0;

					//if (ROO.isole(j))  // read the first column of molecule/radical, returns radical(2) index.
					//{
					//	//std::cout << "Is olefin? = 1" << std::endl;
					//	isHally = 2;
					//	continue;
					//}

					//std::cout << "R_trova=" << R.trova(2) << std::endl;
					if (ROO.isAllylic(j))	isHally = 1;

					reactionComment reacomm = k->v_isom_roo(ROO.tipoROO(pos_o2), ROO.tipoH(j),
						ringSize, ROO.numAbstractableH(j), isHally, isOOally);
					int pos_ooh = pos_o2;
					int pos_r = j;
					Molecola QOOH = ROO;
					QOOH.COOrad_to_COOH(pos_o2);
					QOOH.C_to_Crad(j);
					reactions.push_back(Reaction(ROO, QOOH, new double[3] { k->A, k->n, k->E },
						"Class33: oleROO isomerization", reacomm));
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> QOOHIsomerizationReactions(std::vector<Molecola> QOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& QOOH : QOOHs)
	{
		if (QOOH.kindOfSPecies() != QOOH_)
		{
			UTL::error("QOOHIsomerizationReactions called on a non QOOH species.");
			continue;
		}
		Anello ring = a5;
		int distance = QOOH.dist(QOOH.trova(2), QOOH.trova(4));
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
		default:
			continue;
			break;
		}
		reactionComment reacomm = k->v_isom_qooh(QOOH.tipoROOH(QOOH.trova(4)),
			QOOH.tipoR(QOOH.trova(2)), ring);
		Molecola ROO = QOOH;
		ROO.Crad_to_C(ROO.trova(2));
		ROO.COOH_to_COOrad(ROO.trova(4));
		reactions.push_back(Reaction(QOOH, ROO, new double[3] { k->A, k->n, k->E },
			"Isomerization QOOH", reacomm));
	}
	return reactions;
}

std::vector<Reaction> OOQOOHIsomerizationReactions(std::vector<Molecola> OOQOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& OOQOOH : OOQOOHs)
	{
		if (OOQOOH.kindOfSPecies() != oleO2QOOH_) // Jiaxin change from OOQOOH_ to oleO2QOOH_.
		{
			UTL::error("oleO2QOOHIsomerization called on a non OOQOOH species.");
			continue;
		}
		int j;
		int pos_o2 = OOQOOH.trova(3);
		int pos_ole1 = OOQOOH.posOle(); // position of the first C=C carbon.
		int pos_ole2 = pos_ole1 + 1;

		int trovati[SIZEMAX + 1];
		for (int distance = 1; distance <= 4; distance++)
		{
			//std::cout << "  Distance = " << distance << std::endl;
			OOQOOH.scorri(pos_o2, distance, trovati);
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
				continue;
				break;
			}

			// JIAXIN
			int isOOally = 0;
			//std::cout << "R_trova=" << R.trova(2) << std::endl;
			if (OOQOOH.isAllylic(pos_o2))	isOOally = 1;

			for (j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
			{
				if (trovati[j] == 1 && OOQOOH.numAbstractableH(j) != 0
					&& OOQOOH.tipo(j) != 4)
				{
					// JIAXIN
					int isHally = 0;

					if (OOQOOH.isole(j))  // read the first column of molecule/radical, returns radical(2) index.
					{
						//std::cout << "Is olefin? = 1" << std::endl;
						isHally = 2;
						continue;
					}

					//std::cout << "R_trova=" << R.trova(2) << std::endl;
					if (OOQOOH.isAllylic(j))	isHally = 1;
					int pos_o2dup = pos_o2;
					int j1 = j;
					if (pos_o2dup > j1) std::swap(pos_o2dup, j1);
					if ((pos_ole1 >= pos_o2dup && pos_ole1 <= j1) && (pos_ole2 >= pos_o2dup && pos_ole2 <= j1))
						continue;

					reactionComment reacomm = k->v_isom_ooqooh(OOQOOH.tipoROO(
						OOQOOH.trova(3)), OOQOOH.tipoH(j), ring,
						OOQOOH.numAbstractableH(j), isHally, isOOally);  // JIAXIN
					Molecola POOH2 = OOQOOH;
					POOH2.OO_to_OOH(OOQOOH.trova(3));
					POOH2.C_to_Crad(j);

					reactions.push_back(Reaction(OOQOOH, POOH2,
						new double[3] { k->A, k->n, k->E }, "Class41: oleO2QOOH isomerization",
						reacomm));
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> olefinsFromROOReactions(std::vector<Molecola> ROOs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& ROO : ROOs)
	{
		if (ROO.kindOfSPecies() != oleROO_) // JIAXIN change from ROO_ to oleROO_.
		{
			UTL::error("dienes from oleROO_ called on a non oleROO species.");
			continue;
		}
		int pos_o2 = ROO.trova(3);
		int trovati[SIZEMAX + 1];

		// JIAXIN
		int isOOally = 0;
		//std::cout << "R_trova=" << R.trova(2) << std::endl;
		if (ROO.isAllylic(pos_o2))	isOOally = 1;


		ROO.scorri(pos_o2, 1, trovati);	// find the carbon at distance 1 from the carbon 
		// with the OO
		for (int j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
		{
			if (trovati[j] == 1 && ROO.numAbstractableH(j) != 0)
			{
				int pos_r = j;

				if (ROO.tipoROO(pos_o2) == Rp	// Rate rule not available. However it should
					&& ROO.tipoC(pos_r) == Cp)  // happen only for ethane
					continue;

				// JIAXIN
				int isHally = 0;
				if (ROO.isole(j))  // read the first column of molecule/radical, returns radical(2) index.
				{
					//std::cout << "Is olefin? = 1" << std::endl;
					isHally = 2;
					continue;
				}
				//std::cout << "R_trova=" << R.trova(2) << std::endl;
				if (ROO.isAllylic(j))	isHally = 1;

				Molecola OLE = ROO;
				OLE.removeOO(pos_o2);
				OLE.addole(pos_o2, j);
				Molecola HO2(5);


				reactionComment reacomm = k->v_ole_par_roo(ROO.tipoROO(pos_o2),
					ROO.tipoC(pos_r), ROO.numAbstractableH(pos_r), isOOally, isHally);
				reactions.push_back(Reaction(ROO, std::vector<Molecola>{ OLE, HO2},
					new double[3] { k->A, k->n, k->E }, "Class32: oleROO HO2 elimination", reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> POOH2IsomerizationReactions(std::vector<Molecola> POOH2s, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& POOH2 : POOH2s)
	{
		if (POOH2.kindOfSPecies() != POOH2_)
		{
			UTL::error("POOH2IsomerizationReactions called on a non P(OOH)2 species.");
			continue;
		}
		std::vector<int> posOOH = POOH2.posOOHinPOOH2();
		int posR = POOH2.trova(2);
		for (auto& pos : posOOH)
		{
			int dist = POOH2.dist(posR, pos);
			Anello ring = a5;
			switch (dist)
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
			default:
				continue;
				break;
			}
			Molecola OOQOOH = POOH2;
			OOQOOH.COOH_to_COOrad(pos);
			OOQOOH.Crad_to_C(posR);
			reactionComment reacomm = k->v_isom_pooh2(POOH2.tipoROOH(pos),
				POOH2.tipoR(posR), ring);
			reactions.push_back(Reaction(POOH2, OOQOOH, new double[3] { k->A, k->n, k->E },
				"Isomerization P(OOH)2", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> RBetaDecompositioReactions(std::vector<Molecola> Rs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& R : Rs)
	{
		if (R.kindOfSPecies() != oleR_) // JIAXIN R_ to oleR_
		{
			UTL::error("RBetaDecompositioReactions called on a non oleR species.");
			continue;
		}
		int pos_rad = R.trova(2);

		//std::cout << "JIAXIN_tipo=" << R.tipoC(pos_rad) << std::endl;
		if (R.isole(pos_rad))
			continue;

		int posalfa[SIZEMAX + 1];
		R.scorri(pos_rad, 1, posalfa);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (posalfa[j] == 1) // find the C in alpha
			{
				int posbeta[SIZEMAX + 1];
				R.scorri(j, 1, posbeta);  // find the C in beta
				for (int z = 1; z <= SIZEMAX; z++)
				{
					if (posbeta[z] == 1 && z != pos_rad)
					{
						Molecola m1, m2;
						int isPossible = R.spezza1(j, z, &m1, &m2);
						if (isPossible == 0)
							continue;

						if (!m1.isCRad()) std::swap(m1, m2);

						int isAllorNot = 0;
						if (m1.isole(m1.trova(2)) == 1)
							isAllorNot = 2;
						if (m1.isAllylic(m1.trova(2)))
							isAllorNot = 1;

						reactionComment reacomm = k->v_beta_dec_r(R.tipoR(pos_rad),
							m1.tipoR(m1.trova(2)), isAllorNot);
						if (isAllorNot == 0)
							reactions.push_back(Reaction(std::vector<Molecola>{ m1, m2 }, R, 
								new double[3] { k->A, k->n, k->E }, "Class04: Alkenyl radical beta-scissions", reacomm));
						else
						{
						reactions.push_back(Reaction(R, std::vector<Molecola>{ m1, m2 },
							new double[3] { k->A, k->n, k->E }, "Class04: Alkenyl radical beta-scissions", reacomm));
						}

					}		// end pos beta
				}
			}
		}			// end pos alfa
	}
	return reactions;
}

std::vector<Reaction> olefinsFromRPlusO2Reactions(std::vector<Molecola> Rs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& R : Rs)
	{
		if (R.kindOfSPecies() != R_)
		{
			UTL::error("olefinsFromRPlusO2Reactions called on a non R species.");
			continue;
		}
		int pos_rad = R.trova(2);
		int posalfa[SIZEMAX + 1];
		R.scorri(pos_rad, 1, posalfa);
		Molecola HO2(5);
		Molecola O2(2);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (posalfa[j] == 1 && R.numAbstractableH(j) != 0)
			{
				Molecola OLE = R;
				OLE.Crad_to_C(OLE.trova(2));
				int isPossible = OLE.addole(pos_rad, j);
				if (isPossible == 0)
					continue;
				reactionComment reacomm = k->v_ole_par_r(R.numAbstractableH(j));
				reactions.push_back(Reaction(std::vector<Molecola>{ R, O2 },
					std::vector<Molecola>{ OLE, HO2 }, new double[3] { k->A, k->n, k->E },
					"Olefins from R + O2", reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> QOOHDecompositionReaction(std::vector<Molecola> QOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	Molecola OH(4);
	for (auto& QOOH : QOOHs)
	{
		if (QOOH.kindOfSPecies() != oleQOOH_)
		{
			UTL::error("QOOHDecompositionReaction called on a non oleQOOH species.");
			continue;
		}
		int posR = QOOH.trova(2);
		int posOOH = QOOH.trova(4);
		int dist = QOOH.dist(posR, posOOH);
		if (dist == 1)
		{
			if (QOOH.tipoR(posR) == Rp
				&& QOOH.tipoROOH(posOOH) == Rp)
				continue;
			Molecola OLE = QOOH.parentFuel();
			OLE.addole(posR, posOOH);

			// For the betaoleQOOH scission, use the reverse reaction rate rule,
			// Dienes + HO2 = oleQOOH (same rule as alkanes)
			reactionComment reacomm = k->v_oleqooh_decom(QOOH.tipoROOH(posOOH),
				QOOH.tipoR(posR), 1);
			reactions.push_back(Reaction(std::vector<Molecola>{ OLE, HO2}, QOOH,
				new double[3] { k->A, k->n, k->E },"Class34: oleQOOH decomposition", reacomm));
		}
		if (dist == 2)
		{
			int pos_alfa;   // find position where to break bond
			{
				int alfa_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				QOOH.scorri(posOOH, 1, alfa_ooh);
				QOOH.scorri(posR, 1, alfa_rad);
				int j;
				for (j = 1; j <= SIZEMAX; j++)
					if (alfa_ooh[j] == 1 && alfa_rad[j] == 1)
						break;
				pos_alfa = j;
			}
			// break molecule
			Molecola m1, m2;
			int isPossible = QOOH.spezza(posOOH, pos_alfa, &m1, &m2);
			if (isPossible == 0)
				continue;
			int posRprod = m1.trova(2);
			m1.Crad_to_C(posRprod);
			m1.addcheto(posRprod);
			// same rule as alkane, skip the loop while meets C=C.
			reactionComment reacomm = k->v_oleqooh_decom(QOOH.tipoROOH(posOOH),
				QOOH.tipoR(posR), 2);
			reactions.push_back(Reaction(QOOH,
				std::vector<Molecola>{m1, m2, OH}, new double[3] { k->A, k->n, k->E },
				"Class34: oleQOOH decomposition", reacomm));
		}
		if (dist == 3)
		{
			//                 find the position where to break
			int pos_break1;
			{
				int beta_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				QOOH.scorri(posOOH, 2, beta_ooh);
				QOOH.scorri(posR, 1, alfa_rad);
				int j;
				for (j = 1; j <= SIZEMAX; j++)
					if (beta_ooh[j] == 1 && alfa_rad[j] == 1) break;
				pos_break1 = j;
			};
			int pos_break2;
			{
				int alfa_ooh[SIZEMAX + 1];
				int beta_rad[SIZEMAX + 1];
				QOOH.scorri(posOOH, 1, alfa_ooh);
				QOOH.scorri(posR, 2, beta_rad);
				int j;
				for (j = 1; j <= SIZEMAX; j++)
					if (alfa_ooh[j] == 1 && beta_rad[j] == 1) break;
				pos_break2 = j;
			};
			//           generate decomposition products
			Molecola m1, m2;
			int isPossible = QOOH.spezza(pos_break1, pos_break2, &m1, &m2);
			if (isPossible == 0)
				continue;

			reactionComment reacomm = k->v_oleqooh_decom(QOOH.tipoROOH(posOOH),
				QOOH.tipoR(posR), 3);
			reactions.push_back(Reaction(QOOH, std::vector<Molecola>{m1, m2},
				new double[3] { k->A, k->n, k->E },
				"Class34: oleQOOH decomposition", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> QOOHToCEthReactions(std::vector<Molecola> QOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& QOOH : QOOHs) // 
	{
		if (QOOH.kindOfSPecies() != oleQOOH_) // Jiaxin_ change from QOOH_ to oleQOOH_.
		{
			UTL::error("QOOHToCEthReactions called on a non oleQOOH species.");
			continue;
		}
		std::string correction = "none";
		int posR = QOOH.trova(2);
		int posOOH = QOOH.trova(4);
		int pos_ole1 = QOOH.posOle(); // position of the first C=C carbon.
		int pos_ole2 = pos_ole1 + 1;
 
		// JIAXIN
		int isOOHally = 0;
		int isR = 0;
		//std::cout << "R_trova=" << R.trova(2) << std::endl;
		if (QOOH.isAllylic(posOOH))	isOOHally = 1;
		if (QOOH.isAllylic(posR))	isR = 1;

		//std::cout << "isOOHally: " << isOOHally << std::endl;
		//std::cout << "isR: " << isR << std::endl;

		int dist = QOOH.dist(posR, posOOH);
		AnelloO ring;
		if (dist == 1)
			ring = ao3;
		else if (dist == 2)
		{
			ring = ao4;
			int alpha_ooh[SIZEMAX + 1];
			int alpha_rad[SIZEMAX + 1];
			QOOH.scorri(posOOH, 1, alpha_ooh);
			QOOH.scorri(posR, 1, alpha_rad);
			for (int i = 1; i <= SIZEMAX; i++)
			{
				if (alpha_ooh[i] == 1 && alpha_rad[i] == 1)
				{
					if (QOOH.tipoC(i) == Ct)
						correction = "T";
					if (QOOH.tipoC(i) == Cq)
						correction = "Q";
				}
			}
		}
		else if (dist == 3)
		{
			ring = ao5;
			//int alpha_ooh[SIZEMAX + 1];
			//int alpha_rad[SIZEMAX + 1];
			//int beta_ooh[SIZEMAX + 1];
			//int beta_rad[SIZEMAX + 1];
			//QOOH.scorri(posOOH, 1, alpha_ooh);
			//QOOH.scorri(posR,   1, alpha_rad);
			//QOOH.scorri(posOOH, 2, beta_ooh);
			//QOOH.scorri(posR,   2, beta_rad);
			//
			//int pos1 = -1;
			//int pos2 = -1;
			//for (int i = 1; i <= SIZEMAX; i++)
			//{
			//	if (alpha_ooh[i] == 1 && beta_rad[i] == 1)
			//		pos1 = i;
			//	if (beta_ooh[i] == 1 && alpha_rad[i] == 1)
			//		pos2 = i;
			//}
			//if (pos1 != -1 && pos2 != -1)
			//{
			//	if (QOOH.tipoC(pos1) == Cq && QOOH.tipoC(pos2) == Cs)
			//		correction = "QS";
			//	if (QOOH.tipoC(pos1) == Cs && QOOH.tipoC(pos2) == Cq)
			//		correction = "SQ";
			//}
		}
		else
			continue;
		if (QOOH.tipoROOH(posOOH) == Rp && QOOH.tipoR(posR) == Rp && ring == ao3)
			continue;


		int posOOHdup = posOOH;
		int posRdup = posR;

		if (posOOHdup > posRdup) std::swap(posOOHdup, posRdup);
		if ((pos_ole1 >= posOOHdup && pos_ole1 <= posRdup) && (pos_ole2 >= posOOHdup && pos_ole2 <= posRdup))
			continue;

		// This part is to generate the products based on QOOH;
		Molecola olecEth = QOOH.parentFuel(); // removes/cleans oles, ethers, Os on a QOOH. 
		olecEth.addetero(posR, posOOH);  // add the ether O between the two positions on the "cleaned" molec.

		reactionComment reacomm = k->v_ether_from_qooh(QOOH.tipoROOH(posOOH),
			QOOH.tipoR(posR), ring, correction, isOOHally, isR);
		//// DEBUG
		//bool toMultiply = false;
		//for (int i = 1; i <= QOOH.size(); i++)
		//{
		//	bool toCheck = false;
		//	if (i == posR || i == posOOH)
		//		toCheck = false;
		//	else if (QOOH.dist(i, posR) + QOOH.dist(i, posOOH) == dist)
		//		toCheck = true;
		//	if (toCheck)
		//		if (QOOH.tipoC(i) == Ct || QOOH.tipoC(i) == Cq)
		//			toMultiply = true;
		//}
		//if(toMultiply)
		//	(*k).A *= 10;
		//// DEBUG
		reactions.push_back(Reaction(QOOH, std::vector<Molecola>{ olecEth, OH },
			new double[3] { k->A, k->n, k->E }, "Class35: oleQOOH to cyclic ethers", reacomm));
	}
	return reactions;
}

std::vector<Reaction> POOH2ToCEthOOHReactions(std::vector<Molecola> POOH2s, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& POOH2 : POOH2s)
	{
		if (POOH2.kindOfSPecies() != olePOOH2_)
		{
			UTL::error("POOH2ToCEthOOHReactions called on a non P(OOH)2 species.");
			continue;
		}
		int pos_ole1 = POOH2.posOle(); // position of the first C=C carbon.
		int pos_ole2 = pos_ole1 + 1;
		std::vector<int> posOOHs = POOH2.posOOHinPOOH2();
		int posR = POOH2.trova(2);
		for (auto& posOOH : posOOHs)
		{
			int dist = POOH2.dist(posR, posOOH);
			AnelloO ring = ao3;
			switch (dist)
			{
			case 1:
				ring = ao3;
				break;
			case 2:
				ring = ao4;
				break;
			case 3:
				ring = ao5;
				break;
			default:
				continue;
				break;
			}
			if (POOH2.tipoROOH(posOOH) == Rp && POOH2.tipoR(posR) == Rp && ring == ao3)
				continue;

			int posOOHdup = posOOH;
			int posRdup = posR;
			if ((pos_ole1 >= posOOHdup && pos_ole1 <= posRdup) && (pos_ole2 >= posOOHdup && pos_ole2 <= posRdup))
				continue;
			Molecola cEthOOH = POOH2;
			cEthOOH.removeOOH(posOOH);
			cEthOOH.Crad_to_C(posR);
			cEthOOH.addetero(posOOH, posR);
			reactionComment reacomm = k->v_ether_from_pooh2(POOH2.tipoC(posOOH),
				POOH2.tipoR(posR), ring);
			reactions.push_back(Reaction(POOH2, std::vector<Molecola>{ cEthOOH, OH },
				new double[3] { k->A, k->n, k->E }, "Class43: P(OOH)2 to cyclic ethers-OOH", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> POOH2DecompositionReaction(std::vector<Molecola> POOH2s, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	Molecola OH(4);
	for (auto& POOH2 : POOH2s)
	{
		if (POOH2.kindOfSPecies() != olePOOH2_)
		{
			UTL::error("POOH2DecompositionReaction called on a non olePOOH2 species.");
			continue;
		}
		std::vector<int> posOOHs = POOH2.posOOHinPOOH2();
		int posR = POOH2.trova(2);
		for (auto& posOOH : posOOHs)
		{
			if (POOH2.dist(posOOH, posR) == 1)
			{
				if (POOH2.tipoROOH(posOOH) == Rp && POOH2.tipoR(posR) == Rp)
					continue;
				Molecola OLEOOH = POOH2;
				OLEOOH.removeOOH(posOOH);
				OLEOOH.Crad_to_C(posR);
				OLEOOH.addole(posOOH, posR);
				reactionComment reacomm = k->v_pooh2_decom(POOH2.tipoROOH(posOOH),
					POOH2.tipoR(posR), 1);
				reactions.push_back(Reaction(POOH2, std::vector<Molecola>{ OLEOOH, HO2 },
					new double[3] { k->A, k->n, k->E }, "Class42: oleP(OOH)2 decomposition",
					reacomm));
			}
			if (POOH2.dist(posOOH, posR) == 2)
			{
				//                  find the position where to break
				int pos_alfa;
				{
					int alfa_ooh[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					POOH2.scorri(posOOH, 1, alfa_ooh);
					POOH2.scorri(posR, 1, alfa_rad);
					int j;
					for (j = 1; j <= SIZEMAX; j++)
						if (alfa_ooh[j] == 1 && alfa_rad[j] == 1) break;
					pos_alfa = j;
				}
				//                             generate the broken molecules
				if (POOH2.tipo(pos_alfa) == 4)
					continue;
				Molecola m1, m2;
				int isPossible = POOH2.spezza(posOOH, pos_alfa, &m1, &m2);
				if (isPossible == 0) continue;

				int pos = m1.trova(2);
				m1.tipo(pos, 1);
				m1.addcheto(pos);
				reactionComment reacomm = k->v_pooh2_decom(POOH2.tipoROOH(posOOH),
					POOH2.tipoR(posR), 2);
				reactions.push_back(Reaction(POOH2, std::vector<Molecola>{ m1, m2, OH },
					new double[3] { k->A, k->n, k->E }, "Class42: oleP(OOH)2 decomposition",
					reacomm));
			}
			int pos_rad = -1;		// find the position where the radical will be located
			{
				int beta_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				POOH2.scorri(posOOH, 1, beta_ooh);
				POOH2.scorri(posR, 2, alfa_rad);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					if (beta_ooh[j] == 1 && alfa_rad[j] == 1)
					{
						pos_rad = j;
						break;
					}
				}
			}
			if (pos_rad == -1)
				continue;

			if (POOH2.tipo(pos_rad) == 4)  // if there is a COOH in this position continue
				continue;

			//                  find the position where to break
			int pos_break1 = -1;
			{
				int beta_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				POOH2.scorri(posOOH, 2, beta_ooh);
				POOH2.scorri(posR, 1, alfa_rad);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					if (beta_ooh[j] == 1 && alfa_rad[j] == 1)
					{
						pos_break1 = j;
						break;
					}
				}
			}
			if (pos_break1 == -1)
				continue;

			int pos_break2 = -1;
			{
				int alfa_ooh[SIZEMAX + 1];
				int beta_rad[SIZEMAX + 1];
				POOH2.scorri(posOOH, 1, alfa_ooh);
				POOH2.scorri(posR, 2, beta_rad);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					if (alfa_ooh[j] == 1 && beta_rad[j] == 1)
					{
						pos_break2 = j;
						break;
					}
				}
			}
			if (pos_break2 == -1)
				continue;

			if (POOH2.tipo(pos_break1) == 4)
				continue;

			//                     break the molecule
			if (pos_break1 == pos_break2)
				continue;
			if (POOH2.areBonded(pos_break1, pos_break2) == 0)
				continue;

			Molecola m1, m2;
			int isPossible = POOH2.spezza(pos_break1, pos_break2, &m1, &m2);
			if (isPossible == 0)
				continue;

			reactionComment reacomm = k->v_pooh2_decom(POOH2.tipoROOH(posOOH),
				POOH2.tipoR(posR), 3);
			Reaction reac(POOH2, std::vector<Molecola>{ m1, m2},
				new double[3] { k->A, k->n, k->E }, "Class42: oleP(OOH)2 decomposition", reacomm);
			UTL::addUnique(&reactions, reac);
		}
	}
	return reactions;
}

// JIAXIN0219.
std::vector<Reaction> oleOHConsumption(std::vector<Molecola> oleOHs1, Kinox* k, ChemkinOut* chemOut)
{
	Molecola OH(4);
	Molecola H2O(6);
	Molecola CO(10);
	Molecola HCCO(11);
	std::vector<Reaction> reactions;
	std::vector<std::string> OLEOHs;
	std::string ole_inchi;
	for (auto& oleOH : oleOHs1)
	{
		ole_inchi = oleOH.inchiName();
		if (std::find(OLEOHs.begin(), OLEOHs.end(), ole_inchi) != OLEOHs.end()) continue;
		OLEOHs.push_back(oleOH.inchiName());

		if (oleOH.kindOfSPecies() != oleOH_)
		{
			UTL::error("oleOHConsumption called on a non oleOH_ species.");
			continue;
		}
		std::vector<Molecola> oleOHRs;
		std::vector<int> cToRMap;


		if (oleOH.inchiName() == "InChI=1S/C7H14O/c1-3-5-6-7(8)4-2/h4,7-8H,2-3,5-6H2,1H3")
		{
			std::cout << "here" << std::endl;
		}

		
		generateR(oleOH, &oleOHRs, &cToRMap);
		std::string nameHrad = k->nameHAbsRad(OH);
		//std::string label = "H abstraction";
		Molecola HabsProd = OH;
		if (HabsProd.addH() == 0)
		{
			UTL::error("oleOHConsumption called with OH that cannot host an hydrogen!");
			return std::vector<Reaction> {};
		}
		int num = 0;
		for (auto& oleOHR : oleOHRs)
		{
			if (!oleOHR.isAllylic(oleOHR.trova(2)) && !oleOHR.isAllylic(oleOHR.trova(11))) continue;// && !oleOHR.isAllylic(oleOHR.trova(11))) continue;
			Molecola temp = oleOHR;
			int pos_rad;

			//if (oleOHR.isAllylic(oleOHR.trova(11))) pos_rad = oleOHR.trova(11);		
			if (oleOHR.isAllylic(oleOHR.trova(2)))	pos_rad = oleOHR.trova(2);
			else pos_rad = oleOHR.trova(11);

			int posalfa[SIZEMAX + 1];
			oleOHR.scorri(pos_rad, 1, posalfa);
			for (int j = 1; j <= SIZEMAX; j++)
			{
				if (posalfa[j] == 1 || posalfa[j] == 10) // find the C in alpha
				{
					int posbeta[SIZEMAX + 1];
					oleOHR.scorri(j, 1, posbeta);  // find the C in beta
					for (int z = 1; z <= SIZEMAX; z++)
					{
						if (posbeta[z] == 1 && z != pos_rad)
						{
							Molecola m1, m2;
							int isPossible = oleOHR.spezza1(j, z, &m1, &m2);
							if (isPossible == 0) continue;
							//{
							//	oleOHR.Allylic_structure();
							//	isPossible = oleOHR.spezza1(j, z, &m1, &m2);
							//	if (isPossible == 0) continue;
							//}
								
							std::vector<Molecola> products = { m2 };
							bool IsInClud = isSpeciesIncluded(m1, chemOut);
							if (m1.kindOfSPecies() == dienesOH_ && !IsInClud)
							{
								if (m1.inchiName() == "InChI=1S/C5H8O/c1-2-3-4-5-6/h2-6H,1H3")
									std::cout << "here" << std::endl; //DecomposeoleRO(m1);
								std::vector<Molecola> dienesOH_decom;
								int dien1 = m1.posOles();
								int dien2 = dien1 + 1;
								std::cout << dien1 << std::endl;
								Molecola m3, m4;
								int isPossible = m1.spezza1(dien1, dien2, &m3, &m4);
								if (isPossible != 0)
								{
									products.push_back(m3);
									products.push_back(m4);
								}
								else 
								{
									products.push_back(m1);
								}
							}
							else
								products.push_back(m1);

							products = fullyDecomposeRO1(products);
							products.push_back(H2O);

							reactionComment reacomm = k->v_oleOHsConsum();
							reactions.push_back(Reaction(std::vector<Molecola>{oleOH, OH}, products, new double[3] { k->A, k->n, k->E },
								"Class36: OleOHs (roho2elim) consumption", reacomm));
						}		// end pos beta
					}
				}
			}			// end pos alfa
		}




	}
	return reactions;
}

std::vector<Reaction> oleOHOOHConsumption(std::vector<Molecola> OHoleOOHs, Kinox* k, ChemkinOut* chemOut)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& OHoleOOH : OHoleOOHs)
	{
		if (OHoleOOH.kindOfSPecies() != OHoleOOH_)
		{
			UTL::error("oleOHOOHConsumption called on a non OHoleOOH species.");
			continue;
		}
		int alfa_ooh[SIZEMAX + 1];					// array in which there is 1 if the atom in that position is in alfta to OOH
		OHoleOOH.scorri(OHoleOOH.trova(4), 1, alfa_ooh);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (alfa_ooh[j] == 1)
			{
				Molecola m1;
				Molecola m2;
				int isPossible = OHoleOOH.spezza1(OHoleOOH.trova(4), j, &m1, &m2);
				if (isPossible == 0)
					continue;
				m1.Crad_to_keto(m1.trova(2));
				std::vector<Molecola> prods;
				
				
				bool IsInClud = isSpeciesIncluded(m1, chemOut);

				if (m1.kindOfSPecies() == OHoleketo_ && !IsInClud)  //do not consider reactions that lead
					//to cyclic ethers CO
				{
					if (m1.inchiName() == "InChI=1S/C5H8O2/c6-4-2-1-3-5-7/h2,7H,1,3,5H2")
						std::cout << "here" << std::endl;
					std::vector<Molecola> OHoleKeto_dec = decomposeketo(m1);
					for (auto& pr : OHoleKeto_dec)
						prods.push_back(pr);
					//prods.push_back(m1);
				}
				else
				{
					prods.push_back(m1);
				}
				prods.push_back(m2);
				prods.push_back(OH);


				prods = fullyDecomposeRO1(prods);


				reactionComment reacomm = k->v_oleooh_dec(OHoleOOH.tipoROOH(OHoleOOH.trova(4)));
				reactions.push_back(Reaction(OHoleOOH, prods,
					new double[3] { k->A, k->n, k->E }, "Class40: OHOleOOHs (o2qohoohelim) consumption", reacomm));
			}
		}
	}
	return reactions;
}
std::vector<Reaction> dienesOOHConsumption(std::vector<Molecola> dieneOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& dieneOOH : dieneOOHs)
	{
		if (dieneOOH.kindOfSPecies() != dieneOOH_)
		{
			UTL::error("dieneOOHDecompositionReactions called on a non dieneOOH species.");
			continue;
		}
		int alfa_ooh[SIZEMAX + 1];					// array in which there is 1 if the atom in that position is in alfta to OOH
		dieneOOH.scorri(dieneOOH.trova(4), 1, alfa_ooh);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (alfa_ooh[j] == 1)
			{
				Molecola m1;
				Molecola m2;
				int isPossible = dieneOOH.spezza1(dieneOOH.trova(4), j, &m1, &m2);
				if (isPossible == 0) continue;
				m1.Crad_to_keto(m1.trova(2));
				reactionComment reacomm = k->v_oleooh_dec(dieneOOH.tipoROOH(dieneOOH.trova(4)));

				std::vector<Molecola> prods;
				prods.push_back(m1);
				prods.push_back(m2);
				prods.push_back(OH);
				prods = fullyDecomposeRO1(prods);

				reactions.push_back(Reaction(dieneOOH, prods,
					new double[3] { k->A, k->n, k->E }, "Class49: dieneOOHdecom (oleO2QOOHelim) consumption", reacomm));
			}
		}
	}
	return reactions;
}
std::vector<Reaction> OHcEthOOHDecompositionReactions(std::vector<Molecola> OHcEthOOHs, Kinox* k) 
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	int num = 0;
	for (auto& cEthOOH : OHcEthOOHs)
	{
		if (cEthOOH.kindOfSPecies() != OHcEthOOH_)
		{
			UTL::error("OHcEthOOHDecom called on a non OHcEthOOH species.");
			continue;
		}
		int alpha_ooh[SIZEMAX + 1];	// array in which there is 1 if the atom in that 
		// position is in alpha to OOH
		cEthOOH.scorri(cEthOOH.trova(4), 1, alpha_ooh);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (alpha_ooh[j] == 1 || alpha_ooh[j] == 10)
			{
				Molecola m1;
				Molecola m2;

				//std::cout << "species" << cEthOOH.inchiName() << std::endl;
				//if (cEthOOH.inchiName() == "InChI=1S/C5H10O4/c1-3-5(9-7)4(6)2-8-3/h3-7H,2H2,1H3")
				//{
				int isPossible = cEthOOH.spezza(cEthOOH.trova(4), j, &m1, &m2); // Open the ring.
				if (isPossible == 0)
					continue;
				//std::cout << m1 << std::endl << m2 << std::endl << std::endl;
				reactionComment reacomm = k->v_etherooh_dec(cEthOOH.tipoROOH(cEthOOH.trova(4)));
				//std::cout << "decom_sp" << m1.inchiName() << std::endl;

				std::vector<Molecola> prods;
				//std::cout << etherOOH[i] << std::endl;
				if (m2.size() == 0)    // if the ring was opened
				{
					m1.Crad_to_keto(cEthOOH.trova(4));
					prods = decomposeLinEthRO(m1);
				}
				else if (m2.kindOfSPecies() == cEthR_)  //do not consider reactions that lead
					//to cyclic ethers CO
				{
					m1.Crad_to_keto(m1.trova(2));
					prods = decomposeCEthR(m2);
					prods.push_back(m1);
				}
				else
				{
					continue;
				}

				if (prods.size() < 2)
				{
					UTL::error("In OHcEthOOHDecom less than 2 decomposition products have been found");
					continue;
				}
				prods.push_back(OH);
				prods = fullyDecomposeRO1(prods);
				prods = fullyDecomposeRO1(prods);
			
				reactions.push_back(Reaction(cEthOOH, prods, new double[3] { k->A, k->n, k->E },
					"Class47: OHcEthOOHs (fromPOHOOH2) consumption", reacomm));
				//reactions.push_back(Reaction(cEthOOH, std::vector<Molecola>{ m1, m2, OH }, new double[3] { k->A, k->n, k->E },
				//	"Ether-OOH decomposition", reacomm));
				//}
			}
		}
	}
	return reactions;
}
std::vector<Reaction> OHcEthDecompositionReactions(std::vector<Molecola> OHcEths, Kinox* k) 
{
	std::vector<Reaction> reactions;
	Molecola H2O(6);
	Molecola OH(4);
	for (auto& cEth : OHcEths)
	{
		if (cEth.kindOfSPecies() != OHcEth_)
		{
			UTL::error("OHcEthDecompositionReactions called on a non OHcEth species.");
			continue;
		}
		std::vector<int> etherPos = cEth.posEthero();
		int dist = cEth.dist(etherPos[0], etherPos[1]);
		for (int j = 0; j < 2; j++)			// iterates through the two carbon 
			// the oxygen is bonded with
		{
			int posRad = etherPos[j];
			int posCheto = etherPos[1 - j];
			int numH = cEth.numH(posCheto);
			int numH2 = cEth.numH(posRad);
			int posOH = cEth.posCOH();
			if (posOH == posCheto) continue; // Jiaxin

			if (numH == 0 && numH2 == 0) // if there are no abstractable 
				// hydrogens in alpha to the oxygen use the secondary
			{							 // path
				// look for the hydrogens to abstract
				// if there are methyls abstract it from them
				Molecola interm = cEth;
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
						posAbs = l;			// no matter if it is overwritten 
						// we need only one position
					}
				}
				if (posAbs == 0)	// if no methyl was found look for a secondary 
				{					// carbon in alpha
					carbType = Cs;
					for (int l = 1; l <= SIZEMAX; l++)
					{
						if (alfa_rad[l] == 1 && !interm.isGroup(l) && interm.tipoC(l) == Cs)
						{
							numAbsH += 2;
							posAbs = l;			// no matter if it is overwritten 
							// we need only one position
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
				interm.spezza1(posCheto, posAlfaCheto, &m1, &m2);
				int posRadM1 = m1.posCrad();
				m1.Crad_to_C(posRadM1);
				m1.addcheto(posRadM1);

				std::vector<Molecola> prods;
				prods.push_back(m1);
				prods.push_back(m2);
				prods = fullyDecomposeRO1(prods);
				prods.push_back(H2O);

				reactionComment reacomm = k->v_cyc_eth_dec(); // TO DO replace 0 with allylic or vynilic
				reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
					prods, new double[3] { k->A, k->n, k->E },
					"Class48: OHcEthDec (fromQOHOOH) consumption", reacomm));
			}
			else if (numH == 0 && numH2 != 0)
			{
				continue;
			}
			else
			{
				Molecola interm = cEth;
				interm.removeCycEther(posRad, posCheto);
				interm.addcheto(posCheto);
				interm.C_to_Crad(posRad);
				bool reactionHappened = false;
				Molecola m1, m2;
				switch (dist)
				{
				case 1:				// 3 member ring
				{
					// search if there is a carbon in alfa at the cheto and in beta 
					// at the radical
					int alfa_cheto[SIZEMAX + 1];
					int beta_rad[SIZEMAX + 1];
					interm.scorri(posCheto, 1, alfa_cheto);
					interm.scorri(posRad, 2, beta_rad);
					bool normalPathTaken = false;
					for (int l = 1; l <= SIZEMAX; l++)	// normal path: 
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
										Molecola d1, d2;
										interm.spezza(h, l, &d1, &d2);
										m1 = d1;
										m2 = d2;
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
				// JIAXIN added case 4:
				case 4:		// JIAXIN		// 6 member ring 
				{
					int alfa_cheto[SIZEMAX + 1];
					int beta_cheto[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					int beta_rad[SIZEMAX + 1];
					int gamma_cheto[SIZEMAX + 1];
					interm.scorri(posCheto, 1, alfa_cheto);
					interm.scorri(posCheto, 2, beta_cheto);
					interm.scorri(posCheto, 3, gamma_cheto);
					interm.scorri(posRad, 1, alfa_rad);
					interm.scorri(posRad, 2, beta_rad);
					for (int l = 1; l <= SIZEMAX; l++)
					{
						if (beta_cheto[l] == 1 && beta_rad[l] == 1)
						{
							for (int h = 1; h <= SIZEMAX; h++)
							{
								if (alfa_rad[h] == 1 && gamma_cheto[h] == 1)
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
					std::cout << "Molecule " << cEth
						<< " has the carbon bonding with the oxygen at a distance of "
						<< dist << ", this is not a valid cyclic ether!" << std::endl;
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
						std::cerr << "ERROR: in decomposition of cyclic ethers,"
							<< " two RO are formed, this should not happen!" << std::endl;
						std::cerr << "Press any key to continue..." << std::endl;
						std::string input_keyboard;
						std::cin >> input_keyboard;
					}

					reactionComment reacomm = k->v_cyc_eth_dec();

					if (ism1RO == true && RO_dec_prod.size() == 2)
					{
						std::vector<Molecola> prods1;
						//prods1.push_back(RO_dec_prod[0]);
						//prods1.push_back(RO_dec_prod[1]);
						prods1.push_back(m1);
						prods1.push_back(m2);
						prods1 = fullyDecomposeRO1(prods1);
						prods1.push_back(H2O);

						reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
							prods1, new double[3] { k->A, k->n, k->E },
							"Class48: OHcEthDec (fromQOHOOH) consumption", reacomm));
					}
					else if (ism2RO == true && RO_dec_prod.size() == 2)
					{
						std::vector<Molecola> prods2;
						prods2.push_back(m1);
						prods2.push_back(m2);
						//prods2.push_back(RO_dec_prod[0]);
						//prods2.push_back(RO_dec_prod[1]);
						prods2 = fullyDecomposeRO1(prods2);
						prods2.push_back(H2O);
						reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
							prods2, new double[3] { k->A, k->n, k->E },
							"Class48: OHcEthDec (fromQOHOOH) consumption", reacomm));
					}
					else
					{
						std::vector<Molecola> prods3;
						prods3.push_back(m1);
						prods3.push_back(m2);
						prods3 = fullyDecomposeRO1(prods3);
						prods3.push_back(H2O);
						reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
							prods3, new double[3] { k->A, k->n, k->E },
							"Class48: OHcEthDec (fromQOHOOH) consumption", reacomm));
					}
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> cyclicEtherOOHDecompositionReactions(std::vector<Molecola> olecEthOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	int num = 0;
	for (auto& cEthOOH : olecEthOOHs)
	{
		//num++;
		//std::cout << num << "num" << std::endl;
		if (cEthOOH.kindOfSPecies() != olecEthOOH_)
		{
			UTL::error("cyclicEtherOOHDecompositionReactions called on a non olecEthOOH species.");
			continue;
		}
		int alpha_ooh[SIZEMAX + 1];	// array in which there is 1 if the atom in that 
		// position is in alpha to OOH
		cEthOOH.scorri(cEthOOH.trova(4), 1, alpha_ooh);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (alpha_ooh[j] == 1)
			{
				Molecola m1;
				Molecola m2;
				int isPossible = cEthOOH.spezza1(cEthOOH.trova(4), j, &m1, &m2);
				if (isPossible == 0)
					continue;
				//std::cout << m1 << std::endl << m2 << std::endl << std::endl;
				reactionComment reacomm = k->v_etherooh_dec(cEthOOH.tipoROOH(cEthOOH.trova(4)));

				if (m2.size() != 0)
				{
					if (m2.inchiName() == "InChI=1S/C5H7O/c1-2-3-4-5-6/h2,4-5H,1,3H2")
						std::cout << "here" << std::endl;
				}
				
				//
				//if (m2.inchiName() == "InChI=1S/C5H7O/c1-2-3-4-5-6/h2,4-5H,1,3H2")
				//	std::cout << "here" << std::endl;

				std::vector<Molecola> prods;
				//std::cout << etherOOH[i] << std::endl;
				if (m2.size() == 0)    // if the ring was opened
				{
					m1.Crad_to_keto(cEthOOH.trova(4));
					prods = decomposeLinEthRO(m1);
				}
				else if (m2.kindOfSPecies() == olecEthR_ || m2.kindOfSPecies() == cEthR_)  // JIAXIN
					//do not consider reactions that lead to cyclic ethers CO
				{
					m1.Crad_to_keto(m1.trova(2));
					prods = decomposeCEthR(m2);
					prods.push_back(m1);
				}
				else
				{
					continue;
				}
				prods = fullyDecomposeRO1(prods);

				if (prods.size() < 2)
				{
					UTL::error("In cyclicEtherOOHDecompositionReactions less than 2 decomposition products have been found");
					continue;
				}

				prods.push_back(OH);

				reactions.push_back(Reaction(cEthOOH, prods, new double[3] { k->A, k->n, k->E },
					"Class45: olecEthOOHs (fromoleP(OOH)2) consumption", reacomm));

				//reactions.push_back(Reaction(cEthOOH, std::vector<Molecola>{ m1, m2, OH}, new double[3] { k->A, k->n, k->E },
				//	"Ether-OOH decomposition", reacomm));
			}
		}
	}
	return reactions;
}
std::vector<Reaction> cyclicEthersDecompositionReactions(std::vector<Molecola> olecEths, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola H2O(6);
	Molecola OH(4);
	for (auto& cEth : olecEths)
	{
		if (cEth.kindOfSPecies() != olecEth_)
		{
			UTL::error("cEthDecom called on a non olecEth_ species.");
			continue;
		}
		std::vector<int> etherPos = cEth.posEthero();
		int dist = cEth.dist(etherPos[0], etherPos[1]);
		for (int j = 0; j < 2; j++)			// iterates through the two carbon 
			// the oxygen is bonded with
		{
			int posRad = etherPos[j];
			int posCheto = etherPos[1 - j];
			int numH = cEth.numH(posCheto);
			int numH2 = cEth.numH(posRad);

			if (numH == 0 && numH2 == 0) // if there are no abstractable 
				// hydrogens in alpha to the oxygen use the secondary
			{							 // path
				// look for the hydrogens to abstract
				// if there are methyls abstract it from them
				Molecola interm = cEth;
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
						posAbs = l;			// no matter if it is overwritten 
						// we need only one position
					}
				}
				if (posAbs == 0)	// if no methyl was found look for a secondary 
				{					// carbon in alpha
					carbType = Cs;
					for (int l = 1; l <= SIZEMAX; l++)
					{
						if (alfa_rad[l] == 1 && !interm.isGroup(l) && interm.tipoC(l) == Cs)
						{
							numAbsH += 2;
							posAbs = l;			// no matter if it is overwritten 
							// we need only one position
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
				interm.spezza1(posCheto, posAlfaCheto, &m1, &m2);
				int posRadM1 = m1.posCrad();
				m1.Crad_to_C(posRadM1);
				m1.addcheto(posRadM1);

				reactionComment reacomm = k->v_cyc_eth_dec(); // TO DO replace 0 with allylic or vynilic
				reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
					std::vector<Molecola>{ m1, m2, H2O }, new double[3] { k->A, k->n, k->E },
					"Class46: olecEths (fromoleQOOH) consumption", reacomm));
			}
			else if (numH == 0 && numH2 != 0)
			{
				continue;
			}
			else
			{
				Molecola interm = cEth;
				interm.removeCycEther(posRad, posCheto);
				interm.addcheto(posCheto);
				interm.C_to_Crad(posRad);
				bool reactionHappened = false;
				Molecola m1, m2;
				switch (dist)
				{
				case 1:				// 3 member ring
				{
					// search if there is a carbon in alfa at the cheto and in beta 
					// at the radical
					int alfa_cheto[SIZEMAX + 1];
					int beta_rad[SIZEMAX + 1];
					interm.scorri(posCheto, 1, alfa_cheto);
					interm.scorri(posRad, 2, beta_rad);
					bool normalPathTaken = false;
					for (int l = 1; l <= SIZEMAX; l++)	// normal path: 
					{
						if (alfa_cheto[l] == 1 && beta_rad[l] == 1)
						{
							interm.spezza1(posCheto, l, &m1, &m2);
							normalPathTaken = true;
							reactionHappened = true;
							break;
						}
					}
					if (!normalPathTaken)
					{
						for (int l = 1; l <= SIZEMAX; l++)
						{
							if (beta_rad[l] == 1 && !interm.isole(l)) // JIAXIN
							{
								for (int h = 1; h <= interm.size(); h++)
								{
									if (interm.dist(h, l) == 1 && interm.dist(h, posRad) == 1)
									{
										Molecola d1, d2;
										interm.spezza1(h, l, &d1, &d2);
										m1 = d1;
										m2 = d2;
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
					std::cout << "Molecule " << cEth
						<< " has the carbon bonding with the oxygen at a distance of "
						<< dist << ", this is not a valid cyclic ether!" << std::endl;
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
						std::cerr << "ERROR: in decomposition of cyclic ethers,"
							<< " two RO are formed, this should not happen!" << std::endl;
						std::cerr << "Press any key to continue..." << std::endl;
						std::string input_keyboard;
						std::cin >> input_keyboard;
					}

					reactionComment reacomm = k->v_cyc_eth_dec();

					if (ism1RO == true && RO_dec_prod.size() == 2)
					{
						std::vector<Molecola> prods;
						prods.push_back(RO_dec_prod[0]);
						prods.push_back(RO_dec_prod[1]);
						prods.push_back(m2);
						prods.push_back(H2O);
						prods = fullyDecomposeRO1(prods);

						reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
							prods,new double[3] { k->A, k->n, k->E },
							"Class46: olecEths (fromoleQOOH) consumption", reacomm));
					}
					else if (ism2RO == true && RO_dec_prod.size() == 2)
					{
						std::vector<Molecola> prods;
						prods.push_back(m1);
						prods.push_back(RO_dec_prod[0]);
						prods.push_back(RO_dec_prod[1]);
						prods.push_back(H2O);
						prods = fullyDecomposeRO1(prods);

						reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
							prods, new double[3] { k->A, k->n, k->E },
							"Class46: olecEths (fromoleQOOH) consumption", reacomm));
					}
					else
					{
						std::vector<Molecola> prods;
						prods.push_back(m1);
						prods.push_back(m2);
						prods.push_back(H2O);
						prods = fullyDecomposeRO1(prods);
						reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
							prods, new double[3] { k->A, k->n, k->E },
							"Class46: olecEths (fromoleQOOH) consumption", reacomm));
					}
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> OOQOOHToOLEOOHReactions(std::vector<Molecola>OOQOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	for (auto& OOQOOH : OOQOOHs)
	{
		if (OOQOOH.kindOfSPecies() != oleO2QOOH_) // Jiaxin, change it from OOQOOH_ to oleO2QOOH_
		{
			UTL::error("OOQOOHToOLEOOHReactions called on a non oleO2QOOH species.");
			continue;
		}
		int posOO = OOQOOH.trova(3);
		int trovati[SIZEMAX + 1];

		// JIAXIN
		int isOOally = 0;
		//std::cout << "R_trova=" << R.trova(2) << std::endl;
		if (OOQOOH.isAllylic(posOO))	isOOally = 1;


		OOQOOH.scorri(posOO, 1, trovati);	// find the carbon at distance 1 from the carbon with the OO
		for (int j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
		{
			if (trovati[j] == 1 && OOQOOH.numAbstractableH(j) != 0
				&& OOQOOH.trova(4) != j)
			{
				if (OOQOOH.tipoROO(posOO) == Rp && OOQOOH.tipoC(j) == Cp)
					continue;
				Molecola ole = OOQOOH;
				ole.removeOO(posOO);
				ole.addole(posOO, j);

				int isHally = 0;
				if (OOQOOH.isole(j))  // read the first column of molecule/radical, returns radical(2) index.
				{
					//std::cout << "Is olefin? = 1" << std::endl;
					isHally = 2;
					continue;
				}
				//std::cout << "R_trova=" << R.trova(2) << std::endl;
				if (OOQOOH.isAllylic(j))	isHally = 1;

				reactionComment reacomm = k->v_ole_par_ooqooh(OOQOOH.tipoROO(posOO),
					OOQOOH.tipoC(j), OOQOOH.numAbstractableH(j), isOOally, isHally);
				reactions.push_back(Reaction(OOQOOH, std::vector<Molecola>{ ole, HO2 },
					new double[3] { k->A, k->n, k->E },
					"Class44: oleO2QOOH elimination to dienes-OOH", reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> OLEOOHDecompositionReactions(std::vector<Molecola> OLEOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& OLEOOH : OLEOOHs)
	{
		if (OLEOOH.kindOfSPecies() != oleOOH_)
		{
			UTL::error("OLEOOHDecompositionReactions called on a non oleOOH species.");
			continue;
		}
		int alfa_ooh[SIZEMAX + 1];					// array in which there is 1 if the atom in that position is in alfta to OOH
		OLEOOH.scorri(OLEOOH.trova(4), 1, alfa_ooh);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (alfa_ooh[j] == 1)
			{
				Molecola m1;
				Molecola m2;
				int isPossible = OLEOOH.spezza(OLEOOH.trova(4), j, &m1, &m2);
				if (isPossible == 0)
					continue;
				m1.Crad_to_keto(m1.trova(2));
				reactionComment reacomm = k->v_oleooh_dec(OLEOOH.tipoROOH(OLEOOH.trova(4)));
				reactions.push_back(Reaction(OLEOOH, std::vector<Molecola>{ m1, m2, OH },
					new double[3] { k->A, k->n, k->E }, "Ole-OOH decomposition", reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> KHPFormationReactions(std::vector<Molecola> OOQOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& OOQOOH : OOQOOHs)
	{
		if (OOQOOH.kindOfSPecies() != oleO2QOOH_) // Jiaixn here.
		{
			UTL::error("KHPFormationReactions called on a non OOQOOH species.");
			continue;
		}
		int posOO = OOQOOH.trova(3);		// find the position of the oo group
		int posOOH = OOQOOH.trova(4);		// find the position of the ooh group
		int pos_ole1 = OOQOOH.posOle(); // position of the first C=C carbon.
		int pos_ole2 = pos_ole1 + 1;

		int posOO1 = posOO;
		int posOOH1 = posOOH;

		if (posOOH1 > posOO1) std::swap(posOO1, posOOH1);
		if ((pos_ole1 >= posOOH1 && pos_ole1 <= posOO1) && (pos_ole2 >= posOOH1 && pos_ole2 <= posOO1))
			continue;

		int dist = OOQOOH.dist(posOO, posOOH);	// find the distance between the oo 
		//and ooh group
// from the distance find the type of ring formed during the reaction
		Anello ring = a5;
		switch (dist)
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
		default:		// if the ring is bigger than 7 atoms skip 
			continue;	// since the reaction cannot not happen
			break;
		}

		Molecola KHP = OOQOOH;
		KHP.COOrad_to_COOH(posOO);		// change the oo group in an ooh group
		KHP.tipo(posOOH, 1);					// remove the ooh group and ...
		int isPossible = KHP.addcheto(posOOH);	// ... replace it with the group =o 
		// Check if this step is possible.
		int roo_type = 0;
		int rooh_type = 0;
		if (OOQOOH.isAllylic(posOO))
			roo_type = 1;
		if (OOQOOH.isAllylic(posOOH))
			rooh_type = 1;


		if (isPossible != 0)					// if it is possible save reaction
		{
			reactionComment reacomm = k->v_ooqooh_to_khp(OOQOOH.tipoROO(posOO),
				OOQOOH.tipoROOH(posOOH), ring, OOQOOH.numAbstractableH(posOOH), roo_type, rooh_type);
			reactions.push_back(Reaction(OOQOOH, std::vector<Molecola>{ KHP, OH },
				new double[3] { k->A, k->n, k->E },
				"Class38: oleKHP formation", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> KHPDecompositionReactions(std::vector<Molecola> KHPs, Kinox* k, ChemkinOut* chemOut)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& KHP : KHPs)
	{
		if (KHP.kindOfSPecies() != oleKHP_)
		{
			UTL::error("KHPDecompositionReactions called on a non oleKHP species.");
			continue;
		}
		int posOOH = KHP.trova(4);		// find the position of the ooh group
		int posO = KHP.trova(7);		// find the position of the =o group
		int dist = KHP.dist(posO, posOOH);	// find the distance between the ooh and =o group
		// from the distance find the type of ring formed during the reaction
		Anello ring = a5;
		switch (dist)
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
		default:		// if the ring is bigger than 7 atoms skip
			continue;
			break;
		}

		int pos_alpha;
		{
			int alpha_ooh[SIZEMAX + 1];
			int alpha_rad[SIZEMAX + 1];
			KHP.scorri(posO, dist-1, alpha_ooh);
			KHP.scorri(posOOH, 1, alpha_rad);
			for (int j = 1; j <= SIZEMAX; j++)
			{
				if (alpha_ooh[j] == 1 && alpha_rad[j] == 1)
				{
					pos_alpha = j;
					break;
				}
			}
		}

		Molecola m1, m2;

		int isPossible = KHP.spezza1(pos_alpha, posOOH, &m1, &m2);
		if (isPossible == 0)
			continue;
		int pos_rad_m2 = m2.trova(2);	// find the radical on the molecule m2
		m2.tipo(pos_rad_m2, 1);			// remove the radical and..
		m2.addcheto(pos_rad_m2);		// ... replace it with a cheto group

		std::vector<Molecola> products = { m2 };

		bool IsInClud = isSpeciesIncluded(m1, chemOut);

		if (m1.kindOfSPecies() == RO_)// || (m1.kindOfSPecies() == oleRO_ && m1.size() > KHP.size()-2))
		{
			std::vector<Molecola> RO_dec_prod = fullyDecomposeRO(m1);
			for (auto& pr : RO_dec_prod)
				products.push_back(pr);
		}

		//else if (m1.kindOfSPecies() == oleRO_  && !IsInClud)
		//{
		//	if (m1.inchiName() == "InChI=1S/C6H9O/c1-3-5-6(7)4-2/h3,5H,2,4H2,1H3")
		//		std::cout << "here" << std::endl;
		//	std::vector<Molecola> oleRO_dec_prod = fullyDecomposeRO(m1); //DecomposeoleRO(m1);
		//	for (auto& pr1 : oleRO_dec_prod)
		//		products.push_back(pr1);
		//}
		else
			products.push_back(m1);

		products = fullyDecomposeRO1(products);

		products.push_back(OH);

		reactionComment reacomm = k->v_khp_decomp(KHP.tipoROOH(posOOH), dist);

		reactions.push_back(Reaction(KHP, products, new double[3] { k->A, k->n, k->E },
			"Class39: oleKHP decomposition", reacomm));
	}
	return reactions;
}

std::vector<Reaction> O2AdditionVinyRs(std::vector<Molecola> Rs, Kinox* k) 
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& R : Rs)
	{
		if (R.kindOfSPecies() != oleR_)
		{
			UTL::error("O2AdditionVinyRs called on a non oleR species.");
			continue;
		}
		if (!R.isole(R.trova(2))) continue;
		int pos_ole1 = R.posOle(); // position of the first C=C carbon.
		int pos_ole2 = pos_ole1 + 1;
		Molecola temp = R.parentFuel_complete();

		Molecola* firstRad;			// array of molecules where the first product is stored
		Molecola* secondRad;		// array of molecules where the second product is stored
		Molecola temp1, temp2;		// molecules where to temporary save the products
		firstRad = new Molecola[R.size() + 1]; if (!firstRad) exit(200);
		secondRad = new Molecola[R.size() + 1]; if (!secondRad) exit(201);


		bool alreadyPresent = false;	// tells if the an equivalent reaction 
		// is already present
		int isPossible = temp.spezza(pos_ole1, pos_ole2, &temp1, &temp2);
		if (isPossible == 0) continue; // JIAXIN. this line ignores the rest of single bonds; change from break to continue;

		for (int j = 1; j <= R.size(); j++)
		{
			if ((temp1 == firstRad[j] && temp2 == secondRad[j]) || (temp2 == firstRad[j]
				&& temp1 == secondRad[j])) // check if an equivalent reaction is already saved
			{
				firstRad[j].isomeri++;	// the number of equivalent reaction is stored 
				// in the firsRad.isomeri just for convinience
				alreadyPresent = true;
				break;
			}
			int asd = 1;
		}
		if (!alreadyPresent)
		{
			firstRad[0] = temp1;
			secondRad[0] = temp2;
		}

		firstRad[0].addcheto(firstRad[0].trova(2));
		secondRad[0].addcheto(secondRad[0].trova(2));
		if (R.trova(2) == pos_ole1) firstRad[0].Crad_to_C(firstRad[0].trova(2));
		if (R.trova(2) == pos_ole2) secondRad[0].Crad_to_C(secondRad[0].trova(2));

		reactionComment reacomm = k->v_o2_add_vinyR(R.tipoR(R.trova(2)));
		reactions.push_back(Reaction(std::vector<Molecola>{ R, O2}, std::vector<Molecola>{firstRad[0], secondRad[0]},
			new double[3] { k->A, k->n, k->E }, "Class07: O2 addition to vinylic radicals", reacomm));
		}
	return reactions;
}

std::vector<Reaction> allylicRadicalsFormationReactions(std::vector<Molecola> OLEs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	Molecola H2O(6);
	for (auto& OLE : OLEs)
	{
		if (OLE.kindOfSPecies() != OLE_)
		{
			UTL::error("allylicRadicalsFormationReactions called on a non olefin species.");
			continue;
		}
		for (int j = 1; j <= OLE.numberOle(); j++)	// iterate through all the double bond
		{											// in the molecule
			std::vector<int> oleCarbons = OLE.posOle(j);
			for (int l = 0; l < 2; l++)		// iterate trhough all (two) carbons 
			{								// involved in the double bond
				int posalfa[SIZEMAX + 1];
				OLE.scorri(oleCarbons[l], 1, posalfa);
				// iterate through all the carbon in alfa at the carbon involved in 
				// the double bond
				for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)
				{
					if (OLE.tipo(h) != 1)
						continue;		// the H can be abstracted only on normal carbons
					if (OLE.isole(h))
						continue;
					int numH = OLE.numH(h);
					if (numH == 0)
						continue; // if there are not hydrogen to abstract skip the reaction

					// save the type of substitution degree of the two carbons 
					// (the carbon we abstract the H from and the carbon invovled in the 
					// double bond in beta position to it) in order to decide if the radical
					// isomerizes or not
					Carbonio alfaC = OLE.tipoC(h);
					Carbonio otherOleC;
					otherOleC = OLE.tipoC(oleCarbons[1 - l]);

					int posRad, pos1Ole, pos2Ole;

					Molecola product = OLE;
					if (carbonioToInt(alfaC) < carbonioToInt(otherOleC))
					{		// radical is going to isomerize
						product.removeOle(oleCarbons[l], oleCarbons[1 - l]);
						product.addole(h, oleCarbons[l]);
						product.C_to_Crad(oleCarbons[1 - l]);
					}
					else	// radical is not going to isomerize
					{
						product.C_to_Crad(h);
					}

					reactionComment reacomm = k->v_allylic_rad_form(alfaC, numH);
					reactions.push_back(Reaction(std::vector<Molecola>{ OLE, OH },
						std::vector<Molecola>{ product, H2O },
						new double[3] { k->A, k->n, k->E },
						"Allylic radicals formation", reacomm));
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> alkenylROFormationReactions(std::vector<Molecola> AllRs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	Molecola HO2(5);
	for (auto& AllR : AllRs)
	{
		if (AllR.kindOfSPecies() != oleR_)
		{
			UTL::error("alkenylROFormationReactions called on a non oleR species.");
			continue;
		}
		Molecola product = AllR;
		product.Crad_to_COrad(product.trova(2));
		reactionComment reacomm = k->v_alkenyl_ro_form();
		reactions.push_back(Reaction(std::vector<Molecola>{ AllR, HO2 },
			std::vector<Molecola>{ product, OH }, new double[3] { k->A, k->n, k->E },
			"Alkenyl RO formation", reacomm));
	}
	return reactions;
}

/*
std::vector<Reaction> alkenylRODecompositionReactions(std::vector<Molecola> AlkROs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& AlkRO : AlkROs)
	{
		if (AlkRO.kindOfSPecies() != alkRO_)
		{
			UTL::error("alkenylRODecompositionReactions called on a non alkenyl RO species.");
			continue;
		}
		int posRO = AlkRO.trova(9);
		int posalfa[SIZEMAX + 1];
		AlkRO.scorri(posRO, 1, posalfa);		// find the carbons in alpha
		// iterate through all the carbon in alfa at the CO*
		for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)
		{
			if (AlkRO.isole(h))
				continue;	// if it is the carbon involved in the double bond continue
			Carbonio tipoC = AlkRO.tipoC(h);
			Molecola product = AlkRO;
			Molecola m1, m2;
			int isPossible = AlkRO.spezza(posRO, h, &m1, &m2);
			if (isPossible == 0)
			{
				std::cerr << "WARNING: alkenyl RO decomposition failed in breaking the bond,"
					<< "reaction skipped!" << std::endl;
				continue;
			}

			int posRadM1 = m1.trova(2);
			//m1.Crad_to_C(posRadM1);
			//m1.removeAllKeto();
			if (posRadM1 != 0)
				m1.addcheto(posRadM1);
			reactionComment reacomm = k->v_alkenyl_ro_dec(tipoC);
			reactions.push_back(Reaction(AlkRO, std::vector<Molecola>{ m1, m2 },
				new double[3] { k->A, k->n, k->E }, "Alkenyl RO decomposition", reacomm));
		}
	}
	return reactions;
}
*/

std::vector<Reaction> alkenylRODecompositionReactions(std::vector<Molecola> AlkROs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& AlkRO : AlkROs)
	{
		if (AlkRO.kindOfSPecies() != alkRO_)
		{
			UTL::error("alkenylRODecompositionReactions called on a non alkenyl RO species.");
			continue;
		}
		int posRO = AlkRO.trova(9);
		int posalfa[SIZEMAX + 1];
		AlkRO.scorri(posRO, 1, posalfa);		// find the carbons in alpha
		// iterate through all the carbon in alfa at the CO*
		std::vector<int> alphaCs;
		for (int h = 1; h <= SIZEMAX; h++)
			if (posalfa[h] == 1)
				alphaCs.push_back(h);
		if (alphaCs.size() == 0)
		{
			//std::cout << "DEBUG: Not found alpha C decomposing AlkRO " << AlkRO << std::endl;
			continue;
		}
		std::vector<int> reactingCs;
		for (auto& ind : alphaCs)
			if (AlkRO.isole(ind) == false)
				reactingCs.push_back(ind);
		if (reactingCs.size() == 0)
			reactingCs.push_back(alphaCs[0]);
		for(auto &h : reactingCs)
		{
			// DEBUG ->
			if (AlkRO.isole(h))
			{
				//std::cout << "DEBUG: The decomposition of this AlkRO would have been skipped " << AlkRO << std::endl;
				int posalfa2[SIZEMAX + 1];
				AlkRO.scorri(h, 1, posalfa2);		// find the carbons in alpha
				std::vector<int> alphaCs2;
				for (int g = 1; g <= SIZEMAX; g++)
					if (posalfa2[g] == 1 && g != posRO)
					{
						alphaCs2.push_back(g);
					}
				for (auto& posA : alphaCs2)
				{
					int posBeta[SIZEMAX + 1];
					AlkRO.scorri(posA, 1, posBeta);		// find the carbons in alpha
					for (int g = 1; g <= SIZEMAX; g++)
					{
						if (posBeta[g] == 1 && g!=h)
						{
							Molecola temp = AlkRO;
							temp.COrad_to_CO(posRO);
							temp.removeOle(posRO, h);
							temp.C_to_Crad(h);
							//std::cout << "asd1 " << temp << std::endl;
							Molecola m1a, m2a;
							int isPossible = temp.spezza(posA, g, &m1a, &m2a);
							//std::cout << isPossible << " " << posA << " " << g << std::endl;
							if (isPossible != 0)
							{
								reactionComment reacomm2("DEGUB");
								//std::cout << "DEBUG: alternative AlkRO decomposition path worked!" << std::endl;
								reactions.push_back(Reaction(AlkRO, std::vector<Molecola>{ m1a, m2a },
									new double[3] { 7.84E+10, 0.588, 33300 }, "Alkenyl RO decomposition", reacomm2));
							}
						}
					}
				}

				//continue;	// if it is the carbon involved in the double bond continue
			}
			// <- DEBUG
			Carbonio tipoC = AlkRO.tipoC(h);
			Molecola product = AlkRO;
			Molecola m1, m2;
			int isPossible = AlkRO.spezza(posRO, h, &m1, &m2);
			if (isPossible == 0)
			{
				//std::cerr << "WARNING: alkenyl RO decomposition failed in breaking the bond with species " << AlkRO
				//	<< ", reaction skipped!" << std::endl;
				continue;
			}

			int posRadM1 = m1.trova(2);
			//m1.Crad_to_C(posRadM1);
			//m1.removeAllKeto();
			if (posRadM1 != 0)
				m1.addcheto(posRadM1);
			reactionComment reacomm = k->v_alkenyl_ro_dec(tipoC);
			reactions.push_back(Reaction(AlkRO, std::vector<Molecola>{ m1, m2 },
				new double[3] { k->A, k->n, k->E }, "Alkenyl RO decomposition", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> aldehydesDecompositionReactions(std::vector<Molecola> ALDs, Kinox* k)
{
	Molecola OH(4);
	Molecola H2O(6);
	Molecola CO(10);
	std::vector<Reaction> reactions;
	for (auto& ALD : ALDs)
	{
		if (ALD.isAldehyde() == false)
		{
			UTL::error("aldehydesDecompositionReactions called on a non aldehyde species.");
			continue;
		}
		Molecola m1, m2;
		int posalfa[SIZEMAX + 1];
		ALD.scorri(ALD.posKeto(), 1, posalfa);
		for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)	// find the carbon in alpha
		{
			ALD.spezza(ALD.posKeto(), h, &m1, &m2);
		}

		reactionComment reacomm = k->v_ald_dec();
		reactions.push_back(Reaction(std::vector<Molecola>{ ALD, OH },
			std::vector<Molecola>{ m2, CO, H2O }, new double[3] { k->A, k->n, k->E },
			"Aldehydes decomposition", reacomm));
	}
	return reactions;
}

std::vector<Reaction> aldehydeOlefinsDecompositionReactions(std::vector<Molecola> ALDOLEs, Kinox* k)
{
	Molecola OH(4);
	Molecola H2O(6);
	Molecola CO(10);
	Molecola HCCO(11);
	std::vector<Reaction> reactions;
	for (auto& ALDOLE : ALDOLEs)
	{
		if (ALDOLE.isOleAldehyde() == false)
		{
			UTL::error("aldehydeOlefinsDecompositionReactions called on a non aldehyde olefin species.");
			continue;
		}
		if(ALDOLE.numAbstractableH(ALDOLE.posKeto()) == 1)
		{
			Molecola m1, m2;
			int posalfa[SIZEMAX + 1];
			ALDOLE.scorri(ALDOLE.posKeto(), 1, posalfa);
			bool toSkip = false;
			for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)	// find the carbon in alpha
			{
				int isPossible = ALDOLE.spezza1(ALDOLE.posKeto(), h, &m1, &m2);
				if (isPossible == 0)
				{
					std::cout << "DEBUG: decomposition ALDOLE did not work for " << ALDOLE << std::endl;
					toSkip = true;
				}
			}
			if (toSkip)
				continue;

			reactionComment reacomm = k->v_ald_dec();
			reactions.push_back(Reaction(std::vector<Molecola>{ ALDOLE, OH },
				std::vector<Molecola>{ m2, CO, H2O }, new double[3] { k->A, k->n, k->E },
				"Olefin aldehydes decomposition", reacomm));
		}
		else if (ALDOLE.isole(ALDOLE.posKeto()))
		{
			int posalpha[SIZEMAX + 1];
			ALDOLE.scorri(ALDOLE.posKeto(), 1, posalpha);
			int posbeta[SIZEMAX + 1];
			ALDOLE.scorri(ALDOLE.posKeto(), 2, posbeta);
			int posgamma[SIZEMAX + 1];
			ALDOLE.scorri(ALDOLE.posKeto(), 3, posgamma);
			Molecola m1, m2;
			for (int gamma = 1; gamma <= SIZEMAX; gamma++) if (posgamma[gamma] == 1 && ALDOLE.numAbstractableH(gamma) > 0)
			{
				Molecola temp = ALDOLE;
				int alfatogamma[SIZEMAX + 1];
				int betatogamma[SIZEMAX + 1];
				temp.scorri(gamma, 1, alfatogamma);
				temp.scorri(gamma, 2, betatogamma);
				for (int alpha = 1; alpha <= SIZEMAX; alpha++) if (posalpha[alpha] == 1 && betatogamma[alpha] == 1)
				{
					for (int beta = 1; beta <= SIZEMAX; beta++) if (posbeta[beta] == 1 && alfatogamma[beta] == 1)
					{
						temp.tipo(gamma, 4);
						int isPossible = temp.spezza(alpha, beta, &m1, &m2);
						m2.removeOle(1);
						m2.removeOle(1);
						int pos1 = m2.trova(4);
						int pos2 = m2.trova(2);
						m2.tipo(pos1, 1);
						m2.tipo(pos2, 1);
						m2.addole(pos1, pos2);
						if (isPossible == 0)
							continue;
						//m1.addole(m1.posKeto(), m1.posCrad());
						std::cout << "DEBUG: ALDOLE dec worked for " << ALDOLE << std::endl;
						std::cout << ALDOLE << "  " << temp << "  " << m1 << "  " << m2 << std::endl;
						reactionComment reacomm = k->v_ald_dec();
						if (m1.numberOfC() == 2)
						{
							reactions.push_back(Reaction(std::vector<Molecola>{ ALDOLE, OH },
								std::vector<Molecola>{ HCCO, m2, H2O }, new double[3] { k->A, k->n, k->E },
								"Olefin aldehydes decomposition", reacomm));
						}
						else if (m1.numberOfC() == 3)
						{
							Molecola C2H3;
							C2H3.makeC2H5();
							C2H3.addole(1, 2);
							reactions.push_back(Reaction(std::vector<Molecola>{ ALDOLE, OH },
								std::vector<Molecola>{ C2H3, CO, m2, H2O }, new double[3] { k->A, k->n, k->E },
								"Olefin aldehydes decomposition", reacomm));
						}
					}
				}
			}


			//for (int h = 1; h <= SIZEMAX; h++) if (posalpha[h] == 1)
			//{
			//	for (int g = 1; g <= SIZEMAX; g++) if (posbeta[g] == 1)
			//	{
			//		Molecola temp = ALDOLE;
			//		for (int t = 1; t <= SIZEMAX; t++) if (posgamma[t] == 1)
			//		{
			//			temp.C_to_Crad(t);
			//			break;
			//		}
			//		int isPossible = temp.spezza(h, g, &m1, &m2);
			//		int posARad[SIZEMAX + 1];
			//		//m2.scorri(m2.posCrad(), 1, posARad);
			//		//for (int t = 1; t <= SIZEMAX; t++) if (posARad[t] == 1 && m2.tipo(t) == 1)
			//		//{
			//		//	int posCrad = m2.posCrad();
			//		//	//std::cout << "m2  " << m2 << "   " << posCrad << "  " << t << std::endl;
			//		//	m2.Crad_to_C(posCrad);
			//		//	m2.addole(posCrad, t);
			//		//}
			//		if (isPossible == 0)
			//			continue;
			//		//m1.addole(m1.posKeto(), m1.posCrad());
			//		std::cout << "DEBUG: ALDOLE dec worked for " << ALDOLE << std::endl;
			//		std::cout << ALDOLE << "  " << m1 << "  " << m2 << std::endl;
			//		reactionComment reacomm = k->v_ald_dec();
			//		reactions.push_back(Reaction(std::vector<Molecola>{ ALDOLE, OH },
			//			std::vector<Molecola>{ HCCO, m2, H2O }, new double[3] { k->A, k->n, k->E },
			//			"Olefin aldehydes decomposition", reacomm));
			//	}
			//}
		}
	}
	return reactions;
}

std::vector<Reaction> ketonesDecompositionReactions(std::vector<Molecola> KETOs,
	Kinox* k)
{
	Molecola OH(4);
	Molecola H2O(6);
	std::vector<Reaction> reactions;
	for (auto& KETO : KETOs)
	{
		if (KETO.isKetone() == false)
		{
			UTL::error("ketonesDecompositionReactions called on a non ketone species.");
			continue;
		}
		for (int i = 1; i <= KETO.size(); i++)
		{
			if (KETO.numAbstractableH(i) > 0)
			{
				int posalpha[SIZEMAX + 1];
				int posbeta[SIZEMAX + 1];
				KETO.scorri(i, 1, posalpha);
				KETO.scorri(i, 2, posbeta);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					for (int h = 1; h <= SIZEMAX; h++)
					{
						if (posalpha[j] == 1 && posbeta[h] == 1)
						{
							if (KETO.dist(j, h) == 1)
							{
								Molecola Ket = KETO;
								Ket.C_to_Crad(i);
								Molecola m1, m2;
								int isPossible = Ket.spezza(j, h, &m1, &m2);
								if (isPossible == 0)
									continue;
								std::vector<Molecola> prod = { m1, m2 };
								prod = fullyDecomposeRO(prod);
								prod.push_back(H2O);
								reactionComment reacomm = k->v_h_abstraction(OH,
									KETO.tipoC(i), 0, KETO.numAbstractableH(i),  // TO DO: replace 0 with allylic or vinylic
									1, "none");
								reactions.push_back(Reaction(std::vector<Molecola>
								{ KETO, OH }, prod,
									new double[3] { k->A, k->n, k->E },
									"Ketones decomposition", reacomm));
							}
						}
					}
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> ketonesOlefinsDecompositionReactions(std::vector<Molecola> KETOOLEs,
	Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& KETOOLE : KETOOLEs)
	{
		if (KETOOLE.isOleKetone() == false)
		{
			UTL::error("ketonesOlefinsDecompositionReactions called on a non olefin ketone species.");
			continue;
		}
		int posalpha[SIZEMAX + 1];
		int posKeto = KETOOLE.posKeto();
		KETOOLE.scorri(posKeto, 1, posalpha);
		for (int i = 1; i <= SIZEMAX; i++)
		{
			if (posalpha[i] == 1)
			{
				if (KETOOLE.isole(posKeto, i))
					continue;
				Molecola m1, m2;
				KETOOLE.spezza(posKeto, i, &m1, &m2);
				reactionComment reacomm = k->v_initiation(m1.tipoR(m1.trova(2)), 1, 1); // Jiaxin 25-02-09, added one more param.
				std::vector<Molecola> prod = { m1, m2 };
				if (m1.kindOfSPecies() == unidentified_)
					continue;
				if (m2.kindOfSPecies() == unidentified_)
					continue;
				prod = fullyDecomposeRO(prod);
				reactions.push_back(Reaction(std::vector<Molecola>
				{ KETOOLE }, prod,
					new double[3] { k->A, k->n, k->E },
					"Olefins ketones decomposition", reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Molecola> decomposeRO(Molecola m2)
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
			m2.spezza1(posBeta1, posBeta2, &m3, &m4);
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

std::vector<Molecola> decomposeRO1(Molecola m2)
{
	Molecola m3, m4;
	std::vector<Molecola> vec;
	// find first atom of the beta scission

	int pos2 = m2.posCrad();
	
	int posalpha = 0; // JIAXIN.
	int posBeta1 = 0;
	int foun[SIZEMAX + 1];
	m2.scorri(pos2, 1, foun); // JIAXIN
	for (int j = 1; j <= SIZEMAX; j++)
	{
		if (foun[j] == 1) // find the C in alpha
		{
			if (m2.isole(pos2) && m2.isole(j)) continue;
			int found[SIZEMAX + 1];
			m2.scorri(pos2, 1, found);
			for (int k = 1; k < SIZEMAX + 1; k++)
			{
				if (found[k] == 1)
				{
					//if (m2.dist(m2.posCrad(), k) > 2)
					if (m2.dist(m2.posCrad(), k) == 2)
						posBeta1 = k;
				}
			}
		}
	}
	int found[SIZEMAX + 1];
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
		m2.spezza1(posBeta1, posBeta2, &m3, &m4);
		vec.push_back(m3);
		vec.push_back(m4);
	}
	else     // desired beta scission is not possible
	{
		vec.push_back(m2);
	}
	//std::cout << m2 << std::endl << m3 << std::endl << m4 << std::endl;
	return vec;
}



std::vector<Molecola> decomposeOleRO(Molecola m2)
{
	Molecola m3, m4;
	std::vector<Molecola> vec;
	int pos  = m2.posCrad();
	//int posO = m2.posKeto();
	//int dist = m2.dist(pos, posO);
	//if (dist == 1) 
	//{
	//	int label = 0;
	//	int found[SIZEMAX + 1];
	//	int found1[SIZEMAX + 1];
	//	int found2[SIZEMAX + 1];
	//	m2.scorri(pos, 1, found); // JIAXIN
	//	m2.scorri(posO, 2, found1); // JIAXIN
	//	for (int j = 1; j <= SIZEMAX; j++)
	//	{
	//		if (found[j] == 1 && found1[j] == 1)
	//		{
	//			m2.scorri(pos, 2, found2);
	//			for (int k = 1; k <= SIZEMAX; k++)
	//			{
	//				if (found2[k] == 1 && !m2.isole(j,k))
	//				{
	//					m2.spezza1(j, k, &m3, &m4);
	//					vec.push_back(m3);
	//					vec.push_back(m4);
	//					label = 1;
	//				}
	//			}
	//		}
	//	}
	//	if (label == 0) vec.push_back(m2);
	//}
	//else if (dist == 2)
	//{
	//	int label = 0;
	//	int found[SIZEMAX + 1];
	//	m2.scorri(pos, 1, found); // JIAXIN
	//	for (int j = 1; j <= SIZEMAX; j++)
	//	{
	//		if (found[j] == 1 && m2.areBonded(posO, j))
	//		{
	//			m2.spezza1(posO, j, &m3, &m4);
	//			vec.push_back(m3);
	//			vec.push_back(m4);
	//			label = 1;
	//		}
	//	}
	//	if (label == 0) vec.push_back(m2);
	//}
	//else {
	//	vec.push_back(m2);
	//}

	int label = 0;
	int found[SIZEMAX + 1];
	int found1[SIZEMAX + 1];
	//int found2[SIZEMAX + 1];
	m2.scorri(pos, 1, found); // JIAXIN
	for (int j = 1; j <= SIZEMAX; j++)
	{
		if (found[j] == 1)
		{
			m2.scorri(j, 1, found1);
			for (int k = 1; k <= SIZEMAX; k++)
			{
				if (found1[k] == 1 && k!=pos && !m2.isole(j, k) && m2.areBonded(j,k))
				{
					m2.spezza1(j, k, &m3, &m4);
					vec.push_back(m3);
					vec.push_back(m4);
					return vec;
				}
				//if (label == 1) continue;
			}
			//if (label == 1) continue;
		}
	}
	if (label == 0)
	{
		vec.push_back(m2);
	}
		
	return vec;
}

std::vector<Molecola> decomposeketo(Molecola m2) // JIAXIN decompose CCCC=O
{
	Molecola m3, m4;
	std::vector<Molecola> vec;
	// find first atom of the beta scission
	int found[SIZEMAX + 1];
	int found1[SIZEMAX + 1];
	int pos2 = m2.posKeto();

	m2.scorri(pos2, 1, found);
	int label = 0;
	for (int k = 1; k < SIZEMAX + 1; k++)
	{
		if (found[k] == 1 && !m2.isole(pos2))
		{
			label = 1;
			m2.spezza1(pos2, k, &m3, &m4);
			vec.push_back(m3);
			vec.push_back(m4);
		}
		else if (found[k] == 1 && m2.isole(pos2))
		{
			label = 1;
			m2.scorri(k, 1, found1);
			for (int j = 1; j < SIZEMAX + 1; j++)
				if (found1[j] == 1 && j != pos2)
				{
					m2.spezza1(k, j, &m3, &m4);
					vec.push_back(m3);
					vec.push_back(m4);
				}
		}
		else continue;
	}
	if (label == 0) vec.push_back(m2);

	return vec;
}

std::vector<Molecola> DecomposeoleRO(Molecola m2) // JIAXIN decompose CCCC=O
{
	Molecola m3, m4;
	std::vector<Molecola> vec;
	// find first atom of the beta scission
//try{
	int pos1 = m2.posCrad();
	int pos2 = m2.posKeto();
	if (pos1 == pos2)
	{
		int found[SIZEMAX + 1];
		m2.scorri(pos1, 1, found); // JIAXIN
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (found[j] == 1) // find the C in alpha
			{
				m2.spezza(pos1, j, &m3, &m4);
				vec.push_back(m3);
				vec.push_back(m4);
			}
		}
	}
	else
	{
	
		int posalpha = 0; // JIAXIN.
		int posBeta1 = 0;
		int foun[SIZEMAX + 1];
		m2.scorri(pos1, 1, foun); // JIAXIN
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (foun[j] == 1) // find the C in alpha
			{
				if (m2.isole(pos1) && m2.isole(j)) continue;
				int foun1[SIZEMAX + 1];
				m2.scorri(pos1, 1, foun1);
				for (int k = 1; k < SIZEMAX + 1; k++)
				{
					if (foun1[k] == 1)
					{
						m2.spezza1(j, k, &m3, &m4);
						vec.push_back(m3);
						vec.push_back(m4);
					}
					else
					{
						vec.push_back(m2);
					}
				}
			}
			else {
				vec.push_back(m2);
			}
		}
	
	}
//}
//catch (...) {  //
//	std::cout << "try catch error: " << std::endl;
//	vec.push_back(m2);
//}
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
			if (decProds[i].kindOfSPecies() == RO_ || decProds[i].kindOfSPecies() == OHRO_ )//|| decProds[i].kindOfSPecies() == oleRO_)
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

std::vector<Molecola> fullyDecomposeRO(std::vector<Molecola> vec)
{
	std::vector<Molecola> decProds;
	for (int j = 0; j < vec.size(); j++)
	{
		Molecola m2 = vec[j];
		if (m2.kindOfSPecies() == RO_ || m2.kindOfSPecies() == oleRO_ || m2.kindOfSPecies() == OHRO_)
		{
			decProds = decomposeRO(m2);
			if (decProds.size() > 1)
			{
				for (int i = 0; i < decProds.size(); i++)
				{
					if (decProds[i].kindOfSPecies() == RO_ || decProds[i].kindOfSPecies() == OHRO_ || decProds[i].kindOfSPecies() == oleRO_)
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

std::vector<Molecola> fullyDecomposeRO1(std::vector<Molecola> vec)
{
	std::vector<Molecola> decProds;
	for (int j = 0; j < vec.size(); j++)
	{
		Molecola m2 = vec[j];

		//if (m2.inchiName() == "InChI=1S/C3H5O2/c1-3(5)2-4/h3,5H,1H3")
		//	std::cout << "here" << std::endl;
		if (m2.size() < 3) continue;

		if (m2.kindOfSPecies() == RO_)
		{
			decProds = decomposeRO1(m2);
			if (decProds.size() > 1)
			{
				for (int i = 0; i < decProds.size(); i++)
				{
					if (decProds[i].kindOfSPecies() == RO_ || decProds[i].kindOfSPecies() == OHRO_ || decProds[i].kindOfSPecies() == oleRO_)
					{
						std::vector<Molecola> secDecProds;
						secDecProds = decomposeRO1(decProds[i]);
						decProds.erase(decProds.begin() + i);
						decProds.insert(decProds.end(), secDecProds.begin(), secDecProds.end());
					}
				}
			}
			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}
		else if (m2.kindOfSPecies() == oleRO_)
		{
			decProds = decomposeOleRO(m2);
			if (decProds.size() == 1) 
			{
				if (m2.isAllylic(m2.trova(2)))
				{
				m2.Allylic_structure();
				decProds = decomposeOleRO(m2);
				}
			}
			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}
		else if (m2.kindOfSPecies() == OHRO_)
		{
			//if (m2.inchiName() == "InChI=1S/C3H5O2/c1-3(5)2-4/h3,5H,1H3")
			//	std::cout << "here" << std::endl;
			decProds = decomposeOleRO(m2);

			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}
		else if (m2.kindOfSPecies() == oleROH_)
		{
			//if (m2.inchiName() == "InChI=1S/C3H5O2/c1-3(5)2-4/h3,5H,1H3")
			//	std::cout << "here" << std::endl;
			decProds = decomposeOleRO(m2);
			if (decProds.size() == 1)
			{
				if (m2.inchiName() == "InChI=1S/C4H7O/c1-2-3-4-5/h2-5H,1H3")
					std::cout << "here" << std::endl;
				if (m2.isAllylic(m2.trova(2)))
				{
					m2.Allylic_structure();
					decProds = decomposeOleRO(m2);
				}
				if (m2.isAllylic(m2.trova(11)))
				{
					m2.Allylic_structure();
					decProds = decomposeOleRO(m2);
				}

			}
			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}
		else if (m2.kindOfSPecies() == OHoleketo_)
		{
			//if (m2.inchiName() == "InChI=1S/C3H5O2/c1-3(5)2-4/h3,5H,1H3")
			//	std::cout << "here" << std::endl;
			decProds = decomposeketo(m2);
			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}
		else if (m2.kindOfSPecies() == COHCHO_)
		{
			//if (m2.inchiName() == "InChI=1S/C3H5O2/c1-3(5)2-4/h3,5H,1H3")
			//	std::cout << "here" << std::endl;
			decProds = decomposeketo(m2);
			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}
		else if (m2.kindOfSPecies() == OHdiketo_)
		{
			//if (m2.inchiName() == "InChI=1S/C3H5O2/c1-3(5)2-4/h3,5H,1H3")
			//	std::cout << "here" << std::endl;
			decProds = decomposeketo(m2);
			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}
		else if (m2.kindOfSPecies() == dienesCO_)
		{
			//if (m2.inchiName() == "InChI=1S/C3H5O2/c1-3(5)2-4/h3,5H,1H3")
			//	std::cout << "here" << std::endl;
			decProds = decomposeketo(m2);
			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}
		else if (m2.kindOfSPecies() == oleCO_)
		{
			//if (m2.inchiName() == "InChI=1S/C3H5O2/c1-3(5)2-4/h3,5H,1H3")
			//	std::cout << "here" << std::endl;
			decProds = decomposeketo(m2);
			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}

	}
	return vec;
}

std::vector<Molecola> DecomposeDiene(Molecola m2)
{
	Molecola m3, m4;
	std::vector<Molecola> vec;
	// find first atom of the beta scission
	
	int dien1 = m2.posOles();
	int dien2 = dien1 + 1;

	std::cout << dien1 << std::endl;

	int isPossible = m2.spezza1(dien1, dien2, &m3, &m4);
	if (isPossible != 0)
	{
	vec.push_back(m3);
	vec.push_back(m4);
	}
	else {
		vec.push_back(m2);
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
		mol.spezza1(pos, mol.trova(8), &m1, &m2);
		std::vector<Molecola> prods = decomposeRO1(m1);
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
		mol.spezza1(pos, mol.trova(8), &m1, &m2);
		m2.scorri(m2.posCrad(), 1, found1);

		int pos1 = 0;
		for (int i = 0; i < SIZEMAX + 1; i++)
		{
			if (found1[i] == 1)
				pos = i;
			//m2.removeAtom(m2.posCrad());
			//m2.C_to_COrad(pos);
			//Molecola m3, m4;
			//m2.spezza1(m2.trova(9), pos1, &m3, &m4);
			////std::cout << mol << std::endl << m1 << std::endl << m3 << std::endl << m4 << std::endl;
			//vec.push_back(m1);
			//vec.push_back(m3);
			//vec.push_back(m4);

		}
		m2.removeAtom(m2.posCrad());
		m2.C_to_COrad(pos);
		int foun[SIZEMAX + 1]; // JIAXIN
		m2.scorri(m2.trova(9), 1, foun); // JIAXIN
		for (int z = 1; z <= SIZEMAX; z++) // JIAXIN
		{
			if (foun[z] == 1)
				pos1 = z;
		}
		Molecola m3, m4;
		if (m2.areBonded(m2.trova(9), m2.trova(7))) // JIAXIN
			m2.spezza1(m2.trova(9), m2.trova(7), &m3, &m4);
		else   //JIAXIN
		{
			//std::cout << m2.trova(9) << std::endl;
			//std::cout << pos1 << std::endl;
			m2.spezza1(m2.trova(9), pos1, &m3, &m4);
		}  // JIAXIN
		//m2.spezza1(m2.trova(9), pos1, &m3, &m4); // JIAXIN
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
		if (!m2.ischeto(pos1))
			m2.addcheto(pos1);
		m2.C_to_Crad(pos2);
		std::vector<Molecola> prods = decomposeRO1(m2); // JIAXIN
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
			if (!m2.ischeto(pos3))
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
			m2.spezza1(pos3, pos4, &m3, &m4);
			int posRadm3 = m3.posCrad();
			m3.Crad_to_C(posRadm3);
			if (!m3.ischeto(posRadm3))
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
				m2.spezza1(pos1, pos2, &m3, &m4);
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

// JIAXIN0210
std::vector<Reaction> OHAdditionOnOlefins(Molecola HC, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	if (HC.kindOfSPecies() != OLE_) // JIAXIN
	{
		UTL::error("OH add reactions called with HC that is not an olefin!");
		return std::vector<Reaction> {};
	}
	int pos_ole1 = HC.posOle(); // position of the first C=C carbon.
	int pos_ole2 = pos_ole1 + 1; 
	// next, change pos1 into radical, change pos2 into -C-OH.
	
	Molecola ROH1 = HC;
	Molecola ROH2 = HC;
	ROH1.ole_to_ROH(pos_ole1, pos_ole2); // first pos is the -OH site; second the radical site;
	ROH2.ole_to_ROH(pos_ole2, pos_ole1);

	reactionComment reacomm = k->v_oh_add_ole(ROH1.tipoC(pos_ole1), ROH1.tipoC(pos_ole2)); // JIAXIN
	reactions.push_back(Reaction(std::vector<Molecola>{ HC, OH }, ROH1, new double[3]
		{ k->A, k->n, k->E }, "Class15: OH addition to olefins", reacomm));
	reactionComment reacomm1 = k->v_oh_add_ole(ROH2.tipoC(pos_ole2), ROH2.tipoC(pos_ole1)); // JIAXIN
	reactions.push_back(Reaction(std::vector<Molecola>{ HC, OH }, ROH2, new double[3]
		{ k->A, k->n, k->E }, "Class15: OH addition to olefins", reacomm1));
	
	return reactions;
}
std::vector<Reaction> O2AdditionToROH(std::vector<Molecola> ROHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& R : ROHs)
	{
		//JIAXIN std::cout << "R.kindOfSPecies   = " << R.kindOfSPecies() << std::endl;
		if (R.kindOfSPecies() != ROH_) // Jiaxin Convert R_ to oleR_.
		{
			UTL::error("O2AddToROH called on a non ROH species");
			continue;
		}


		Molecola ROHOO = R;
		ROHOO.Crad_to_ROO(ROHOO.trova(2));
		Radicale radicalType = R.tipoR(R.trova(2));
		reactionComment reacomm = k->v_o2_add_roh(radicalType); // JIAXIN
		reactions.push_back(Reaction(std::vector<Molecola>{ R, O2 }, ROHOO, new double[3]
			{ k->A, k->n, k->E }, "Class17: O2 addition to ROHs", reacomm));
		//chemout.wrireaDetailed( 3, k.A, k.n, k.E, r[i], roo[i]);
	}
	return reactions;
}
std::vector<Reaction> betaROHIsomerization(std::vector<Molecola> ROHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& R_reac : ROHs)
	{
		if (R_reac.kindOfSPecies() != ROH_)
		{
			UTL::error("betaROHIsomerization called on a non ROH species");
			continue;
		}

		for (int dist = 3; dist < 6; dist++)  // look for 5,6 and 7 atom ring isomerizations
		{
			Anello ring = a5;
			switch (dist)
			{
			case 3:
				ring = a5;
				break;
			case 4:
				ring = a6;
				break;
			case 5:
				ring = a7;
				break;
			default:
				continue;
				break;
			}
			int pos_rad = R_reac.trova(2);
			int trovati[SIZEMAX + 1];
			//-----------------------------isomerization involving 5 atom ring
			R_reac.scorri(pos_rad, dist, trovati);		// flag in the vector "trovati" all the carbons that are at 3 atom distance from the radical
			for (int j = 1; j <= SIZEMAX; j++)
			{
				if (trovati[j] == 1 && R_reac.numAbstractableH(j) != 0)
				{
					Molecola R_prod = R_reac;
					R_prod.Crad_to_C(R_reac.trova(2));
					R_prod.C_to_Crad(j);
					if (R_prod == R_reac)
						continue;

					reactionComment reacomm = k->v_roh_isom_r(R_reac.tipoR(pos_rad), R_reac.tipoH(j),
						ring, R_reac.numAbstractableH(j));
					reactions.push_back(Reaction(R_reac, R_prod, new double[3] { k->A, k->n, k->E },
						"Class16: betaROH isomerization", reacomm));
				}
			}
		}
	}
	return reactions;
}
// JIAXIN
std::vector<Reaction> ROHO2IsomerizationReac(std::vector<Molecola> ROHOOs, Kinox* k)  // JIAXIN ROHO2 isomerization.
{
	{
		std::vector<Reaction> reactions;
		for (auto& ROO : ROHOOs)
		{
			if (ROO.kindOfSPecies() != ROHOO_)
			{
				UTL::error("ROHO2IsomerizationReac called on a non ROHOO species");
				continue;
			}
			int pos_o2 = ROO.trova(3);
			int pos_oh = ROO.trova(10);
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
				case 4: // Jiaxin add one more
					ringSize = a8;

					break;
				default:
					break;
				}
				ROO.scorri(pos_o2, dist, trovati);


				if (ROO.inchiName() == "InChI=1S/C7H15O3/c1-2-3-4-5-7(6-8)10-9/h7-8H,2-6H2,1H3"  && dist == 1)
				{
					std::cout << "here" << pos_oh << std::endl;
					
				}
				
				for (int j = 1; j <= SIZEMAX; j++)          // j iterates trough the found H
				{
					if (trovati[j] == 1 && ROO.numAbstractableH(j) != 0)
					{
						int isOHRing = 1;
						if (dist != 4) 
						{
							if ((pos_oh > pos_o2 && pos_oh < j) || (pos_oh > j && pos_oh < pos_o2))
							{ // check if -OH is on the ring.
								isOHRing = 2;
							}
							std::cout << "here1 " << pos_oh << std::endl;
							std::cout << "here2 " << j << std::endl;
							if (pos_oh == j)
							{
								std::cout << "here" << std::endl;
								isOHRing = 2;
							}
							
						}



						//if (pos_oh == j) continue; // JLMOD0213.

						reactionComment reacomm = k->v_isom_rohoo(ROO.tipoROO(pos_o2), ROO.tipoH(j),
							ringSize, ROO.numAbstractableH(j), isOHRing);
						int pos_ooh = pos_o2;
						int pos_r = j;
						Molecola QOOH = ROO;
						QOOH.COOrad_to_COOH(pos_o2);
						QOOH.C_to_Crad(j);
						reactions.push_back(Reaction(ROO, QOOH, new double[3] { k->A, k->n, k->E },
							"Class21: ROHOO isomerization", reacomm));
					}
				}
			}
		}
		return reactions;
	}
}
// JIAXIN H addition to olefins
std::vector<Reaction> HAdditionToOlefins(Molecola HC, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola H(8);
	if (HC.kindOfSPecies() != OLE_) // JIAXIN
	{
		UTL::error("H add reactions called with HC that is not an olefin!");
		return std::vector<Reaction> {};
	}
	int pos_ole1 = HC.posOle(); // position of the first C=C carbon.
	int pos_ole2 = pos_ole1 + 1;
	// next, change pos1 into radical, change pos2 into -C-OH.

	Molecola R1 = HC;
	Molecola R2 = HC;
	R1.ole_to_R(pos_ole1, pos_ole2); // first pos is the -OH site; second the radical site;
	R2.ole_to_R(pos_ole2, pos_ole1);

	reactionComment reacomm = k->v_h_add_ole(R1.tipoC(pos_ole1), R1.tipoC(pos_ole2)); // JIAXIN
	reactions.push_back(Reaction(std::vector<Molecola>{ HC, H }, R1, new double[3]
		{ k->A, k->n, k->E }, "Class05: H addition to olefins", reacomm));
	reactionComment reacomm1 = k->v_h_add_ole(R2.tipoC(pos_ole2), R2.tipoC(pos_ole1)); // JIAXIN
	reactions.push_back(Reaction(std::vector<Molecola>{ HC, H }, R2, new double[3]
		{ k->A, k->n, k->E }, "Class05: H addition to olefins", reacomm1));

	return reactions;
}
// JIAXIN HO2 addition to allyRs
std::vector<Reaction> HO2AddToAllyRrecom(std::vector<Molecola> Rs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	for (auto& R : Rs)
	{
		//JIAXIN std::cout << "R.kindOfSPecies   = " << R.kindOfSPecies() << std::endl;
		if (R.kindOfSPecies() != oleR_) // Jiaxin Convert R_ to oleR_.
		{
			UTL::error("HO2AddToAllyRrecom called on a non radical species");
			continue;
		}
		if (!R.isAllylic(R.trova(2))) continue;

		Molecola ROOH = R;
		ROOH.Crad_to_ROOH(R.trova(2));
		Radicale radicalType = R.tipoR(R.trova(2));
		reactionComment reacomm = k->v_ho2_add_allyR_rec(radicalType); // JIAXIN
		reactions.push_back(Reaction(std::vector<Molecola>{ R, HO2 }, ROOH, new double[3]
			{ k->A, k->n, k->E }, "Class08: HO2 addition to allyRs (RA+HO2=RAOOH)", reacomm));
		//chemout.wrireaDetailed( 3, k.A, k.n, k.E, r[i], roo[i]);
	}
	return reactions;
}
// JIAXIN HO2 addition to allyRs forming RAO+OH
std::vector<Reaction> HO2AddToAllyRdecom(std::vector<Molecola> Rs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	Molecola OH(4);

	for (auto& R : Rs)
	{
		//JIAXIN std::cout << "R.kindOfSPecies   = " << R.kindOfSPecies() << std::endl;
		if (R.kindOfSPecies() != oleR_) // Jiaxin Convert R_ to oleR_.
		{
			UTL::error("HO2AddToAllyRrecom called on a non radical species");
			continue;
		}
		if (!R.isAllylic(R.trova(2))) continue;

		Molecola RAO = R;
		RAO.Crad_to_COrad(R.trova(2));

		Radicale radicalType = R.tipoR(R.trova(2));
		reactionComment reacomm = k->v_ho2_add_allyR_dec(radicalType); // JIAXIN
		reactions.push_back(Reaction(std::vector<Molecola>{ R, HO2 }, std::vector<Molecola>{RAO, OH}, new double[3]
			{ k->A, k->n, k->E }, "Class09: HO2 addition to allyRs (RA+HO2=alkROs+OH)", reacomm));
		//chemout.wrireaDetailed( 3, k.A, k.n, k.E, r[i], roo[i]);
	}
	return reactions;
}

// JIAXIN CH3 addition to allyRs
std::vector<Reaction> CH3AddToAllyR(std::vector<Molecola> Rs, Kinox* k)
{
	std::vector<Reaction> reactions;
	//Molecola CH3;
	Molecola CH3;
	CH3.makeCH3();
	for (auto& R : Rs)
	{
		//JIAXIN std::cout << "R.kindOfSPecies   = " << R.kindOfSPecies() << std::endl;
		if (R.kindOfSPecies() != oleR_) // Jiaxin Convert R_ to oleR_.
		{
			UTL::error("CH3AddToAllyR called on a non radical species");
			continue;
		}
		if (!R.isAllylic(R.trova(2))) continue;

		Molecola largerHC = R;
		largerHC.addch3(R.trova(2));
		//Radicale radicalType = R.tipoR(R.trova(2));
		reactionComment reacomm = k->v_ch3_add_allyR(); // JIAXIN
		reactions.push_back(Reaction(std::vector<Molecola>{ R, CH3 }, largerHC, new double[3]
			{ k->A, k->n, k->E }, "Class12: CH3 addition to allyRs (RACH3)", reacomm));
		//chemout.wrireaDetailed( 3, k.A, k.n, k.E, r[i], roo[i]);
	}
	return reactions;
}
// JIAXIN HO2 addition to olefins.
std::vector<Reaction> HO2AddOlefinToQOOH(Molecola HC, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	if (HC.kindOfSPecies() != OLE_) // JIAXIN
	{
		UTL::error("HO2AddOlefinToQOOH called with HC that is not an olefin!");
		return std::vector<Reaction> {};
	}
	int pos_ole1 = HC.posOle(); // position of the first C=C carbon.
	int pos_ole2 = pos_ole1 + 1;
	// next, change pos1 into radical, change pos2 into -C-OH.

	Molecola QOOH1 = HC;
	Molecola QOOH2 = HC;
	QOOH1.ole_to_QOOH(pos_ole1, pos_ole2); // first pos is the -OH site; second the radical site;
	QOOH2.ole_to_QOOH(pos_ole2, pos_ole1);

	reactionComment reacomm = k->v_ho2_add_ole(QOOH1.tipoC(pos_ole1), QOOH1.tipoC(pos_ole2)); // JIAXIN
	reactions.push_back(Reaction(std::vector<Molecola>{ HC, HO2 }, QOOH1, new double[3]
		{ k->A, k->n, k->E }, "Class13: HO2 addition to olefins", reacomm));
	reactionComment reacomm1 = k->v_ho2_add_ole(QOOH2.tipoC(pos_ole2), QOOH2.tipoC(pos_ole1)); // JIAXIN
	reactions.push_back(Reaction(std::vector<Molecola>{ HC, HO2 }, QOOH2, new double[3]
		{ k->A, k->n, k->E }, "Class13: HO2 addition to olefins", reacomm1));

	return reactions;
}

/* Here only consider two possible channels : beta - CC break and beta - CH break.
Pathways here:
P1: addition and producing a H;
P2: If found, break the corresponding C-C bond;
*/
std::vector<Reaction> OAtomAddToOlefins(Molecola HC, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	Molecola O(3);
	Molecola H(8);
	if (HC.kindOfSPecies() != OLE_) // JIAXIN
	{
		UTL::error("OAtomAddToOlefins called with HC that is not an olefin!");
		return std::vector<Reaction> {};
	}
	int pos_ole1 = HC.posOle(); // position of the first C=C carbon.
	int pos_ole2 = pos_ole1 + 1;
	// next, change pos1 into radical, change pos2 into -C-OH.

	Molecola RO1 = HC;
	Molecola RO2 = HC;
	RO1.ole_to_QO(pos_ole1, pos_ole2); //producing the intermediate species first.
	RO2.ole_to_QO(pos_ole2, pos_ole1);
	/* Algo here:
	For R, 
	P1: addcheto(pos1) and another H;
	P2: If found, break the common beta bond of R.
	*/
	Molecola P11 = RO1; //channel1 of RO1
	Molecola P21 = RO2; //channel1 of RO2
	P11.addcheto1(pos_ole1); // replace -CO* with C=O.
	P21.addcheto1(pos_ole2); // replace -CO* with C=O.
	
	reactionComment reacomm = k->v_o_add_ole(P11.tipoC(pos_ole1), "H"); // JIAXIN
	reactions.push_back(Reaction(std::vector<Molecola>{ HC, O }, std::vector<Molecola> {P11, H}, new double[3]
		{ k->A, k->n, k->E }, "Class06: O addition to olefins", reacomm));
	
	reactionComment reacomm1 = k->v_o_add_ole(P21.tipoC(pos_ole2), "H"); // JIAXIN
	reactions.push_back(Reaction(std::vector<Molecola>{ HC, O }, std::vector<Molecola> {P21, H}, new double[3]
		{ k->A, k->n, k->E }, "Class06: O addition to olefins", reacomm1));


	Molecola P12 = RO1; //channel2 of RO1
	int posalfa[SIZEMAX + 1];
	P12.scorri(pos_ole2, 1, posalfa);
	for (int j = 1; j <= SIZEMAX; j++)
	{
		if ((posalfa[j] == 1) && (j == pos_ole1)) // find the C in alpha
		{
			int posbeta[SIZEMAX + 1];
			P12.scorri(j, 1, posbeta);  // find the C in beta
			for (int z = 1; z <= SIZEMAX; z++)
			{
				if (posbeta[z] == 1 && z != pos_ole2)
				{
					Molecola m1, m2;
					int isPossible = P12.spezza1(j, z, &m1, &m2);
					if (isPossible == 0)
						continue;
					if (m1.kindOfSPecies() == OLE_) {
						for (int k = 1; k <= SIZEMAX; k++)
						{
							if (m1.isole(k) && m1.tipoC(k) == Cp)
							{
								m1.tipo(k, 9);
							}
						}
					}
					if (m2.kindOfSPecies() == OLE_) {
						for (int k = 1; k <= SIZEMAX; k++)
						{
							if (m2.isole(k) && m2.tipoC(k) == Cp)
							{
								m2.tipo(k, 9);
							}
						}
					}

					reactionComment reacomm = k->v_o_add_ole(P12.tipoC(pos_ole1), "Other"); // JIAXIN
					reactions.push_back(Reaction(std::vector<Molecola>{ HC, O }, std::vector<Molecola>{ m1, m2 }, new double[3] 
						{ k->A, k->n, k->E },	"Class06: O addition to olefins", reacomm));
				}		// end pos beta
			}
		}
	}

	Molecola P22 = RO2; //channel2 of RO2
	int posalfa1[SIZEMAX + 1];
	P22.scorri(pos_ole1, 1, posalfa1);
	for (int j = 1; j <= SIZEMAX; j++)
	{
		if ((posalfa1[j] == 1) && (j == pos_ole2)) // find the C in alpha
		{
			int posbeta1[SIZEMAX + 1];
			P22.scorri(j, 1, posbeta1);  // find the C in beta
			for (int z = 1; z <= SIZEMAX; z++)
			{
				if (posbeta1[z] == 1 && z != pos_ole1)
				{
					Molecola m3, m4;
					int isPossible = P22.spezza1(j, z, &m3, &m4);
					if (isPossible == 0)
						continue;

					if (m3.kindOfSPecies() == R_) 
					{
						int pos = m3.trova(2);
						m3.Crad_to_C(pos);
						m3.addcheto(pos);

					}
					if (m4.kindOfSPecies() == R_) 
					{
						int pos1 = m4.trova(2);
						m4.Crad_to_C(pos1);
						m4.addcheto(pos1);
					}
					reactionComment reacomm = k->v_o_add_ole(P22.tipoC(pos_ole2), "Other"); // JIAXIN
					reactions.push_back(Reaction(std::vector<Molecola>{ HC, O }, std::vector<Molecola>{ m3, m4 }, new double[3]
						{ k->A, k->n, k->E }, "Class06: O addition to olefins", reacomm));
				}		// end pos beta
			}
		}
	}
	return reactions;
}

std::vector<Reaction> ROOHdecomposition(std::vector<Molecola> ROOHs, Kinox* k)      // JIAXIN RAOOH=RAO+OH;
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& ROOH : ROOHs)
	{
		if (ROOH.kindOfSPecies() != oleOOH_)
		{
			UTL::error("ROOHdecomposition called on a non ROOH species.");
			continue;
		}
		int posROOH = ROOH.trova(4);
		Molecola RO = ROOH;
		RO.ROOH_to_RO(posROOH);

		reactionComment reacomm = k->v_rooh_to_ro();
		reactions.push_back(Reaction(ROOH, std::vector<Molecola>{ RO, OH }, new double[3] { k->A, k->n, k->E },
			"Class10: RAOOH=RAO+OH", reacomm));
	}
	return reactions;
}

std::vector<Reaction> alkROdecomposition(std::vector<Molecola> alkROs, Kinox* k)     // JIAXIN RAO to dienes;
{
	std::vector<Reaction> reactions;

	for (auto& alkRO : alkROs)
	{
		if (alkRO.kindOfSPecies() != alkRO_)
		{
			UTL::error("alkROdecomposition called on a non ROOH species.");
			continue;
		}

		int pos_RO = alkRO.trova(9);
		int posalfa[SIZEMAX + 1];
		alkRO.scorri(pos_RO, 1, posalfa);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (posalfa[j] == 1) // find the C in alpha
			{
				Molecola m1, m2;
				int isPossible = alkRO.spezza1(pos_RO, j, &m1, &m2);
				if (isPossible == 0)
					continue;
				std::string IsViny = "No";
				if (alkRO.isole(j)) IsViny = "Yes";

				reactionComment reacomm = k->v_ro_beta(IsViny);
				reactions.push_back(Reaction(alkRO,	std::vector<Molecola>{ m1, m2 },
					new double[3] { k->A, k->n, k->E },
					"Class11: RAO decomposition", reacomm));
							// end pos beta
				
			}
		}			// end pos alfa
	}
	return reactions;
}
/* Jiaxin RACH3 oxidation. wait to be added.
std::vector<Reaction> higherOleconsumption(std::vector<Molecola> higheroles, Kinox* k)      // JIAXIN RAOOH=RAO+OH;
{}
*/
//std::vector<Reaction> O2AdditionVinyRs(std::vector<Molecola> Rs, Kinox* k)      // JIAXIN RAOOH=RAO+OH;
//{
//
//}
std::vector<Reaction> O2allyRHabstraction(std::vector<Molecola> Rs, Kinox* k)      // JIAXIN RAOOH=RAO+OH; change here
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	Molecola O2(2);
	for (auto& R : Rs)
	{
		if (R.kindOfSPecies() != oleR_)
		{
			UTL::error("O2allyRHabstraction called on a non oleR species.");
			continue;
		}
		if (!R.isAllylic(R.trova(2))) continue;
		int pos_ally = R.trova(2);
		Molecola Diene = R;
		int posalfa[SIZEMAX + 1];
		R.scorri(pos_ally, 1, posalfa);
		int pos_abs = 0;
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if ((posalfa[j] == 1) && (!R.isole(j))) // find the C in alpha
			{
				//if (!Diene.isole(j)) continue;
				Diene.Crad_to_C(pos_ally);
				Diene.addole(pos_ally, j);
				pos_abs = j;
				reactionComment reacomm = k->v_o2_habs_allyR(R.tipoC(pos_abs), R.numAbstractableH(pos_abs), R.isomeri);
				reactions.push_back(Reaction(std::vector<Molecola>{ R, O2 },
					std::vector<Molecola>{ Diene, HO2 }, new double[3] { k->A, k->n, k->E },
					"Class14: allyRs HAA by O2 reactions", reacomm));
			}
		}
	}
	return reactions;
}

// Jiaxin: Only for beta-ROHOO. see Zhou2022.
std::vector<Reaction> ROHO2ToWaddington(std::vector<Molecola> ROHOOs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	for (auto& ROHOO : ROHOOs)
	{
		if (ROHOO.kindOfSPecies() != ROHOO_)
		{
			UTL::error("ROHO2ToWaddington called on a non ROHOO species.");
			continue;
		}
		int posOH = ROHOO.trova(10);
		int posOO = ROHOO.trova(3);
		if (std::abs(posOH - posOO) != 1) continue;

		Molecola OROOH = ROHOO;
		OROOH.rohoo_to_orooh(posOH, posOO);

		reactionComment reacomm = k->v_beta_rohoo_wad(ROHOO.tipoC(posOO), ROHOO.tipoC(posOH));
		reactions.push_back(Reaction(ROHOO, OROOH, new double[3] { k->A, k->n, k->E },
			"Class18: betaROHO2 Waddington mechanism", reacomm));
	}
	return reactions;
}
// Jiaxin: Only for beta-ROHOO. see Zhou2022.
std::vector<Reaction> WadOROOHdecomposition(std::vector<Molecola> OROOHs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& OROOH : OROOHs)
	{
		if (OROOH.kindOfSPecies() != OROOHs_) //
		{
			UTL::error("WadOROOHdecomposition called on a non OROOH species.");
			continue;
		}
		int posO = OROOH.trova(9);
		int posOOH = OROOH.trova(4);
		if (std::abs(posO - posOOH) != 1) continue;
		//std::cout << "JIAXIN_tipo=" << R.tipoC(pos_rad) << std::endl;
		Molecola ORO = OROOH;
		ORO.ROOH_to_RO(posOOH);

		Molecola m1, m2;
		int isPossible = ORO.spezza(posO, posOOH, &m1, &m2);
		if (isPossible == 0)
			continue;

		reactionComment reacomm = k->v_wad_decom(OROOH.tipoC(posOOH),
			OROOH.tipoC(posO));
		reactions.push_back(Reaction(OROOH, std::vector<Molecola>{ m1, m2, OH },
			new double[3] { k->A, k->n, k->E },
			"Class19: Waddington product OROOH decompositions", reacomm));
			// end pos beta
	}
	return reactions;
}
std::vector<Reaction> betaROHO2Elimination(std::vector<Molecola> ROHOOs, Kinox* k)     // JIAXIN
{
	std::vector<Reaction> reactions;
	for (auto& ROHOO : ROHOOs)
	{
		if (ROHOO.kindOfSPecies() != ROHOO_) // 
		{
			UTL::error("betaROHO2Elimination called on a non ROHOO species.");
			continue;
		}
		int pos_o2 = ROHOO.trova(3);
		int pos_oh = ROHOO.trova(10);
		int trovati[SIZEMAX + 1];

		ROHOO.scorri(pos_o2, 1, trovati);	// find the carbon at distance 1 from the carbon 
		// with the OO
		for (int j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
		{
			if (trovati[j] == 1 && ROHOO.numAbstractableH(j) != 0)
			{
				int pos_r = j;

				if (ROHOO.tipoROO(pos_o2) == Rp	// Rate rule not available. However it should
					&& ROHOO.tipoC(pos_r) == Cp)  // happen only for ethane
					continue;

				Molecola oleOH = ROHOO;
				oleOH.removeOO(pos_o2);
				oleOH.addole(pos_o2, j);
				Molecola HO2(5);

				std::string oh_to_r = "No"; // OH and expelled H are not on the same C.
				if (pos_r == pos_oh) oh_to_r = "Yes";



				std::cout << ROHOO.numAbstractableH(pos_r) << std::endl;
				reactionComment reacomm = k->v_beta_rohoo_elim(ROHOO.tipoC(pos_o2),
					ROHOO.tipoC(pos_r), ROHOO.numAbstractableH(pos_r), oh_to_r);
				reactions.push_back(Reaction(ROHOO, std::vector<Molecola>{ oleOH, HO2},
					new double[3] { k->A, k->n, k->E }, "Class20: ROHO2 elimination reactions added", reacomm));
			}
		}
	}
	return reactions;
}
std::vector<Reaction> betaOHQOOHScission(std::vector<Molecola> QOHOOHs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	for (auto& QOHOOH : QOHOOHs)
	{
		if (QOHOOH.kindOfSPecies() != QOHOOH_) //
		{
			UTL::error("betaOHQOOHScission called on a non QOHOOH species.");
			continue;
		}
		int pos_o2h = QOHOOH.trova(4);
		int pos_oh  = QOHOOH.trova(10);
		int pos_rad = QOHOOH.trova(2);
		if (QOHOOH.trova(11))
		{
			pos_rad = QOHOOH.trova(11);
			pos_oh = QOHOOH.trova(11);

		}

		if (std::abs(pos_rad - pos_o2h) != 1) continue;

		Molecola oleOH = QOHOOH;
		oleOH.removeOOH(pos_o2h);
		//if (pos_rad != pos_oh) 
		oleOH.Crad_to_C(pos_rad); //JLMOD0213
		oleOH.addole(pos_o2h, pos_rad);
		Molecola HO2(5);

		std::string oh_to_r = "No"; // OH and expelled H are not on the same C.
		if (pos_rad == pos_oh) oh_to_r = "Yes";

		reactionComment reacomm = k->v_beta_qohooh_decom(QOHOOH.tipoROOH(pos_o2h),
			QOHOOH.tipoR(pos_rad), oh_to_r);
		std::string RR = reacomm.rateRule;
		if (!RR.empty() && RR.substr(RR.find_first_not_of(" \t")).front() == '-')
		{
			reactions.push_back(Reaction(std::vector<Molecola>{ oleOH, HO2}, QOHOOH,
				new double[3] { k->A, k->n, k->E }, "Class22: beta-QOHOOH scissions", reacomm));
		}
		else
		reactions.push_back(Reaction(QOHOOH, std::vector<Molecola>{ oleOH, HO2},
			new double[3] { k->A, k->n, k->E }, "Class22: beta-QOHOOH scissions", reacomm));
	}

	return reactions;
}
std::vector<Reaction> OHQOOHToCycEther(std::vector<Molecola> QOHOOHs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& QOHOOH : QOHOOHs)
	{
		if (QOHOOH.kindOfSPecies() != QOHOOH_) // Jiaxin_ change from QOOH_ to oleQOOH_.
		{
			UTL::error("OHQOOHToCycEther called on a non QOHOOH species.");
			continue;
		}
		std::string correction = "none";
		int posOOH = QOHOOH.trova(4);
		int posR = QOHOOH.trova(2);
		int posoh = QOHOOH.trova(10);
		if (QOHOOH.trova(11))
		{
			posR = QOHOOH.trova(11);
			posoh = QOHOOH.trova(11);
		}

		int dist = QOHOOH.dist(posR, posOOH);
		AnelloO ring;
		if (dist == 1)
			ring = ao3;
		else if (dist == 2)
		{
			ring = ao4;
			int alpha_ooh[SIZEMAX + 1];
			int alpha_rad[SIZEMAX + 1];
			QOHOOH.scorri(posOOH, 1, alpha_ooh);
			QOHOOH.scorri(posR, 1, alpha_rad);
			for (int i = 1; i <= SIZEMAX; i++)
			{
				if (alpha_ooh[i] == 1 && alpha_rad[i] == 1)
				{
					if (QOHOOH.tipoC(i) == Ct)
						correction = "T";
					if (QOHOOH.tipoC(i) == Cq)
						correction = "Q";
				}
			}
		}
		else if (dist == 3)
		{
			ring = ao5;
		}
		else if (dist == 4)
		{
			ring = ao6;
		}
		else
			continue;
		if (QOHOOH.tipoROOH(posOOH) == Rp && QOHOOH.tipoR(posR) == Rp && ring == ao3)
			continue;

		// This part is to generate the products based on QOOH;

		Molecola OHcEth = QOHOOH.parentFuel(); // removes/cleans oles, ethers, Os on a QOOH. 
		OHcEth.addetero(posR, posOOH);  // add the ether O between the two positions on the "cleaned" molec.

		std::string oh_ring = "No";
		if ((posoh >= posOOH && posoh <= posR) || (posoh >= posR && posoh <= posOOH))
			oh_ring = "Yes";

		reactionComment reacomm = k->v_qohooh_to_ether(QOHOOH.tipoC(posOOH),
			QOHOOH.tipoC(posR), ring, oh_ring);
		reactions.push_back(Reaction(QOHOOH, std::vector<Molecola>{ OHcEth, OH },
			new double[3] { k->A, k->n, k->E }, "Class23: QOHOOH to cyclic ethers reactions", reacomm));
	}
	return reactions;
}

std::vector<Reaction> alphaQOHOOHToKHP(std::vector<Molecola> QOHOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	Molecola HO2(5);
	for (auto& QOHOOH : QOHOOHs)
	{
		if (QOHOOH.kindOfSPecies() != QOHOOH_) // Jiaxin_ change from QOOH_ to oleQOOH_.
		{
			UTL::error("alphaQOHOOHToKHP called on a non QOHOOH species.");
			continue;
		}

		if (!QOHOOH.trova(11)) continue;
		int posR = QOHOOH.trova(11);
		int posoh = QOHOOH.trova(11);
		int posooh = QOHOOH.trova(4);

		Molecola KHP = QOHOOH;
		KHP.COHrad_to_C(posR);
		KHP.addcheto(posR);

		reactionComment reacomm = k->v_alpha_qohooh(QOHOOH.tipoR(posR));
		reactions.push_back(Reaction(std::vector<Molecola>{ QOHOOH, O2 }, std::vector<Molecola>{ KHP, HO2 },
			new double[3] { k->A, k->n, k->E }, "Class50: alpha-QOHOOH To KHPs", reacomm));
	}
	return reactions;
}


std::vector<Reaction> OHQOOHdecompositions(std::vector<Molecola> QOHOOHs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& QOHOOH : QOHOOHs)
	{
		if (QOHOOH.kindOfSPecies() != QOHOOH_)
		{
			UTL::error("OHQOOHdecompositions called on a non QOHOOH species.");
			continue;
		}
		int posR = QOHOOH.trova(2);
		if (QOHOOH.trova(11)) posR = QOHOOH.trova(11);
		int posOOH = QOHOOH.trova(4);
		int dist = QOHOOH.dist(posR, posOOH);
		if (dist == 2)
		{
			int pos_alfa;   // find position where to break bond
			{
				int alfa_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				QOHOOH.scorri(posOOH, 1, alfa_ooh);
				QOHOOH.scorri(posR, 1, alfa_rad);
				int j;
				for (j = 1; j <= SIZEMAX; j++)
					if (alfa_ooh[j] == 1 && alfa_rad[j] == 1)
						break;
				pos_alfa = j;
			}
			// break molecule
			Molecola m1, m2;
			int isPossible = QOHOOH.spezza1(posOOH, pos_alfa, &m1, &m2);
			if (isPossible == 0)
				continue;
			int posRprod = m1.trova(2);
			if (m1.trova(11)) posRprod = m1.trova(11);
			m1.Crad_to_C(posRprod);
			m1.addcheto(posRprod);

			std::vector<Molecola> prods;
			prods.push_back(m1);
			prods.push_back(m2);
			prods.push_back(OH);
			prods = fullyDecomposeRO1(prods);

			reactionComment reacomm = k->v_qohooh_decom("gamma");
			reactions.push_back(Reaction(QOHOOH, prods, new double[3] { k->A, k->n, k->E },
				"Class22_combine: QOHOOH decomposition reactions", reacomm));
		}

		if (dist == 3)
		{
			//                 find the position where to break
			int pos_break1;
			{
				int beta_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				QOHOOH.scorri(posOOH, 2, beta_ooh);
				QOHOOH.scorri(posR, 1, alfa_rad);
				int j;
				for (j = 1; j <= SIZEMAX; j++)
					if (beta_ooh[j] == 1 && alfa_rad[j] == 1) break;
				pos_break1 = j;
			};
			int pos_break2;
			{
				int alfa_ooh[SIZEMAX + 1];
				int beta_rad[SIZEMAX + 1];
				QOHOOH.scorri(posOOH, 1, alfa_ooh);
				QOHOOH.scorri(posR, 2, beta_rad);
				int j;
				for (j = 1; j <= SIZEMAX; j++)
					if (alfa_ooh[j] == 1 && beta_rad[j] == 1) break;
				pos_break2 = j;
			};
			//           generate decomposition products
			Molecola m1, m2;
			int isPossible = QOHOOH.spezza1(pos_break1, pos_break2, &m1, &m2);
			if (isPossible == 0)
				continue;
			
			reactionComment reacomm = k->v_qohooh_decom("delta");
			reactions.push_back(Reaction(QOHOOH, std::vector<Molecola>{m1, m2},
				new double[3] { k->A, k->n, k->E },
				"Class22_combine: QOHOOH decomposition reactions", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> O2AdditionToOHQOOH(std::vector<Molecola> QOHOOHs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& QOHOOH : QOHOOHs)
	{
		if (QOHOOH.kindOfSPecies() != QOHOOH_) // JIAXIN change from QOOH_ to oleQOOH_
		{
			UTL::error("O2AdditionToOHQOOH called on a non QOHOOH species");
			continue;
		}

		if (QOHOOH.trova(10) == QOHOOH.trova(2)) continue;

		int num_C = 0;
		Radicale radicalType;
		Molecola O2QOHOOH;
		if (QOHOOH.trova(2))
		{
			if (QOHOOH.tipoR(QOHOOH.trova(2)) == Rs) num_C = QOHOOH.numberOfC();
			O2QOHOOH = QOHOOH;
			O2QOHOOH.Crad_to_ROO(O2QOHOOH.trova(2));
			radicalType = QOHOOH.tipoR(QOHOOH.trova(2));
		}
		if (QOHOOH.trova(11))
		{
			if (QOHOOH.tipoR(QOHOOH.trova(11)) == Rs) num_C = QOHOOH.numberOfC();
			O2QOHOOH = QOHOOH;
			O2QOHOOH.Crad_to_ROO(O2QOHOOH.trova(11));
			radicalType = QOHOOH.tipoR(QOHOOH.trova(11));
		}
		reactionComment reacomm = k->v_o2_add_qohooh(radicalType, num_C); // JIAXIN, the C number dependency.
		reactions.push_back(Reaction(std::vector<Molecola>{ QOHOOH, O2 }, O2QOHOOH,
			new double[3] { k->A, k->n, k->E }, "Class24: O2 addition to QOHOOH reactions", reacomm));
	}
	return reactions;
}
std::vector<Reaction> O2QOHOOHElimination(std::vector<Molecola> O2QOHOOHs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	for (auto& O2QOHOOH : O2QOHOOHs)
	{
		if (O2QOHOOH.kindOfSPecies() != O2QOHOOH_)
		{
			UTL::error("O2QOHOOHElimination called on a non O2QOHOOH species.");
			continue;
		}
		int trovati[SIZEMAX + 1];
		int posOO = O2QOHOOH.trova(3);
		int pos_oh = O2QOHOOH.trova(10);
		//if (O2QOHOOH.trova(12))
		//{
		//	posOO = O2QOHOOH.trova(12);
		//	pos_oh = O2QOHOOH.trova(12);
		//}

		O2QOHOOH.scorri(posOO, 1, trovati);	// find the carbon at distance 1 from the carbon with the OO
		for (int j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
		{
			if (trovati[j] == 1 && O2QOHOOH.numAbstractableH(j) != 0
				&& O2QOHOOH.trova(4) != j)
			{
				if (O2QOHOOH.tipoROO(posOO) == Rp && O2QOHOOH.tipoC(j) == Cp)
					continue;
				Molecola OHoleOOH = O2QOHOOH;
				OHoleOOH.removeOO(posOO);
				OHoleOOH.addole(posOO, j);

				std::string oh_to_r = "No"; // OH and expelled H are not on the same C.
				if (j == pos_oh) oh_to_r = "Yes";

				reactionComment reacomm = k->v_o2qohooh_elim(O2QOHOOH.tipoC(posOO),
					O2QOHOOH.tipoC(j), oh_to_r, O2QOHOOH.numAbstractableH(j));
				reactions.push_back(Reaction(O2QOHOOH, std::vector<Molecola>{ OHoleOOH, HO2 },
					new double[3] { k->A, k->n, k->E },
					"Class25: O2QOHOOH HO2 elimination reactions", reacomm));
			}
		}
	}
	return reactions;
}
std::vector<Reaction> O2QOHOOHToKHP(std::vector<Molecola> O2QOHOOHs, Kinox* k)     // JIAXIN
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& O2QOHOOH : O2QOHOOHs)
	{
		if (O2QOHOOH.kindOfSPecies() != O2QOHOOH_)
		{
			UTL::error("O2QOHOOHToKHP called on a non O2QOHOOH species.");
			continue;
		}
		int posOO = O2QOHOOH.trova(3);		// find the position of the oo group
		int posOOH = O2QOHOOH.trova(4);		// find the position of the ooh group

		int dist = O2QOHOOH.dist(posOO, posOOH);	// find the distance between the oo 
		//and ooh group
		//from the distance find the type of ring formed during the reaction
		Anello ring = a5;
		switch (dist)
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
		default:		// if the ring is bigger than 7 atoms skip 
			continue;	// since the reaction cannot not happen
			break;
		}
		Molecola OHKHP = O2QOHOOH;
		OHKHP.COOrad_to_COOH(posOO);		// change the oo group in an ooh group
		OHKHP.tipo(posOOH, 1);					// remove the ooh group and ...
		int isPossible = OHKHP.addcheto(posOOH);	// ... replace it with the group =o 
		// Check if this step is possible.
		if (isPossible != 0)					// if it is possible save reaction
		{
			reactionComment reacomm = k->v_o2qohooh_to_khp(O2QOHOOH.tipoC(posOO),
				O2QOHOOH.tipoC(posOOH), ring, O2QOHOOH.numAbstractableH(posOOH));
			reactions.push_back(Reaction(O2QOHOOH, std::vector<Molecola>{ OHKHP, OH },
				new double[3] { k->A, k->n, k->E },
				"Class26: OH-KHP formation from O2QOHOOH reactions", reacomm));
		}
	}
	return reactions;
}
std::vector<Reaction> OHKHPdecomposition(std::vector<Molecola> OHKHPs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& OHKHP : OHKHPs)
	{
		if (OHKHP.kindOfSPecies() != OHKHP_)
		{
			UTL::error("OHKHPdecomposition called on a non OHKHP species.");
			continue;
		}
		int posOOH = OHKHP.trova(4);		// find the position of the ooh group
		int posO = OHKHP.trova(7);		// find the position of the =o group
		int dist = OHKHP.dist(posO, posOOH);	// find the distance between the ooh and =o group
		// from the distance find the type of ring formed during the reaction
		Anello ring = a5;
		switch (dist)
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
		default:		// if the ring is bigger than 7 atoms skip
			continue;
			break;
		}

		int pos_alpha;
		{
			int alpha_ooh[SIZEMAX + 1];
			int alpha_rad[SIZEMAX + 1];
			OHKHP.scorri(posO, dist - 1, alpha_ooh);
			OHKHP.scorri(posOOH, 1, alpha_rad);
			for (int j = 1; j <= SIZEMAX; j++)
			{
				if (alpha_ooh[j] == 1 && alpha_rad[j] == 1)
				{
					pos_alpha = j;
					break;
				}
			}
		}

		Molecola m1, m2;

		int isPossible = OHKHP.spezza1(pos_alpha, posOOH, &m1, &m2); // JIAXIN
		if (isPossible == 0)
			continue;
		int pos_rad_m2 = m2.trova(2);	// find the radical on the molecule m2
		m2.tipo(pos_rad_m2, 1);			// remove the radical and..
		//if (m2.trova(11))
		//{
		//	pos_rad_m2 = m2.trova(11);
		//	m2.tipo(pos_rad_m2, 10);
		//}
		
		m2.addcheto(pos_rad_m2);		// ... replace it with a cheto group
		
		
		std::vector<Molecola> products;

		if (m2.kindOfSPecies() == COHCHO_ && m2.size() > 2)
		{
			std::vector<Molecola> OHKHP_dec_prod = decomposeketo(m2);
			for (auto& pr1 : OHKHP_dec_prod)
				products.push_back(pr1);
		}
		else products = { m2 };

		if (m1.kindOfSPecies() == RO_ || m1.kindOfSPecies() == OHRO_ || m1.kindOfSPecies() == COHCHO_)
		{
			std::vector<Molecola> RO_dec_prod = fullyDecomposeRO(m1);
			for (auto& pr : RO_dec_prod)
				products.push_back(pr);
		}
		else
			products.push_back(m1);

		products.push_back(OH);

		reactionComment reacomm = k->v_ohkhp_decom(OHKHP.tipoROOH(posOOH), dist);

		reactions.push_back(Reaction(OHKHP, products, new double[3] { k->A, k->n, k->E },
			"Class27: OH-KHP decomposition reactions", reacomm));
	}
	return reactions;
}
std::vector<Reaction> O2QOHOOHIsomerization(std::vector<Molecola> O2QOHOOHs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	for (auto& O2QOHOOH : O2QOHOOHs)
	{
		if (O2QOHOOH.kindOfSPecies() != O2QOHOOH_)
		{
			UTL::error("O2QOHOOH called on a non O2QOHOOH species.");
			continue;
		}
		int j;
		int pos_o2 = O2QOHOOH.trova(3);
		int pos_oh = O2QOHOOH.trova(10);
		int trovati[SIZEMAX + 1];
		for (int distance = 1; distance <= 4; distance++)
		{
			//std::cout << "  Distance = " << distance << std::endl;
			O2QOHOOH.scorri(pos_o2, distance, trovati);
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
				continue;
				break;
			}

			for (j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
			{
				if (trovati[j] == 1 && O2QOHOOH.numAbstractableH(j) != 0
					&& O2QOHOOH.tipo(j) != 4)
				{
					//if (pos_oh == j) continue; //JLMOD0213
					reactionComment reacomm = k->v_o2qohooh_isom(O2QOHOOH.tipoC(O2QOHOOH.trova(3)), O2QOHOOH.tipoC(j), 
						ring, O2QOHOOH.numAbstractableH(j));  // JIAXIN
					Molecola POOH2 = O2QOHOOH;
					POOH2.OO_to_OOH(O2QOHOOH.trova(3));
					POOH2.C_to_Crad(j);

					reactions.push_back(Reaction(O2QOHOOH, POOH2,
						new double[3] { k->A, k->n, k->E }, "Class28: O2QOHOOH isomerization reactions",
						reacomm));
				}
			}
		}
	}
	return reactions;

}

std::vector<Reaction> POHOOHdecomposition(std::vector<Molecola> POHOOHs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	Molecola OH(4);
	for (auto& POHOOH : POHOOHs)
	{
		if (POHOOH.kindOfSPecies() != POHOOH_)
		{
			UTL::error("POHOOHdecomposition called on a non POHOOH  species.");
			continue;
		}
		std::vector<int> posOOHs = POHOOH.posOOHinPOOH2();
		int posR = POHOOH.trova(2);
		for (auto& posOOH : posOOHs)
		{
			if (POHOOH.dist(posOOH, posR) == 1)
			{
				if (POHOOH.tipoROOH(posOOH) == Rp && POHOOH.tipoR(posR) == Rp)
					continue;
				Molecola OHOLEOOH = POHOOH;
				OHOLEOOH.removeOOH(posOOH);
				OHOLEOOH.Crad_to_C(posR);
				OHOLEOOH.addole(posOOH, posR);
				reactionComment reacomm = k->v_ohpooh_decom(POHOOH.tipoC(posOOH),
					POHOOH.tipoR(posR), 1);
				reactions.push_back(Reaction(POHOOH, std::vector<Molecola>{ OHOLEOOH, HO2 },
					new double[3] { k->A, k->n, k->E }, "Class29: POHOOH decomposition reactions",
					reacomm));
			}
			if (POHOOH.dist(posOOH, posR) == 2)
			{
				//                  find the position where to break
				int pos_alfa;
				{
					int alfa_ooh[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					POHOOH.scorri(posOOH, 1, alfa_ooh);
					POHOOH.scorri(posR, 1, alfa_rad);
					int j;
					for (j = 1; j <= SIZEMAX; j++)
						if (alfa_ooh[j] == 1 && alfa_rad[j] == 1) break;
					pos_alfa = j;
				}
				//                             generate the broken molecules
				if (POHOOH.tipo(pos_alfa) == 4)
					continue;
				Molecola m1, m2;
				int isPossible = POHOOH.spezza(posOOH, pos_alfa, &m1, &m2);
				if (isPossible == 0) continue;

				int pos = m1.trova(2);
				m1.tipo(pos, 1);
				m1.addcheto(pos);

				std::vector<Molecola> prods;
				prods.push_back(m1);
				prods.push_back(m2);
				prods.push_back(OH);
				prods = fullyDecomposeRO1(prods);


				reactionComment reacomm = k->v_ohpooh_decom(POHOOH.tipoC(posOOH),
					POHOOH.tipoR(posR), 2);
				reactions.push_back(Reaction(POHOOH, prods,
					new double[3] { k->A, k->n, k->E }, "Class29: POHOOH decomposition reactions",
					reacomm));
			}
			int pos_rad = -1;		// find the position where the radical will be located
			{
				int beta_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				POHOOH.scorri(posOOH, 1, beta_ooh);
				POHOOH.scorri(posR, 2, alfa_rad);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					if (beta_ooh[j] == 1 && alfa_rad[j] == 1)
					{
						pos_rad = j;
						break;
					}
				}
			}
			if (pos_rad == -1)
				continue;

			if (POHOOH.tipo(pos_rad) == 4)  // if there is a COOH in this position continue
				continue;

			//                  find the position where to break
			int pos_break1 = -1;
			{
				int beta_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				POHOOH.scorri(posOOH, 2, beta_ooh);
				POHOOH.scorri(posR, 1, alfa_rad);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					if (beta_ooh[j] == 1 && alfa_rad[j] == 1)
					{
						pos_break1 = j;
						break;
					}
				}
			}
			if (pos_break1 == -1)
				continue;

			int pos_break2 = -1;
			{
				int alfa_ooh[SIZEMAX + 1];
				int beta_rad[SIZEMAX + 1];
				POHOOH.scorri(posOOH, 1, alfa_ooh);
				POHOOH.scorri(posR, 2, beta_rad);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					if (alfa_ooh[j] == 1 && beta_rad[j] == 1)
					{
						pos_break2 = j;
						break;
					}
				}
			}
			if (pos_break2 == -1)
				continue;

			if (POHOOH.tipo(pos_break1) == 4)
				continue;

			//                     break the molecule
			if (pos_break1 == pos_break2)
				continue;
			if (POHOOH.areBonded(pos_break1, pos_break2) == 0)
				continue;

			Molecola m1, m2;
			int isPossible = POHOOH.spezza(pos_break1, pos_break2, &m1, &m2);
			if (isPossible == 0)
				continue;

			reactionComment reacomm = k->v_ohpooh_decom(POHOOH.tipoC(posOOH),
				POHOOH.tipoR(posR), 3);
			Reaction reac(POHOOH, std::vector<Molecola>{ m1, m2},
				new double[3] { k->A, k->n, k->E }, "Class29: POHOOH decomposition reactions", reacomm);
			UTL::addUnique(&reactions, reac);

		}
	}
	return reactions;

}

std::vector<Reaction> POHOOHToEthers(std::vector<Molecola> POHOOHs, Kinox* k)      // JIAXIN
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& POHOOH : POHOOHs)
	{
		if (POHOOH.kindOfSPecies() != POHOOH_)
		{
			UTL::error("POHOOHToEthers called on a non POHOOH species.");
			continue;
		}

		std::vector<int> posOOHs = POHOOH.posOOHinPOOH2();
		int posR = POHOOH.trova(2);
		for (auto& posOOH : posOOHs)
		{
			int dist = POHOOH.dist(posR, posOOH);
			AnelloO ring = ao3;
			switch (dist)
			{
			case 1:
				ring = ao3;
				break;
			case 2:
				ring = ao4;
				break;
			case 3:
				ring = ao5;
				break;
			case 4:
				ring = ao6;
				break;
			default:
				continue;
				break;
			}
			if (POHOOH.tipoROOH(posOOH) == Rp && POHOOH.tipoR(posR) == Rp && ring == ao3)
				continue;
			Molecola OHcEthOOH = POHOOH;
			OHcEthOOH.removeOOH(posOOH);
			OHcEthOOH.Crad_to_C(posR);
			OHcEthOOH.addetero(posOOH, posR);
			OHcEthOOH.kindOfSPecies();
			reactionComment reacomm = k->v_ohpooh_cEther(POHOOH.tipoC(posOOH),
				POHOOH.tipoR(posR), ring);
			reactions.push_back(Reaction(POHOOH, std::vector<Molecola>{ OHcEthOOH, OH },
				new double[3] { k->A, k->n, k->E }, "Class30: POHOOH to ether-OHOOHs reactions", reacomm));
		}
	}
	return reactions;
}


































int getAdditionalProducts(std::vector<Molecola>* molVec,
	std::vector<Reaction> reactions, species kind)
{
	int newElemCount = 0;
	for (auto& reac : reactions)
	{
		std::vector<Molecola*> singProd = reac.productList();
		for (auto& prod : singProd)
		{
			if (prod->kindOfSPecies() == kind)
				if (UTL::addUnique(molVec, *prod) == -1)
					newElemCount++;
		}
	}
	return newElemCount;
}

std::vector<Molecola> getProducts(std::vector<Reaction> reactions, species kind)
{
	std::vector<Molecola> products = {};
	getAdditionalProducts(&products, reactions, kind);
	return products;
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

void addNewSpecies(std::vector<Molecola>* molVec, std::vector<Reaction>* reacVec)
{
	for (auto& reac : *reacVec)
	{
		std::vector<Molecola*> reactants = reac.reactantList();
		for (auto& spec : reactants)
			UTL::addUnique(molVec, *spec);
		std::vector<Molecola*> products = reac.productList();
		for (auto& spec : products)
			UTL::addUnique(molVec, *spec);
	}
}

int processDuplicateReactions(std::vector<Reaction>* reacVec)
{
	std::vector<Reaction> newReac;
	std::vector<bool> toSkip(reacVec->size(), false);
	int mergedReactions = 0;
	for (int i = 0; i < reacVec->size(); i++)
	{
		if (toSkip[i])
			continue;
		Reaction reac = (*reacVec)[i];
		int numDuplicates = 0;
		for (int j = i + 1; j < reacVec->size(); j++)
		{
			if (toSkip[j])
				continue;
			if ((*reacVec)[i] == (*reacVec)[j])
			{
				numDuplicates++;
				toSkip[j] = true;
			}
		}
		if (numDuplicates > 0)
		{
			reac.setMultiplePathMultiplier(numDuplicates + 1);
			reac.setA(reac.A() * double(numDuplicates + 1));
			mergedReactions += numDuplicates;
		}
		newReac.push_back(reac);
	}
	for (int i = 0; i < newReac.size(); i++)
	{
		for (int j = i + 1; j < newReac.size(); j++)
		{
			if (newReac[i].weakEquality(newReac[j]))
			{
				newReac[i].setDuplicate();
				newReac[j].setDuplicate();
			}
		}
	}
	(*reacVec) = newReac;
	return mergedReactions;
}

void printSpeciesInFile(std::ofstream* outfile, std::vector<Molecola> mols,
	std::string label, ChemkinOut* chemOut)
{
	if (mols.size() > 0)
	{
		*outfile << "! " << label << std::endl;
		for (int i = 0; i < mols.size(); i++)
		{
			*outfile << std::left;
			*outfile << std::setw(16) << chemOut->molToName(mols[i])
				<< std::setw(3) << " ";
			if ((i + 1) % 4 == 0 && i != mols.size() - 1)
				*outfile << std::endl;
		}
		*outfile << std::endl;
		*outfile << std::endl;
	}
}

void printSpeciesInFile(std::ofstream* outfile, std::vector<std::string> mols,
	std::string label)
{
	if (mols.size() > 0)
	{
		*outfile << "! " << label << std::endl;
		for (int i = 0; i < mols.size(); i++)
		{
			*outfile << std::left;
			*outfile << std::setw(16) << mols[i]
				<< std::setw(3) << " ";
			if ((i + 1) % 4 == 0 && i != mols.size() - 1)
				*outfile << std::endl;
		}
		*outfile << std::endl;
		*outfile << std::endl;
	}
}

void printReaction(std::ofstream* outfile, Reaction reac, ChemkinOut* chemOut)
{
	//std::cout << reac << std::endl;
	std::vector<Molecola*> reactants = reac.reactantList();
	std::vector<Molecola*> products = reac.productList();

	std::stringstream reactantsSide;
	for (int i = 0; i < reactants.size(); i++)
	{
		reactantsSide << chemOut->molToName(*(reactants[i]));
		if (i != reactants.size() - 1)
			reactantsSide << " + ";
	}

	std::stringstream productsSide;
	for (int i = 0; i < products.size(); i++)
	{
		productsSide << chemOut->molToName(*(products[i]));
		if (i != products.size() - 1)
			productsSide << " + ";
	}

	char cost[50];
	sprintf_s(cost, 50, "   %9.3e %8.3f %8.0f ", reac.A(), reac.n(), reac.E());

	*outfile << std::left << std::setw(23) << reactantsSide.str();
	if(reac.isReversible() == false)
		*outfile << " => ";
	else
		*outfile << " <=> ";
	*outfile << std::setw(40) << productsSide.str() << " "
		<< std::setw(30) << cost << " "
		<< "!" << reac.printComment() << std::endl;
	if (reac.isDuplicate())
		*outfile << "DUP" << std::endl;

}

void printReactions(std::ofstream* outfile, std::vector<Reaction> reacs,
	ChemkinOut* chemOut)
{
	for (auto& rea : reacs)
		printReaction(outfile, rea, chemOut);
}

bool isThereDecompositionPath(Molecola spec,
	std::vector<baseReaction>* baseMechReacs,
	std::vector<Reaction>* totReactions, ChemkinOut* chemOut)
{
	std::string specName = chemOut->molToName(spec);
	for (auto& baseReac : *baseMechReacs)
	{
		for (auto& reactant : baseReac.reactants)
			if (specName == reactant)
				return true;
		if (baseReac.isReversible)
			for (auto& product : baseReac.products)
				if (specName == product)
					return true;
	}
	for (auto& reaction : *totReactions)
	{
		std::vector<Molecola*> reactants = reaction.reactantList();
		for (auto& reac : reactants)
			if (spec == *reac)
				return true;
		//if (reaction.isReversible())
		//{
		//	std::vector<Molecola*> products = reaction.productList();
		//
		//	for (auto& product : products)
		//		if (spec == *product)
		//			return true;
		//}
	}
	return false;
}

std::string speciesToText(species kind)
{
	switch (kind)
	{
	case fuel_:
		return "HC";
		break;
	case R_:
		return "R";
		break;
	case ROO_:
		return "ROO";
		break;
	case QOOH_:
		return "QOOH";
		break;
	case OOQOOH_:
		return "OOQOOH";
		break;
	case OLE_:
		return "OLE";
		break;
	case CO_:
		return "ALD/KETO";
		break;
	case cEth_:
		return "cEth";
		break;
	case RO_:
		break;
		return "RO";
	case KHP_:
		break;
		return "KHP";
	case ROOH_:
		return "ROO";
		break;
	case POOH2_:
		return "P(OOH)2";
		break;
	case ROOH2_:
		return "R(OOH)2";
		break;
	case oleOOH_:
		return "OLE-OOH";
		break;
	case cEthOOH_:
		return "cEth-OOH";
		break;
	case oleR_:
		return "OLE-R";
		break;
	case oleCO_:
		return "OLE-CO";
		break;
	case cEthR_:
		return "cEth-R";
		break;
	case cEthCO_:
		return "cEth-CO";
		break;
	case lEthRO_:
		return "lEth-RO";
		break;
	case alkRO_:
		return "alkRO";
		break;
		// JIAXIN
		//case dienes_:
		//	oleRO2_elim.push_back(mol);
		//	break;
		//case dieneOOH_:
		//	oleO2QOOH_elim.push_back(mol);
		//	break;
	case oleKHP_:
		return "oleKHP";
		break;
	case dienes_:
		return "dienes";
		break;
	case oleOH_:
		return "oleOH";
		break;
	case oleRO_:
		return "oleRO";
		break;
	case OHoleOOH_:
		return "OHoleOOH";
		break;
	case OHoleketo_:
		return "OHoleketo";
		break;
	case oleROH_:
		return "oleROH";
		break;
	case dienesOH_:
		return "dienesOH";
		break;
	case ROH_:
		return "ROH";
		break;
	case dienesR_:
		return "dienesR";
		break;
	case diketo_:
		return "diketo";
		break;
	case dienesCO_:
		return "dienesCO";
		break;
	case OHRO_:
		return "OHRO";
		break;
	case QOHOOH_:
		return "QOHOOH";
		break;
	case COHCHO_:
		return "COHCHO";
		break;
	case OHKHP_:
		return "OHKHP";
		break;
	case special_:
		return "special";
		break;
	case unidentified_:
		return "unidentified";
		break;

	}
	return "ERR";
}

bool isIncluded(Reaction reac, std::vector<baseReaction>* reacList,
	ChemkinOut* chemOut)
{
	std::vector<Molecola*> reactants = reac.reactantList();
	std::vector<Molecola*> products = reac.productList();
	std::vector<std::string> reactantsNames;
	std::vector<std::string> productsNames;
	for (auto& spec : reactants)
		reactantsNames.push_back(chemOut->molToName(*spec));
	for (auto& spec : products)
		productsNames.push_back(chemOut->molToName(*spec));

	for (auto& baseReac : *reacList)
	{
		if (baseReac.isReversible == false)
		{
			if (UTL::areEquivalent(&reactantsNames, &(baseReac.reactants))
				&& UTL::areEquivalent(&productsNames, &(baseReac.products)))
				return true;
		}
		else
		{
			if ((UTL::areEquivalent(&reactantsNames, &(baseReac.reactants))
				&& UTL::areEquivalent(&productsNames, &(baseReac.products)))
				|| (UTL::areEquivalent(&productsNames, &(baseReac.reactants))
					&& UTL::areEquivalent(&reactantsNames, &(baseReac.products))))
				return true;
		}
	}
	return false;
}

bool isSpeciesIncluded(Molecola spec, ChemkinOut* chemOut)
{
	std::string molInchi = spec.inchiName();

	std::string InC = chemOut->checkIfNamed(&spec);
	if (InC == "not found") return false;
	return true;
}


std::string expandEnvVar(std::string inString)
{
	const char * inCh = inString.c_str();
	char outCh[500];
	ExpandEnvironmentStringsA(inCh, outCh, 500);
	std::string outSt(outCh);
	return outSt;
}