#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdarg.h>
#include "molec.h"
#include <vector>
#include "RateRule.h"
#include "ReactionComment.h"

//using namespace std;

#ifndef ENUM
#define ENUM
enum Carbonio { Cp, Cs, Ct, Cq };
enum Idrogeno { Hp, Hs, Ht, Hcooh };
enum Radicale { Rpmet, Rpet, Rp, Rs, Rt };
enum Anello { a5, a6, a7 };
enum AnelloO { ao3, ao4, ao5 };
enum Direz { dir, inv };
enum Output { MATRICE, FORMULA };
enum HAbsRad { o2, oh, h, o, ho2, ch3, c2h5, ch3o, ch3o2, ch3oh, ch3o2h}; // JIAXIN
#endif

class Kinox
{
	
public:

	RateRule initiation{ 3, true };  // 2 to 3. // first value is the number of dependencies (or parameters) in csv, the second value is whether the rate rule is symmetrical or not
	RateRule hAbstraction{ 3 ,false };
	RateRule isomerizationR{ 4 ,false }; // 3 to 4
	RateRule betaDecR{ 3 ,false }; // Jiaxin 2 to 3
	RateRule oleFromR{ 0 ,false };
	RateRule O2AdditionR{ 2 ,false }; // 2 parameters. JIAXIN Change it to fit Olefin rate rule.
	RateRule OHAddOle{ 2, false }; // Jiaxin 2 params.
	RateRule O2AddROH{ 1, false }; // Jiaxin 1 param.
	RateRule ROHIsomReac{ 3, false }; // Jiaxin 1 param.
	RateRule ROHO2IsomReac{ 4, false }; // Jiaxin 1 param.

	RateRule HAddOle{2, false}; //
	RateRule HO2AddAllyR_rec{1, false};
	RateRule HO2AddAllyR_dec{ 1, false };
	RateRule ch3AddAllyR{ 1, false };
	RateRule ho2AddOle{ 2, false };

	RateRule OatomAddOle{ 2, false };
	RateRule rooh_to_ro{ 1, false };
	RateRule ro_betadecom{ 1, false };
	RateRule o2addvinyR{ 1, false };
	RateRule o2habsallyR{ 1, false };

	RateRule rohoo_wad{ 2, false };
	RateRule wad_decom{ 2, false };
	RateRule rohoo_elim{ 3, false };
	RateRule bqohooh_decom{ 3, false };
	RateRule qohoohtocEther{ 4, false };
	RateRule qohooh_decom{ 1, false };
	RateRule o2addqohooh{ 2, false };
	RateRule o2qohoohelim{ 3, false };
	RateRule o2qohoohTokhp{ 3, false };
	RateRule ohkhpdecom{ 2, false };
	RateRule o2qohoohIsom{ 3, false };
	RateRule ohpoohdecom{ 3, false };
	RateRule ohpoohtoEther{ 3, false };
	RateRule O2RemovalROO{ 1 ,false };
	RateRule isomROO{ 5 ,false }; // JIAXIN. 3 to 5 params.
	RateRule isomOOQOOH{ 5 ,false }; // JIAXIN, 3 to 5 params.
	RateRule isomQOOH{ 3 ,false };
	RateRule isomPOOH2{ 3 ,false };
	RateRule O2AdditionQOOH{ 3 ,false }; // JIAXIN change from 1 to 3.
	RateRule O2RemovalOOQOOH{ 1 ,false };
	RateRule OOQOOHToKHP{ 5 ,false }; // Jiaxin no change. change from 3 to 5. 
	RateRule KHPDecomp{ 2 ,false }; // change from 1 to 2.
	RateRule oleFromROO{ 4 ,false };  // JIAXIN 2 to 4.
	RateRule oleFromBetaQOOH{ 2 ,false };
	RateRule oleFromGammaQOOH{ 2 ,false };
	RateRule oleFromDeltaQOOH{ 2 ,false };
	RateRule POOH2Dec1{ 2 ,false };
	RateRule POOHDec2{ 2 ,false };
	RateRule POOHDec3{ 2 ,false };
	RateRule etherFromQOOH{ 6 ,false }; // Jiaxin 4 to 6
	RateRule oleqooh_decom{ 3, false };
	RateRule etherFromPOOH2{ 3 ,false };
	RateRule oleFromOOQOOH{ 4 ,false }; // JIAXIN 3 to 4.
	RateRule oleOOHDec{ 1 ,false };
	RateRule etherOOHDec{ 1 ,false };
	RateRule cycEthDec{ 0 ,false };
	RateRule allylicRadForm{ 1 ,false };
	RateRule alkenylROForm{ 0 ,false };
	RateRule alkenylRODec{ 1 ,false };
	RateRule aldDec{ 0 ,false };
	RateRule oleOHsConsum{ 0, false };
	RateRule Alpha_qohooh { 1, false };



public:

	double A;
	double n;
	double E;
	int sizeReactions = 0;
	Kinox(void);    // default constructor
	Kinox(std::string nome);    // constructor
	void leggi(char[80]);

	friend std::ostream& operator<<(std::ostream&, Kinox&);
	friend std::istream& operator>>(std::istream&, Kinox&);

	reactionComment v_initiation(Radicale r1, int isomers, int IsRAlly); // JIAXIN.
	reactionComment v_h_abstraction(Molecola r, Carbonio c, int isAllorVyn, int numH, int isomers, std::string OH_label);
	reactionComment v_isomerization_r(Radicale r, Idrogeno h, Anello a, int numH, int Hally);
	// JIAXIN
	reactionComment v_roh_isom_r(Radicale r, Idrogeno h, Anello a, int numH); //betaROH isomerization.
	reactionComment v_beta_dec_r(Radicale r1, Radicale r2, int isAllyOrNot);
	reactionComment v_ole_par_r(int numH);
	reactionComment v_o2_add_r(Radicale r, int isAllyOrNot); // JIAXIN
	reactionComment v_o2_add_roh(Radicale r); // JIAXIN
	reactionComment v_alpha_qohooh(Radicale r);
	reactionComment v_oh_add_ole(Carbonio OHAdd_site, Carbonio radical_site); // JIAXIN OH add

	reactionComment v_h_add_ole(Carbonio HAdd_site, Carbonio radical_site); // JIAXIN OH add
	reactionComment v_ho2_add_allyR_rec(Radicale radical_site); // JIAXIN
	reactionComment v_ho2_add_allyR_dec(Radicale radical_site); // JIAXIN
	reactionComment v_ch3_add_allyR(); // JIAXIN
	reactionComment v_ho2_add_ole(Carbonio HO2Add_site, Carbonio radical_site); // JIAXIN OH add

	reactionComment v_o_add_ole(Carbonio OAdd_site, std::string R_expel); // JIAXIN OH add
	reactionComment v_rooh_to_ro(); // JIAXIN OH add
	reactionComment v_ro_beta(std::string IsViny); 
	reactionComment v_o2_add_vinyR(Radicale radical_site); // JIAXIN OH add
	reactionComment v_o2_habs_allyR(Carbonio c1, int numH, int isomers);

	reactionComment v_beta_rohoo_wad(Carbonio COO_site, Carbonio COH_site);
	reactionComment v_wad_decom(Carbonio COOH_site, Carbonio CO_site);
	reactionComment v_beta_rohoo_elim(Carbonio COO_site, Carbonio H_site, int numH, std::string oh_to_r);
	reactionComment v_beta_qohooh_decom(Radicale COOH_site, Radicale r_site, std::string oh_to_r);
	reactionComment v_qohooh_to_ether(Carbonio r1, Carbonio r2, AnelloO a, std::string oh_ring);
	reactionComment v_qohooh_decom(std::string ooh_to_r);
	reactionComment v_o2_add_qohooh(Radicale r, int numC);
	reactionComment v_o2qohooh_elim(Carbonio COO_site, Carbonio H_site, std::string oh_to_hsite, int numH);
	reactionComment v_o2qohooh_to_khp(Carbonio COO_site, Carbonio COOH_site, Anello ring_size, int numH);
	reactionComment v_ohkhp_decom(Radicale COOH_site, int dist);
	reactionComment v_o2qohooh_isom(Carbonio COO_site, Carbonio H_site, Anello ring_size, int numH);
	reactionComment v_ohpooh_decom(Carbonio COOH_site, Radicale radical_site, int dist);
	reactionComment v_ohpooh_cEther(Carbonio COOH_site, Radicale radical_site, AnelloO ring_size);

	reactionComment v_o2_rem_roo(Radicale r);
	reactionComment v_isom_roo(Radicale r, Idrogeno h, Anello a, int numH, int Hally, int OOally); // Jiaxin, add 2 more conditions.
	reactionComment v_isom_rohoo(Radicale r, Idrogeno h, Anello a, int numH, int OHonRing); // Jiaxin, add 2 more conditions.
	reactionComment v_isom_ooqooh(Radicale r, Idrogeno h, Anello a, int numH, int Hally, int OOally); // Jiaxin, add 2 more conditions.
	reactionComment v_isom_qooh(Radicale r1, Radicale r2, Anello a);
	reactionComment v_isom_pooh2(Radicale r1, Radicale r2, Anello a);
	reactionComment v_o2_add_qooh(Radicale r, int isAllyOrNot, int numC);  // Jiaxin, add the C number dependency.
	reactionComment v_o2_rem_ooqooh(Radicale r);
	reactionComment v_ooqooh_to_khp(Radicale r1, Radicale r2, Anello a, int numH, int roo_type, int rooh_type);
	reactionComment v_khp_decomp(Radicale r, int dist);
	reactionComment v_ole_par_roo(Radicale r1, Carbonio c, int numH, int OOally, int Hally);
	reactionComment v_ole_par_ooqooh(Radicale r1, Carbonio c, int numH, int OOally, int Hally);
	reactionComment v_oleqooh_decom(Radicale r1, Radicale r2, int dist);
	//reactionComment v_ole_from_gamma_qooh(Radicale r1, Radicale r2);
	//reactionComment v_ole_from_delta_qooh(Radicale r1, Radicale r2);
	reactionComment v_pooh2_decom(Radicale r1, Radicale r2, int dist);
	//reactionComment v_pooh2_dec_2(Radicale r1, Radicale r2);
	//reactionComment v_pooh2_dec_3(Radicale r1, Radicale r2);
	reactionComment v_ether_from_qooh(Radicale r1, Radicale r2, AnelloO a, 
		std::string correction, int OOHally, int Rally);
	reactionComment v_ether_from_pooh2(Carbonio r1, Radicale r2, AnelloO a);
	reactionComment v_ohpooh_to_Ether(Carbonio r1, Radicale r2, AnelloO a);
	reactionComment v_oleooh_dec(Radicale r);
	reactionComment v_etherooh_dec(Radicale r);
	reactionComment v_cyc_eth_dec();
	reactionComment v_allylic_rad_form(Carbonio c, int numH);
	reactionComment v_alkenyl_ro_form();
	reactionComment v_alkenyl_ro_dec(Carbonio c);
	reactionComment v_ald_dec();
	reactionComment v_oleOHsConsum();

	//void setT(double Temp);
	//void wrireaLump(std::ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	//void wrireaLump1(std::ofstream& stream, Molecola reag, HAbsRad rad, int numC, double A, double n, double E);
	//void wrireaLump3(std::ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	//void wrireaLump4(std::ofstream& stream, int numC, double A, double n, double E);     // R'numC' + O2 -> R'numC'OO
	//void wrireaLump5(std::ofstream& stream, int numC, double A, double n, double E);     // R'numC'OO -> R'numC' + O2
	//void wrireaLump6(std::ofstream& stream, int numC, double A, double n, double E);     // R'numC' + O2 -> OLE'numC' + HO2
	//void wrireaLump7(std::ofstream& stream, int numC, double A, double n, double E);     // R'numC'OO -> Q'numC'OOH
	//void wrireaLump8(std::ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH -> R'numC'OO
	//void wrireaLump9(std::ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH -> ETER'numC' + OH
	//void wrireaLump10(std::ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH -> OLE'numC' + HO2
	//void wrireaLump11(std::ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	//void wrireaLump11b(std::ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff1, double A1, double n1, double E1,
	//	std::vector<double> stoicCoeff2, double A2, double n2, double E2);
	//void wrireaLump12(std::ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH + O2 -> OOQ'numC'OOH
	//void wrireaLump13(std::ofstream& stream, int numC, double A, double n, double E);     // OOQ'numC'OOH -> Q'numC'OOH + O2
	//void wrireaLump14(std::ofstream& stream, int numC, double A, double n, double E);     // OOQ'numC'OOH -> OQ'numC'OOH + OH
	//void wrireaLump15(std::ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);

	std::string nameHAbsRad(Molecola r);				// return a string containing the name of the radical rad
	std::string nameHAbsRadPlusH(Molecola r);			// return a string containing the name of the radical rad when an H is added to it
	void wrirea(std::ofstream& stream, int tiporeaz, double A, double n, double E, Molecola m1, ...);

};




