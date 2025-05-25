#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdarg.h>
#include "molec.h"
#include <vector>

using namespace std;

#ifndef ENUM
#define ENUM
enum Carbonio { Cp, Cs, Ct, Cq };
enum Idrogeno { Hp, Hs, Ht, Hcooh };
enum Radicale { Rpmet, Rpet, Rp, Rs, Rt };
enum Anello { a5, a6, a7 };
enum AnelloO { ao3, ao4, ao5 };
enum Direz { dir, inv };
enum Output { MATRICE, FORMULA };
enum HAbsRad { o2, oh, h, o, ho2, ch3, c2h5 };
#endif

class Kinox
{
	//protected:
public:
	// costanti da leggere da file
	double A_initiation[5][5];				// initiation
	double n_initiation[5][5];
	double E_initiation[5][5];
	double A_h_abstraction[8][3][15];			// h abstraction
	double n_h_abstraction[8][3][15];
	double E_h_abstraction[8][3][15];
	double A_isometization_R[3][3][3];		// R isomerization
	double n_isometization_R[3][3][3];
	double E_isometization_R[3][3][3];
	double A_beta_dec_R[3][5];				// beta decomposition R
	double n_beta_dec_R[3][5];
	double E_beta_dec_R[3][5];
	double A_ole_from_R;					// olefins from R + O2
	double n_ole_from_R;
	double E_ole_from_R;
	double A_O2_addition_R[3];				// O2 addition to R
	double n_O2_addition_R[3];
	double E_O2_addition_R[3];
	double A_O2_removal_ROO[3][4];			// O2 removal from ROO
	double n_O2_removal_ROO[3][4];
	double E_O2_removal_ROO[3][4];
	double A_isom_ROO[3][3][4];				// isomerization ROO
	double n_isom_ROO[3][3][4];
	double E_isom_ROO[3][3][4];
	double A_isom_OOQOOH[3][3][4];			// isomerization OOQOOH
	double n_isom_OOQOOH[3][3][4];
	double E_isom_OOQOOH[3][3][4];
	double A_isom_QOOH[3][3][3][4];			// isomerization QOOH
	double n_isom_QOOH[3][3][3][4];
	double E_isom_QOOH[3][3][3][4];
	double A_isom_POOH2[3][3][3];			// isomerization POOH2
	double n_isom_POOH2[3][3][3];
	double E_isom_POOH2[3][3][3];
	double A_O2_addition_QOOH[3];			// O2 addition to QOOH
	double n_O2_addition_QOOH[3];
	double E_O2_addition_QOOH[3];
	double A_O2_removal_OOQOOH[3][4];		// O2 removal from OOQOOH
	double n_O2_removal_OOQOOH[3][4];
	double E_O2_removal_OOQOOH[3][4];
	double A_OOQOOH_to_KHP[3][3][4];		// conversion of OOQOOH to KHP
	double n_OOQOOH_to_KHP[3][3][4];
	double E_OOQOOH_to_KHP[3][3][4];
	double A_KHP_decomp[3][3];				// KHP decomposition
	double n_KHP_decomp[3][3];
	double E_KHP_decomp[3][3];
	double A_ole_from_ROO[3][3];			// olefins from ROO
	double n_ole_from_ROO[3][3];
	double E_ole_from_ROO[3][3];
	double A_ole_from_beta_QOOH[3][3];		// olefins from beta QOOH
	double n_ole_from_beta_QOOH[3][3];
	double E_ole_from_beta_QOOH[3][3];
	double A_ole_from_gamma_QOOH[3][3];		// olefins from gamma QOOH
	double n_ole_from_gamma_QOOH[3][3];
	double E_ole_from_gamma_QOOH[3][3];
	double A_ole_from_delta_QOOH[3][3];		// olefins from delta QOOH
	double n_ole_from_delta_QOOH[3][3];
	double E_ole_from_delta_QOOH[3][3];
	double A_POOH2_dec_1[3][3];				// decomposition P(OOH)2 1
	double n_POOH2_dec_1[3][3];
	double E_POOH2_dec_1[3][3];
	double A_POOH2_dec_2[3][3];				// decomposition P(OOH)2 2
	double n_POOH2_dec_2[3][3];
	double E_POOH2_dec_2[3][3];
	double A_POOH2_dec_3[3][3];				// decomposition P(OOH)2 3
	double n_POOH2_dec_3[3][3];
	double E_POOH2_dec_3[3][3];
	double A_ether_from_QOOH[3][3][4][3];	// ethers from QOOH
	double n_ether_from_QOOH[3][3][4][3];
	double E_ether_from_QOOH[3][3][4][3];
	double A_ether_from_POOH2[3][3][4];		// ethers from POOH2
	double n_ether_from_POOH2[3][3][4];
	double E_ether_from_POOH2[3][3][4];
	double A_QOOH_beta_dec[3][3];			// QOOH beta decomposition
	double n_QOOH_beta_dec[3][3];
	double E_QOOH_beta_dec[3][3];
	double A_OH_abs[3];						// OH abstraction (reaction used for R species velocity)
	double n_OH_abs[3];
	double E_OH_abs[3];
	double A_ole_from_OOQOOH[3][3];			// olefin from OOQOOH
	double n_ole_from_OOQOOH[3][3];
	double E_ole_from_OOQOOH[3][3];
	double A_OOQOOH_dec_tert;				// OOQOOH homolytic decomposition tertiary sites
	double n_OOQOOH_dec_tert;
	double E_OOQOOH_dec_tert;
	double A_OOQOOH_dec_quat;				// OOQOOH homolytic decomposition quaternary sites
	double n_OOQOOH_dec_quat;
	double E_OOQOOH_dec_quat;
	double A_ROO_h_abs;						// h abstraction by generic ROO
	double n_ROO_h_abs;
	double E_ROO_h_abs;
	double A_R_h_abs;						// h abstraction by generic R
	double n_R_h_abs;
	double E_R_h_abs;
	double A_oleOOH_dec[3];					// oleOOH decomposition
	double n_oleOOH_dec[3];
	double E_oleOOH_dec[3];
	double A_etherOOH_dec[3];				// etherOOH decomposition
	double n_etherOOH_dec[3];
	double E_etherOOH_dec[3];
	double A_cyc_eth_dec;					// cyc-ether decomposition
	double n_cyc_eth_dec;
	double E_cyc_eth_dec;
	double A_allylic_rad_form[3];			// allylic radical formation
	double n_allylic_rad_form[3];
	double E_allylic_rad_form[3];
	double A_alkenyl_RO_form;				// alkenyl RO formation
	double n_alkenyl_RO_form;
	double E_alkenyl_RO_form;
	double A_alkenyl_RO_dec[4];				// alkenyl RO decomposition
	double n_alkenyl_RO_dec[4];
	double E_alkenyl_RO_dec[4]; 
	double A_ald_dec;						// aldehydes decomposition
	double n_ald_dec;
	double E_ald_dec;


	//double T;

public:

	double A;
	double n;
	double E;
	int sizeReactions = 0;
	Kinox(void);    // default constructor
	Kinox(std::string nome);    // constructor
	void leggi(char[80]);

	friend ostream& operator<<(ostream&, Kinox&);
	friend istream& operator>>(istream&, Kinox&);

	double v_initiation(Radicale r1, Radicale r2, int isomers, double Temp);
	double v_h_abstraction(HAbsRad r, Carbonio c, int numH, int isomers, std::string corr, double Temp);
	double v_isomerization_r(Radicale r, Idrogeno h, Anello a, int numH, double Temp);
	double v_beta_dec_r(Radicale r1, Radicale r2, double Temp);
	double v_ole_par_r(int numH, double Temp);
	double v_o2_add_r(Radicale r, double Temp);
	double v_o2_rem_roo(Radicale r, int numC, double Temp);
	double v_isom_roo(Radicale r, Idrogeno h, Anello a, int numH, double Temp);
	double v_isom_ooqooh(Radicale r, Idrogeno h, Anello a, int numH, double Temp);
	double v_isom_qooh(Radicale r1, Radicale r2, Anello a, int numC, double Temp);
	double v_isom_pooh2(Radicale r1, Radicale r2, Anello a, double Temp);
	double v_o2_add_qooh(Radicale r, double Temp);
	double v_o2_rem_ooqooh(Radicale r, int numC, double Temp);
	double v_ooqooh_to_khp(Radicale r1, Radicale r2, Anello a, int numH, double Temp);
	double v_khp_decomp(Radicale r, int dist, double Temp);
	double v_ole_par_roo(Radicale r1, Carbonio c, int numH, double Temp);
	double v_ole_par_ooqooh(Radicale r1, Carbonio c, int numH, double Temp);
	double v_ole_from_beta_qooh(Radicale r1, Radicale r2, double Temp);
	double v_ole_from_gamma_qooh(Radicale r1, Radicale r2, double Temp);
	double v_ole_from_delta_qooh(Radicale r1, Radicale r2, double Temp);
	double v_pooh2_dec_1(Radicale r1, Radicale r2, double Temp);
	double v_pooh2_dec_2(Radicale r1, Radicale r2, double Temp);
	double v_pooh2_dec_3(Radicale r1, Radicale r2, double Temp);
	double v_ether_from_qooh(Radicale r1, Radicale r2, AnelloO a, Carbonio corr, double Temp);
	double v_ether_from_pooh2(Radicale r1, Radicale r2, AnelloO a, double Temp);
	double v_beta_dec_qooh(Radicale r1, Radicale r2, double Temp);
	double v_OH_abs(Carbonio c, int numH, int isom, double Temp);
	double v_ooqooh_dec_tert(double Temp);
	double v_ooqooh_dec_quat(double Temp);
	double v_roo_h_abs(double Temp);
	double v_r_h_abs(double Temp);
	double v_oleooh_dec(Radicale r, double Temp);
	double v_etherooh_dec(Radicale r, double Temp);
	double v_cyc_eth_dec(double Temp);
	double v_allylic_rad_form(Carbonio c, int numH, double Temp);
	double v_alkenyl_ro_form(double Temp);
	double v_alkenyl_ro_dec(Carbonio c, double Temp);
	double v_ald_dec(double Temp);


	//void setT(double Temp);
	void wrireaLump(ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	void wrireaLump1(ofstream& stream, Molecola reag, HAbsRad rad, int numC, double A, double n, double E);
	void wrireaLump3(ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	void wrireaLump4(ofstream& stream, int numC, double A, double n, double E);     // R'numC' + O2 -> R'numC'OO
	void wrireaLump5(ofstream& stream, int numC, double A, double n, double E);     // R'numC'OO -> R'numC' + O2
	void wrireaLump6(ofstream& stream, int numC, double A, double n, double E);     // R'numC' + O2 -> OLE'numC' + HO2
	void wrireaLump7(ofstream& stream, int numC, double A, double n, double E);     // R'numC'OO -> Q'numC'OOH
	void wrireaLump8(ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH -> R'numC'OO
	void wrireaLump9(ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH -> ETER'numC' + OH
	void wrireaLump10(ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH -> OLE'numC' + HO2
	void wrireaLump11(ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	void wrireaLump11b(ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff1, double A1, double n1, double E1,
		std::vector<double> stoicCoeff2, double A2, double n2, double E2);
	void wrireaLump12(ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH + O2 -> OOQ'numC'OOH
	void wrireaLump13(ofstream& stream, int numC, double A, double n, double E);     // OOQ'numC'OOH -> Q'numC'OOH + O2
	void wrireaLump14(ofstream& stream, int numC, double A, double n, double E);     // OOQ'numC'OOH -> OQ'numC'OOH + OH
	void wrireaLump15(ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);

	std::string nameHAbsRad(HAbsRad rad);				// return a string containing the name of the radical rad
	std::string nameHAbsRadPlusH(HAbsRad rad);			// return a string containing the name of the radical rad when an H is added to it
	void wrirea(ofstream& stream, int tiporeaz, double A, double n, double E, Molecola m1, ...);

};




