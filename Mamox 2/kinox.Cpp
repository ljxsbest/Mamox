
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include "kinox.h"
//using namespace std;

//#define FILE_READING_DEBUG
#define R 1.987204


std::istream& eol(std::istream& stream);

Kinox::Kinox(void)
{
	leggi((char*)"kinox.kin");
};

void readParameters(std::ifstream& input, int entries, RateRule& RR_)
{
	std::string cell;
	std::string line;
	std::string identifier;
	std::vector<std::string> dimension(RR_.dimensionality());
	int maxDimensionality = 6; // maximum number of possible parameters for the rate rule.
	double A;
	double n;
	double E;
	for (int i = 0; i < entries; i++)
	{
		std::getline(input, line);
		std::stringstream linS(line);
		std::getline(linS, cell, ',');
		identifier = cell;
		std::getline(linS, cell, ','); std::getline(linS, cell, ','); // discard two cells

		for (int j = 0; j < RR_.dimensionality(); j++)
		{
			std::getline(linS, cell, ',');
			dimension[j] = cell;
		}
		for (int j = RR_.dimensionality(); j < maxDimensionality; j++)	//discard remaining cells
			std::getline(linS, cell, ',');

		std::getline(linS, cell, ',');		// get A
		A = std::stod(cell);
		std::getline(linS, cell, ',');		// get n
		n = std::stod(cell);
		std::getline(linS, cell, ',');		// get E
		E = std::stod(cell);

		RR_.addRate(identifier, dimension, A, n, E);
	}
}


void readParameters(std::ifstream& input, RateRule& RR_)
{

	std::streampos pos = input.tellg(); // 保存当前位置
	std::string l;
	int count = 0;

	while (std::getline(input, l)) {
		std::stringstream ss(l);
		std::string firstCell;
		if (!std::getline(ss, firstCell, ',')) break; // 没有第一列就结束
		// 去掉空格
		firstCell.erase(0, firstCell.find_first_not_of(" \t\r\n"));
		firstCell.erase(firstCell.find_last_not_of(" \t\r\n") + 1);

		if (firstCell.empty()) break; // 遇到空cell就停止
		count++;
	}

	input.clear(); // 清除 EOF 状态
	input.seekg(pos); // 回到计数前的位置
	int entries = count; // Get entries automatically instead of manually input.

	
	//std::cout << "Num of Class "<< ": " << entries << std::endl;

	std::string cell;
	std::string line;
	std::string identifier;
	std::vector<std::string> dimension(RR_.dimensionality());
	int maxDimensionality = 6; // maximum number of possible parameters for the rate rule.
	double A;
	double n;
	double E;

	for (int i = 0; i < entries; i++)
	{
		std::getline(input, line);
		std::stringstream linS(line);
		std::getline(linS, cell, ',');
		identifier = cell;
		std::getline(linS, cell, ','); std::getline(linS, cell, ','); // discard two cells

		for (int j = 0; j < RR_.dimensionality(); j++)
		{
			std::getline(linS, cell, ',');
			dimension[j] = cell;
		}
		for (int j = RR_.dimensionality(); j < maxDimensionality; j++)	//discard remaining cells
			std::getline(linS, cell, ',');

		std::getline(linS, cell, ',');		// get A
		A = std::stod(cell);
		std::getline(linS, cell, ',');		// get n
		n = std::stod(cell);
		std::getline(linS, cell, ',');		// get E
		E = std::stod(cell);

		RR_.addRate(identifier, dimension, A, n, E);
	}
}


//int Kinox::ClassCount(std::istream& stream, int index_class) {
//	std::streampos pos = stream.tellg(); // 保存当前位置
//	std::string line;
//	int count = 0;
//
//	while (std::getline(stream, line)) {
//		std::stringstream ss(line);
//		std::string firstCell;
//		if (!std::getline(ss, firstCell, ',')) break; // 没有第一列就结束
//		// 去掉空格
//		firstCell.erase(0, firstCell.find_first_not_of(" \t\r\n"));
//		firstCell.erase(firstCell.find_last_not_of(" \t\r\n") + 1);
//
//		if (firstCell.empty()) break; // 遇到空cell就停止
//		count++;
//	}
//
//	stream.clear(); // 清除 EOF 状态
//	stream.seekg(pos); // 回到计数前的位置
//	return count;
//}



Kinox::Kinox(std::string nome)  // TO DO: make the software read automatically the number of rate rules (so we don't have to define them by hand)
{
	//leggi(nome);
	std::ifstream input(nome);

	std::string line;
	std::getline(input, line);				// discard first line
	while (std::getline(input, line))		// get all the lines of the file one by one
	{
		std::stringstream lineStream(line);
		std::string index;
		std::getline(lineStream, index, ',');	// get the first cell

		if (index.size() == 0) index = "0";
		std::string cell;

		int value = std::stoi(index);
		switch (value)
		{
		case 1:					// read initiation reaction parameters
		{
			readParameters(input, initiation); // JIAXIN.
			break;
		}
		case 2:
		{
			readParameters(input, hAbstraction); // JIAXIN
			break;
		}

		case 3:
		{
			readParameters(input, isomerizationR); // JIAXIN.
			break;
		}
		case 4:
		{
			readParameters(input, betaDecR); // JIAXIN 15 to 16.
			break;
		}


		case 5:
		{
			readParameters(input, HAddOle);
			break;
		}

		case 6:
		{
			readParameters(input, OatomAddOle);
			break;
		}
		case 7:
		{
			readParameters(input, o2addvinyR);
			break;
		}
		case 8:
		{
			readParameters(input, HO2AddAllyR_rec);
			break;
		}
		case 9:
		{
			readParameters(input, HO2AddAllyR_dec);
			break;
		}
		case 10:
		{
			readParameters(input, rooh_to_ro);
			break;
		}
		case 11:
		{
			readParameters(input, ro_betadecom);
			break;
		}
		case 12:
		{
			readParameters(input, ch3AddAllyR);
			break;
		}
		case 13:
		{
			readParameters(input, ho2AddOle);
			break;
		}
		case 14:
		{
			readParameters(input, o2habsallyR);
			break;
		}
		case 15:
		{
			readParameters(input, OHAddOle); // JIAXIN OH addition to olefins.
			break;
		}
		case 16:
		{
			readParameters(input, ROHIsomReac); // JIAXIN betaROH isomerization.
			break;
		}
		case 17:
		{
			readParameters(input, O2AddROH); // JIAXIN O2 addition to ROH.
			break;
		}
		case 18:
		{
			readParameters(input, rohoo_wad);
			break;
		}
		case 19:
		{
			readParameters(input, wad_decom);
			break;
		}
		case 20:
		{
			readParameters(input, rohoo_elim);
			break;
		}
		case 21:
		{
			readParameters(input, ROHO2IsomReac); // JIAXIN betaROH isomerization.
			break;
		}
		case 22:
		{
			readParameters(input, bqohooh_decom);
			break;
		}
		case 23:
		{
			readParameters(input, qohoohtocEther);
			break;
		}
		case 24:
		{
			readParameters(input, o2addqohooh);
			break;
		}
		case 25:
		{
			readParameters(input, o2qohoohelim);
			break;
		}
		case 26:
		{
			readParameters(input, o2qohoohTokhp);
			break;
		}
		case 27:
		{
			readParameters(input, ohkhpdecom);
			break;
		}
		case 28:
		{
			readParameters(input, o2qohoohIsom);
			break;
		}
		case 29:
		{
			readParameters(input, ohpoohdecom);
			break;
		}
		case 30:
		{
			readParameters(input, ohpoohtoEther);
			break;
		}
		case 31:
		{
			readParameters(input, O2AdditionR);  // Jiaxin Change it to olefin format. 4 RRS totally in the rate rule file.
			break;
		}
		//case -6:
		//{
		//	readParameters(input, O2RemovalROO);
		//	break;
		//}
		case 32:
		{
			readParameters(input, oleFromROO); // Jiaxin. change to 7.
			break;
		}
		case 33:
		{
			readParameters(input, isomROO); // Jiaxin Change it from 27 to 29.
			break;
		}
		//case -7:
		//{
		//	readParameters(input, isomQOOH);
		//	break;
		//}
		case 34:
		{
			readParameters(input, oleqooh_decom);
			//std::getline(input, line);
			//readParameters(input, 9, oleFromGammaQOOH);
			//std::getline(input, line);
			//readParameters(input, 9, oleFromDeltaQOOH);
			break;
		}
		case 35:
		{
			readParameters(input, etherFromQOOH); // Jiaxin. 
			break;
		}
		case 36:
		{
			readParameters(input, oleOHsConsum);
			break;
		}
		case 37:
		{
			readParameters(input, O2AdditionQOOH); // Jiaxin Change it from 3 to 6.
			break;
		}
		//case -8:
		//{
		//	readParameters(input, O2RemovalOOQOOH);
		//	break;
		//}
		case 38:
		{
			readParameters(input, OOQOOHToKHP); // Jiaxin 15 to 45
			break;
		}
		case 39:
		{
			readParameters(input, KHPDecomp);
			break;
		}
		case 40:
		{
			readParameters(input, oleOOHDec);
			break;
		}
		case 41:
		{
			readParameters(input, isomOOQOOH);  // change it from 27 to 29.
			break;
		}
		//case -14:
		//{
		//	readParameters(input, isomPOOH2);
		//	break;
		//}
		case 42:
		{
			readParameters(input, POOH2Dec1);
			//std::getline(input, line);
			//readParameters(input, 9, POOHDec2);
			//std::getline(input, line);
			//readParameters(input, 9, POOHDec3);
			break;
		}
		case 43:
		{
			readParameters(input, etherFromPOOH2);
			break;
		}
		case 44:
		{
			std::vector<std::string> namesRows{ "P", "S", "T" };
			std::vector<std::string> namesColumns{ "P", "S", "T" };
			readParameters(input, oleFromOOQOOH);    // JIAXIN 
			break;
		}
		case 45:
		{
			readParameters(input, etherOOHDec);
			break;
		}
		case 46:
		{
			readParameters(input, 1, cycEthDec);
			break;
		}

		case 47:
		{
			readParameters(input, etherOOHDec);
			break;
		}
		case 48:
		{
			readParameters(input,1, cycEthDec1);
			break;
		}
		case 49:
		{
			readParameters(input, oleOOHDec);
			break;
		}
		case 50:
		{
			readParameters(input, Alpha_qohooh);
			break;
		}
		case 51:
		{
			readParameters(input, cycloadd);
			break;
		}
		case 52:
		{
			readParameters(input, cycoodec);
			break;
		}
		case 53:
		{
			readParameters(input, qohooh_decom);
			break;
		}

		case 54:
		{
			readParameters(input, AllyRrecom);
			break;
		}

		case 55:
		{
			readParameters(input, dieneDec);
			break;
		}
		
		
		//OHcEthOOH consumption
		// 51, 52: cycloaddition to be done.
		//case 21:
		//{
		//	readParameters(input, allylicRadForm);
		//	break;
		//}
		//case 22:
		//{
		//	readParameters(input, alkenylROForm);
		//	break;
		//}
		//case 23:
		//{
		//	readParameters(input, alkenylRODec);
		//	break;
		//}
		//case 24:
		//{
		//	readParameters(input, aldDec);
		//	break;
		//}




		default:
		{
			break;
		}
		}
	}
};

void Kinox::leggi(char nome[80])
{
	std::ifstream myfile(nome);
	if (!myfile)
	{
		std::cerr << "Error in opening the kinetic constants file: " << nome;
		exit(1);
	};
	myfile >> *this;
	myfile.close();
};

std::istream& eol(std::istream& stream)
/////////////////////////////////////////////////////////////
// It is a stream manipulator:
//     stream >> eol      flushes until the end of the line
/////////////////////////////////////////////////////////////
{
	char ch;
	do { stream.get(ch); } while (ch != '\n');
	return stream;
};

std::ostream& operator<<(std::ostream& stream, Kinox& obj)
/////////////////////////////////////////////////////////////
// Write the object into a stream
/////////////////////////////////////////////////////////////
{
	char buff[40];

	return stream;
};

std::istream& operator>>(std::istream& stream, Kinox& obj)
/////////////////////////////////////////////////////////////
// Read the object from a stream
/////////////////////////////////////////////////////////////
{
	char ch;

	return stream;
};

std::string MeEtPST(Radicale r)
{
	switch (r)
	{
	case Rpmet:
		return "me";
		break;
	case Rpet:
		return "et";
		break;
	case Rp:
		return "P";
		break;
	case Rs:
		return "S";
		break;
	case Rt:
		return "T";
		break;
	default:
		std::cout << "Error in MeEtPST: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

std::string MeEtPST_comm(Radicale r)
{
	switch (r)
	{
	case Rpmet:
		return "methyl";
		break;
	case Rpet:
		return "ethyl";
		break;
	case Rp:
		return "primary";
		break;
	case Rs:
		return "secondary";
		break;
	case Rt:
		return "tertiary";
		break;
	default:
		std::cout << "Error in MeEtPST_comm: " << r << " is not a valid value!" << std::endl;
		break;
	}
}


std::string HAbs(HAbsRad r)
{
	switch (r)
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
	case ch3o:
		return "CH3O";
		break;
	case ch3o2:
		return "CH3O2";
		break;
	default:
		std::cout << "Error in HAbs: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

//std::string HAbs(Molecola r)
//{
//	Molecola CH3;
//	CH3.makeCH3();
//	Molecola C2H5;
//	C2H5.makeC2H5();
//	
//	if (r.isSpecialMolecule() == 2)
//		return "O2";
//	else if (r.isSpecialMolecule() == 3)
//		return "O";
//	else if (r.isSpecialMolecule() == 4)
//		return "OH";
//	else if (r.isSpecialMolecule() == 5)
//		return "HO2";
//	else if (r.isSpecialMolecule() == 8)
//		return "H";
//	else if (r == CH3)
//		return "CH3";
//	else if (r == C2H5)
//		return "C2H5";
//	else
//		UTL::error("in HAbs: molecule has not been recognized!");
//	return "";
//}

std::string PST(Carbonio c)
{
	switch (c)
	{
	case Cp:
		return "P";
		break;
	case Cs:
		return "S";
		break;
	case Ct:
		return "T";
		break;
	default:
		std::cout << "Error in PST(C): " << c << " is not a valid value!" << std::endl;
		break;
	}
}

std::string n_all_vyn(int c)
{
	switch (c)
	{
	case 0:
		return "none";
		break;
	case 1:
		return "ally";
		break;
	case 2:
		return "viny";
		break;
	default:
		std::cout << "Error in n_all_vyn: " << c << " is not a valid value!" << std::endl;
		break;
	}
}

std::string n_all_vyn1(int c, std::string l)
{
	std::string type;
	switch (c)
	{
	case 0:
		type = "none";
		break;
	case 1:
		type = "ally";
		break;
	case 2:
		type = "viny";
		break;
	default:
		std::cout << "Error in n_all_vyn: " << c << " is not a valid value!" << std::endl;
		break;
	}
	if (l != "none")
	{
		type = type + l;
		std::cout << "type: " << type << std::endl;
	}
	return type;
}

std::string n_oh_ring(int c)
{
	switch (c)
	{
	case 0:
		return "none"; // ring size is 8.
		break;
	case 1:
		return "G1"; // outside the ring
		break;
	case 2:
		return "G2"; // on the ring
		break;
	default:
		std::cout << "Error in n_oh_ring: " << c << " is not a valid value!" << std::endl;
		break;
	}
}

std::string rr_sec_C(int c)
{
	switch (c)
	{
	case 0:
		return "none"; // here, 0/none means the C type is not secondary.
		break;
	case 12:
		return "12";
		break;
	case 11:
		return "11";
		break;
	case 10:
		return "10";
		break;
	case 9:
		return "9";
		break;
	case 8:
		return "8";
		break;
	case 7:
		return "7";
		break;
	case 6:
		return "6";
		break;
	case 5:
		return "5";
		break;
	case 4:
		return "4";
		break;
	default:
		std::cout << "Error in rr_sec_C: " << c << " the number exceeds range 5-12!" << std::endl;
		break;
	}
}


std::string PST_comm(Carbonio c)
{
	switch (c)
	{
	case Cp:
		return "primary";
		break;
	case Cs:
		return "secondary";
		break;
	case Ct:
		return "tertiary";
		break;
	default:
		std::cout << "Error in PST_comm: " << c << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PSTQ(Carbonio c)
{
	switch (c)
	{
	case Cp:
		return "P";
		break;
	case Cs:
		return "S";
		break;
	case Ct:
		return "T";
		break;
	case Cq:
		return "Q";
		break;
	default:
		std::cout << "Error in PSTQ: " << c << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PSTQ_comm(Carbonio c)
{
	switch (c)
	{
	case Cp:
		return "primary";
		break;
	case Cs:
		return "secondary";
		break;
	case Ct:
		return "tertiary";
		break;
	case Cq:
		return "quaternary";
		break;
	default:
		std::cout << "Error in PSTQ_comm: " << c << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PST(Radicale r)
{
	switch (r)
	{
	case Rp:
		return "P";
		break;
	case Rs:
		return "S";
		break;
	case Rt:
		return "T";
		break;
	default:
		std::cout << "Error in PST here: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PST_comm(Radicale r)
{
	switch (r)
	{
	case Rp:
		return "primary";
		break;
	case Rs:
		return "secondary";
		break;
	case Rt:
		return "tertiary";
		break;
	default:
		std::cout << "Error in PST_comm: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PST(Idrogeno h)
{
	switch (h)
	{
	case Hp:
		return "P";
		break;
	case Hs:
		return "S";
		break;
	case Ht:
		return "T";
		break;
	default:
		std::cout << "Error in PST(h): " << h << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PST_comm(Idrogeno h)
{
	switch (h)
	{
	case Hp:
		return "primary";
		break;
	case Hs:
		return "secondary";
		break;
	case Ht:
		return "tertiary";
		break;
	default:
		std::cout << "Error in PST_comm: " << h << " is not a valid value!" << std::endl;
		break;
	}
}

std::string ring(Anello a)
{
	switch (a)
	{
	case a5:
		return "5";
		break;
	case a6:
		return "6";
		break;
	case a7:
		return "7";
		break;
	case a8:
		return "8";
		break;
	default:
		std::cout << "Error in ring: " << a << " is not a valid value!" << std::endl;
		break;
	}
}


std::string ring_comm(Anello a)
{
	switch (a)
	{
	case a5:
		return "5 members";
		break;
	case a6:
		return "6 members";
		break;
	case a7:
		return "7 members";
		break;
	case a8:
		return "8 members";
		break;
	default:
		std::cout << "Error in ring_comm: " << a << " is not a valid value!" << std::endl;
		break;
	}
}

std::string ringO(AnelloO a)
{
	switch (a)
	{
	case ao3:
		return "3";
		break;
	case ao4:
		return "4";
		break;
	case ao5:
		return "5";
		break;
	case ao6:
		return "6";
		break;
	default:
		std::cout << "Error in ringO: " << a << " is not a valid value!" << std::endl;
		break;
	}
}

std::string ringO_comm(AnelloO a)
{
	switch (a)
	{
	case ao3:
		return "3 members";
		break;
	case ao4:
		return "4 members";
		break;
	case ao5:
		return "5 members";
		break;
	default:
		std::cout << "Error in ringO_comm: " << a << " is not a valid value!" << std::endl;
		break;
	}
}

std::string numCtolabel(int numC)
{
	if (numC == 5)
	{
		return "5";
	}
	if (numC > 5 && numC < 12)
	{
		return "6 to 11";
	}
	if (numC > 11 && numC < 15)
	{
		return "12 to 14";
	}
	if (numC > 14)
	{
		return "15+";
	}
	return "err";
}

reactionComment Kinox::v_initiation(Radicale r1, int isomers, int IsRAlly)
{
	std::string strR;
	if (MeEtPST(r1) == "me") strR = "me";
	else strR = "other";
	std::string RRID;
	if (IsRAlly == 2) strR = "--"; // for vinylic bond scission, the type of R does not matter.
	if (IsRAlly == 1) strR = "other"; // JL251007
	std::vector<std::string> param = { strR , n_all_vyn(IsRAlly)};
	std::vector<double> kinPar = initiation.returnRates(param, &RRID);
	A = kinPar[0] * double(isomers);
	n = kinPar[1];
	E = kinPar[2];
	std::string HRComm = "Type radical 1: ";
	HRComm.append(MeEtPST_comm(r1));
	//HRComm.append("; Type radical 2: ");
	//HRComm.append(MeEtPST_comm(r2));

	reactionComment comment(RRID, HRComm, 1, isomers);
	return comment;
}

reactionComment Kinox::v_initiation1(int numH_posr)
{
	std::string RRID;
	std::vector<std::string> param = { "H" , "ally" };
	std::vector<double> kinPar = initiation.returnRates(param, &RRID);
	A = kinPar[0] * double(numH_posr);
	n = kinPar[1];
	E = kinPar[2];

	std::string HRComm = "Type radical 1: H";

	reactionComment comment(RRID, HRComm, 1, numH_posr);
	return comment;
}

reactionComment Kinox::v_h_abstraction(Molecola r, Carbonio c, int isAllOrVyn, int numH, 
	int isomers, std::string OH_label)  // isAllOrVyn  =0 if not vyl nor all, =1 if all, =2 if vyn
{

	std::string RRID;
	std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn1(isAllOrVyn, OH_label)};
	std::vector<double> kinPar = hAbstraction.returnRates(param, &RRID);
	A = kinPar[0] * double(isomers) * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	// JIAXIN
	//std::cout << "radical = " << r << std::endl;
	//std::cout << "A = " << A << std::endl;

	std::string comm = "Type of attacked carbon: ";
	comm.append(PST_comm(c));
	comm.append("; correction: ");

	reactionComment comment(RRID, comm, numH, isomers);
	return comment;
}

reactionComment Kinox::v_dienes_dec(Molecola r, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { nameHAbsRad(r) };
	std::vector<double> kinPar = dieneDec.returnRates(param, &RRID);
	A = kinPar[0] * numH;
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "unity ";
	//comm.append(PST_comm(r));
	reactionComment comment(RRID, comm);
	return comment;
}


reactionComment Kinox::v_isomerization_r(Radicale r, Idrogeno h, Anello a, int numH, int Rally, int Hally)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) , PST(h) , ring(a), n_all_vyn(Rally), n_all_vyn(Hally)};
	std::vector<double> kinPar = isomerizationR.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));
	comm.append("; Type of hydrogen: ");
	comm.append(PST_comm(h));
	comm.append("; Size of ring: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm, numH);
	return comment;
}
// JIAXIN
reactionComment Kinox::v_roh_isom_r(Radicale r, Idrogeno h, Anello a, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) , PST(h) , ring(a)};
	std::vector<double> kinPar = ROHIsomReac.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));
	comm.append("; Type of hydrogen: ");
	comm.append(PST_comm(h));
	comm.append("; Size of ring: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm, numH);
	return comment;
}



reactionComment Kinox::v_beta_dec_r(Radicale r1, Radicale r2, int isAllyOrNot)
{

	std::string strR;
	if (MeEtPST(r2) == "me") strR = "me";
	else strR = "other";


	std::string RRID;
	std::vector<std::string> param;  // 先定义在外面
	if (n_all_vyn(isAllyOrNot) != "none")
		param = { "none", "none", n_all_vyn(isAllyOrNot)};
	else
		param = { PST(r1), strR, "none" };
	std::vector<double> kinPar = betaDecR.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of decomposing radical: ";
	comm.append(PST_comm(r1));
	comm.append("; Type of expelled radical: ");
	comm.append(MeEtPST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ole_par_r(int numH)
{
	std::string RRID;
	std::vector<std::string> param = { };
	std::vector<double> kinPar = oleFromR.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID, numH);
	return comment;
}

reactionComment Kinox::v_o2_add_r(Radicale r, int isAllyOrNot)
{
	//std::cout << "JL isAllyOrNot=" << isAllyOrNot << std::endl;

	std::string RRID;
	std::vector<std::string> param = { PST(r), n_all_vyn(isAllyOrNot) };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = O2AdditionR.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));

	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_o2_add_roh(Radicale r)
{
	//std::cout << "JL isAllyOrNot=" << isAllyOrNot << std::endl;

	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = O2AddROH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_alpha_qohooh(Radicale r)
{
	//std::cout << "JL isAllyOrNot=" << isAllyOrNot << std::endl;

	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = Alpha_qohooh.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));

	reactionComment comment(RRID, comm);
	return comment;
}

// JIAXIN 
reactionComment Kinox::v_oh_add_ole(Carbonio c1, Carbonio c2)
{
	//std::cout << "JL isAllyOrNot=" << isAllyOrNot << std::endl;

	std::string RRID;
	std::vector<std::string> param = { PST(c1), PST(c2) };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = OHAddOle.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	std::string comm = "OHAdd_site: ";
	comm.append(PST_comm(c1));
	comm.append("radical_site: ");
	comm.append(PST_comm(c2));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_cycloadd(int ring, std::string label)
{
	//std::cout << "JL isAllyOrNot=" << isAllyOrNot << std::endl;

	std::string RRID;
	std::vector<std::string> param = { std::to_string(ring), label };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = cycloadd.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	std::string comm = "cycloadd: ";
	//comm.append(PST_comm(ring));
	//comm.append("radical_site: ");
	//comm.append(PST_comm(label));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_h_add_ole(Carbonio c1, Carbonio c2) // JIAXIN OH add
{
	std::string RRID;
	std::vector<std::string> param = { PST(c1), PST(c2) };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = HAddOle.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	std::string comm = "HAdd_site: ";
	comm.append(PST_comm(c1));
	comm.append("radical_site: ");
	comm.append(PST_comm(c2));

	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_ho2_add_allyR_rec(Radicale r1) // JIAXIN
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = HO2AddAllyR_rec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	std::string comm = "radical_site: ";
	comm.append(PST_comm(r1));

	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_ho2_add_allyR_dec(Radicale r1) // JIAXIN
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = HO2AddAllyR_dec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "radical_site: ";
	comm.append(PST_comm(r1));
	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_ch3_add_allyR() // JIAXIN
{
	std::string RRID;
	std::vector<std::string> param = { "none"};
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = ch3AddAllyR.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "unity";
	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_allyR_recom() // JIAXIN
{
	std::string RRID;
	std::vector<std::string> param = { "none" };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = AllyRrecom.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "unity";
	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ho2_add_ole(Carbonio c1, Carbonio c2) // JIAXIN 
{
	std::string RRID;
	std::vector<std::string> param = { PST(c1), PST(c2) };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = ho2AddOle.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "HO2Add_site: ";
	comm.append(PST_comm(c1));
	comm.append("radical_site: ");
	comm.append(PST_comm(c2));
	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_o_add_ole(Carbonio c1, std::string s1) // JIAXIN 
{
	std::string RRID;
	std::vector<std::string> param = { PST(c1), s1 };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = OatomAddOle.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "OAdd_site: ";
	comm.append(PST_comm(c1));
	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_rooh_to_ro() // JIAXIN
{
	std::string RRID;
	std::vector<std::string> param = { "none" };
	std::vector<double> kinPar = rooh_to_ro.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "unity";
	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_ro_beta(std::string IsViny) // JIAXIN
{
	std::string RRID;
	std::vector<std::string> param = { IsViny };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = ro_betadecom.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "beta C type: ";
	comm.append(IsViny);
	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_o2_add_vinyR(Radicale r1) // JIAXIN
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = o2addvinyR.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "beta C type: ";
	comm.append(PST_comm(r1));
	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_o2_habs_allyR(Carbonio c1, int numH, int isomers)
{
	std::string RRID;
	std::vector<std::string> param = { PST(c1) };
	// std::vector<std::string> param = { nameHAbsRad(r) , PST(c), n_all_vyn(isAllOrVyn)};
	std::vector<double> kinPar = o2habsallyR.returnRates(param, &RRID);
	A = kinPar[0] * double(isomers) * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "C type: ";
	comm.append(PST_comm(c1));
	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_beta_rohoo_wad(Carbonio r1, Carbonio r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2) };
	std::vector<double> kinPar = rohoo_wad.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; type of C-OH:");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_wad_decom(Carbonio r1, Carbonio r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2) };
	std::vector<double> kinPar = wad_decom.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; type of CO*:");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_beta_rohoo_elim(Carbonio r1, Carbonio r2, int numH, std::string oh_to_r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2), oh_to_r };
	std::vector<double> kinPar = rohoo_elim.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; extracted H:");
	comm.append(PST_comm(r2));
	reactionComment comment(RRID, comm, numH);
	return comment;
}
reactionComment Kinox::v_beta_qohooh_decom(Radicale r1, Radicale r2, std::string oh_to_r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2), oh_to_r };

	std::vector<double> kinPar = bqohooh_decom.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; radical site:");
	comm.append(PST_comm(r2));
	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_qohooh_to_ether(Carbonio r1, Carbonio r2, AnelloO a, std::string oh_ring)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2), ringO(a), oh_ring };
	std::vector<double> kinPar = qohoohtocEther.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; type of radical:");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_qohooh_decom(std::string ooh_to_r)
{
	std::string RRID;
	std::vector<std::string> param = { ooh_to_r };
	std::vector<double> kinPar = qohooh_decom.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "ooh_to_r: ";
	comm.append(ooh_to_r);

	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_o2_add_qohooh(Radicale r, int numC)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r), rr_sec_C(numC) };
	std::vector<double> kinPar = o2addqohooh.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "radical type: ";
	comm.append(PST_comm(r));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_o2qohooh_elim(Carbonio r1, Carbonio r2, std::string oh_to_hsite, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2), oh_to_hsite};
	std::vector<double> kinPar = o2qohoohelim.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; expelled H:");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm, numH);
	return comment;
}
reactionComment Kinox::v_o2qohooh_to_khp(Carbonio r1, Carbonio r2, Anello ring_size, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2), ring(ring_size) };
	std::vector<double> kinPar = o2qohoohTokhp.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; Type of C-OOH:");
	comm.append(PST_comm(r2));
	comm.append("; size of ring: ");
	comm.append(ring_comm(ring_size));

	reactionComment comment(RRID, comm, numH);
	return comment;
}
reactionComment Kinox::v_ohkhp_decom(Radicale COOH_site, int dist)
{
	std::string posStr;
	if (dist == 1)
		posStr = "beta";
	if (dist == 2)
		posStr = "gamma";
	if (dist > 2)
		posStr = "other";

	std::string RRID;

	std::vector<std::string> param = { PST(COOH_site), posStr };
	std::vector<double> kinPar = ohkhpdecom.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OOH: ";
	comm.append(PST_comm(COOH_site));
	comm.append("; distance of keto_site:");
	comm.append(posStr);;

	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_o2qohooh_isom(Carbonio COO_site, Carbonio H_site, Anello ring_size, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(COO_site), PST(H_site), ring(ring_size) };
	std::vector<double> kinPar = o2qohoohIsom.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OO: ";
	comm.append(PST_comm(COO_site));
	comm.append("; type of H:");
	comm.append(PST_comm(H_site));
	comm.append("; size of ring: ");
	comm.append(ring_comm(ring_size));

	reactionComment comment(RRID, comm, numH);
	return comment;
}
reactionComment Kinox::v_ohpooh_decom(Carbonio COOH_site, Radicale radical_site, int dist)
{
	std::string posStr;
	if (dist == 1)
		posStr = "beta";
	if (dist == 2)
		posStr = "gamma";
	if (dist > 2)
		posStr = "delta";

	std::string RRID;
	std::vector<std::string> param = { PST(COOH_site), PST(radical_site), posStr };
	std::vector<double> kinPar = ohpoohdecom.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OOH: ";
	comm.append(PST_comm(COOH_site));

	reactionComment comment(RRID, comm);
	return comment;
}
reactionComment Kinox::v_ohpooh_cEther(Carbonio COOH_site, Radicale radical_site, AnelloO ring_size)
{
	std::string RRID;
	std::vector<std::string> param = { PST(COOH_site), PST(radical_site), ringO(ring_size) };
	std::vector<double> kinPar = ohpoohtoEther.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OOH: ";
	comm.append(PST_comm(COOH_site));
	comm.append("; type of radical:");
	comm.append(PST_comm(radical_site));
	//comm.append("; size of ring: ");
	//comm.append(ring_comm(ring_size));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_o2_rem_roo(Radicale r)
{
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));
	comm.append("; number of carbons: ");
	
	std::string RRID;
	std::vector<std::string> param = {PST(r)};
	std::vector<double> kinPar = O2RemovalROO.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID, comm);
	return comment;

}

reactionComment Kinox::v_isom_roo(Radicale r, Idrogeno h, Anello a, int numH, int Hally, int OOally) // Jiaxin, add 2 more conditions.
{
	std::string RRID;

	//std::cout << ring(a) << std::endl;
	//std::cout << n_all_vyn(OOally) << std::endl;
	//std::cout << n_all_vyn(Hally) << std::endl;

	std::vector<std::string> param = { PST(r), PST(h), ring(a) , n_all_vyn(OOally), n_all_vyn(Hally) };
	std::vector<double> kinPar = isomROO.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OO: ";
	comm.append(PST_comm(r));
	comm.append("; type of H:");
	comm.append(PST_comm(h));
	comm.append("; size of ring: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_isom_rohoo(Radicale r, Idrogeno h, Anello a, int numH, int OHonRing) // Jiaxin, add 2 more conditions.
{
	std::string RRID;
	std::vector<std::string> param = { PST(r), PST(h), ring(a) , n_oh_ring(OHonRing) };
	std::vector<double> kinPar = ROHO2IsomReac.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OO: ";
	comm.append(PST_comm(r));
	comm.append("; type of H:");
	comm.append(PST_comm(h));
	comm.append("; size of ring: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_isom_ooqooh(Radicale r, Idrogeno h, Anello a, int numH, int Hally, int OOally)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r), PST(h), ring(a), n_all_vyn(OOally),  n_all_vyn(Hally)};
	std::vector<double> kinPar = isomOOQOOH.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OO: ";
	comm.append(PST_comm(r));
	comm.append("; type of H:");
	comm.append(PST_comm(h));
	comm.append("; size of ring: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_isom_qooh(Radicale r1, Radicale r2, Anello a)
{
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; C-H type: ");
	comm.append(PST_comm(r2));
	comm.append("; ring size : ");
	comm.append(ring_comm(a));
	comm.append("; number of carbons: ");
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2), ring(a)};
	std::vector<double> kinPar = isomQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_isom_pooh2(Radicale r1, Radicale r2, Anello a)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2), ring(a) };
	std::vector<double> kinPar = isomPOOH2.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; type C-H: ");
	comm.append(PST_comm(r2));
	comm.append("; ring size: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_o2_add_qooh(Radicale r, int isAllyOrNot, int numC)
{
	std::string RRID, para1, para2;
	para1 = n_all_vyn(isAllyOrNot);
	para2 = rr_sec_C(numC);
	if (para1 == "ally") para2 = "none";
	std::vector<std::string> param = { PST(r),  para1, para2};

	
	//std::cout << "PST=" << PST(r) << std::endl;
	//std::cout << "n_all_vyn=" << n_all_vyn(isAllyOrNot) << std::endl;
	//std::cout << "rr_sec_C=" << rr_sec_C(numC) << std::endl;

	std::vector<double> kinPar = O2AdditionQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_o2_rem_ooqooh(Radicale r)
{
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));
	comm.append("; number of carbons: ");
	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	std::vector<double> kinPar = O2RemovalOOQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ooqooh_to_khp(Radicale r1, Radicale r2, Anello a, int numH, int roo_type, int rooh_type)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) , ring(a), n_all_vyn(roo_type), n_all_vyn(rooh_type)};
	std::vector<double> kinPar = OOQOOHToKHP.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; type C-OOH: ");
	comm.append(PST_comm(r2));
	comm.append("; ring size: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_khp_decomp(Radicale r, int dist)
{
	std::string posStr;
	if (dist == 1)
		posStr = "beta";
	if (dist == 2)
		posStr = "gamma";
	if (dist > 2)
		posStr = "other";

	std::string RRID;
	std::vector<std::string> param = { PST(r) , posStr };
	std::vector<double> kinPar = KHPDecomp.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r));
	comm.append("; CO position: ");
	switch (dist) {
	case 1:
		comm.append("beta");
		break;
	case 2:
		comm.append("gamma");
		break;
	default:
		comm.append("other");
		break;
	}

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ole_par_roo(Radicale r1, Carbonio c, int numH, int OOally, int Hally)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(c) , n_all_vyn(OOally), n_all_vyn(Hally)};
	std::vector<double> kinPar = oleFromROO.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; type C-H: ");
	comm.append(PST_comm(c));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_ole_par_ooqooh(Radicale r1, Carbonio c, int numH, int OOally, int Hally)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(c) , n_all_vyn(OOally), n_all_vyn(Hally)};
	std::vector<double> kinPar = oleFromOOQOOH.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; type C-H: ");
	comm.append(PST_comm(c));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_oleqooh_decom(Radicale r1, Radicale r2, int dist)
{
	std::string posStr;
	if (dist == 1)
		posStr = "beta";
	if (dist == 2)
		posStr = "gamma";
	if (dist == 3)
		posStr = "delta";

	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) , posStr};
	std::vector<double> kinPar = oleqooh_decom.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; radical type: ");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

//reactionComment Kinox::v_ole_from_gamma_qooh(Radicale r1, Radicale r2)
//{
//	std::string RRID;
//	std::vector<std::string> param = { PST(r1) , PST(r2) };
//	std::vector<double> kinPar = oleFromGammaQOOH.returnRates(param, &RRID);
//	A = kinPar[0];
//	n = kinPar[1];
//	E = kinPar[2];;
//	//return A * pow(Temp, n) * exp(-E / R / Temp);
//	std::string comm = "Type C-OOH: ";
//	comm.append(PST_comm(r1));
//	comm.append("; radical type: ");
//	comm.append(PST_comm(r2));
//
//	reactionComment comment(RRID, comm);
//	return comment;
//}
//reactionComment Kinox::v_ole_from_delta_qooh(Radicale r1, Radicale r2)
//{
//	std::string RRID;
//	std::vector<std::string> param = { PST(r1) , PST(r2) };
//	std::vector<double> kinPar = oleFromDeltaQOOH.returnRates(param, &RRID);
//	A = kinPar[0];
//	n = kinPar[1];
//	E = kinPar[2];
//	//return A * pow(Temp, n) * exp(-E / R / Temp);
//	std::string comm = "Type C-OOH: ";
//	comm.append(PST_comm(r1));
//	comm.append("; radical type: ");
//	comm.append(PST_comm(r2));
//
//	reactionComment comment(RRID, comm);
//	return comment;
//}

reactionComment Kinox::v_pooh2_decom(Radicale r1, Radicale r2, int dist)
{
	std::string posStr;
	if (dist == 1)
		posStr = "beta";
	if (dist == 2)
		posStr = "gamma";
	if (dist == 3)
		posStr = "delta";

	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2), posStr };
	std::vector<double> kinPar = POOH2Dec1.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	///return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; radical type: ");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

//reactionComment Kinox::v_pooh2_dec_2(Radicale r1, Radicale r2)
//{
//	std::string RRID;
//	std::vector<std::string> param = { PST(r1) , PST(r2) };
//	std::vector<double> kinPar = POOHDec2.returnRates(param, &RRID);
//	A = kinPar[0];
//	n = kinPar[1];
//	E = kinPar[2];
//	//return A * pow(Temp, n) * exp(-E / R / Temp);
//	std::string comm = "Type C-OOH: ";
//	comm.append(PST_comm(r1));
//	comm.append("; radical type: ");
//	comm.append(PST_comm(r2));
//
//	reactionComment comment(RRID, comm);
//	return comment;
//}
//reactionComment Kinox::v_pooh2_dec_3(Radicale r1, Radicale r2)
//{
//	std::string RRID;
//	std::vector<std::string> param = { PST(r1) , PST(r2) };
//	std::vector<double> kinPar = POOHDec3.returnRates(param, &RRID);
//	A = kinPar[0];
//	n = kinPar[1];
//	E = kinPar[2];
//	//return A * pow(Temp, n) * exp(-E / R / Temp);
//	std::string comm = "Type C-OOH: ";
//	comm.append(PST_comm(r1));
//	comm.append("; radical type: ");
//	comm.append(PST_comm(r2));
//
//	reactionComment comment(RRID, comm);
//	return comment;
//}

reactionComment Kinox::v_ether_from_qooh(Radicale r1, Radicale r2, AnelloO a, 
	std::string correction, int OOHally, int Rally)
{
	std::string RRID;

	//std::cout << PST(r1) << std::endl;
	//std::cout << PST(r2) << std::endl;
	//std::cout << ringO(a) << std::endl;
	//std::cout << n_all_vyn(OOHally) << std::endl;
	//std::cout << n_all_vyn(Rally) << std::endl;


	std::vector<std::string> param = { PST(r1) , PST(r2), ringO(a), correction, n_all_vyn(OOHally) , n_all_vyn(Rally)};
	std::vector<double> kinPar = etherFromQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; type of radical: ");
	comm.append(PST_comm(r2));
	comm.append("; size of ring: ");
	comm.append(ringO_comm(a));

	reactionComment comment(RRID, comm);
	return comment;
}


reactionComment Kinox::v_ether_from_pooh2(Carbonio r1, Radicale r2, AnelloO a)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2), ringO(a) };
	std::vector<double> kinPar = etherFromPOOH2.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; type of radical: ");
	comm.append(PST_comm(r2));
	comm.append("; ring size: ");
	
	reactionComment comment(RRID, comm);
	return comment;
}


reactionComment Kinox::v_oleooh_dec(Radicale r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	std::vector<double> kinPar = oleOOHDec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r));
	reactionComment comment(RRID, comm);
	return comment;
}


reactionComment Kinox::v_etherooh_dec(Radicale r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	std::vector<double> kinPar = etherOOHDec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_cyc_eth_dec()
{
	std::string RRID;
	std::vector<std::string> param = { };
	std::vector<double> kinPar = cycEthDec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	reactionComment comment(RRID);
	return comment;
}

reactionComment Kinox::v_cyc_eth_dec1() // OHcyceth
{
	std::string RRID;
	std::vector<std::string> param = { };
	std::vector<double> kinPar = cycEthDec1.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	reactionComment comment(RRID);
	return comment;
}

reactionComment Kinox::v_cycoodec(std::string label)
{
	std::string RRID;
	std::vector<std::string> param = { label };
	std::vector<double> kinPar = cycoodec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	reactionComment comment(RRID);
	return comment;
}


reactionComment Kinox::v_oleOHsConsum()
{
	std::string RRID;
	std::vector<std::string> param = { };
	std::vector<double> kinPar = oleOHsConsum.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	reactionComment comment(RRID);
	return comment;
}

reactionComment Kinox::v_allylic_rad_form(Carbonio c, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(c) };
	std::vector<double> kinPar = allylicRadForm.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(c));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_alkenyl_ro_form()
{
	std::string RRID;
	std::vector<std::string> param = { };
	std::vector<double> kinPar = alkenylROForm.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];


	reactionComment comment(RRID);
	return comment;
}

reactionComment Kinox::v_alkenyl_ro_dec(Carbonio c)
{
	std::string RRID;
	std::vector<std::string> param = { PSTQ(c) };
	std::vector<double> kinPar = alkenylRODec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of carbon: ";
	comm.append(PSTQ_comm(c));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ald_dec()
{
	std::string RRID;
	std::vector<std::string> param = { };
	std::vector<double> kinPar = aldDec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	reactionComment comment(RRID);
	return comment;
}


std::string Kinox::nameHAbsRad(Molecola r)
{

	Molecola CH3;
	CH3.makeCH3();
	Molecola C2H5;
	C2H5.makeC2H5();
	// JIAXIN
	//Molecola CH3O;
	//CH3O.makeCH3O();
	//Molecola CH3O2;
	//CH3O2.makeCH3O2();

	// JIAXIN
	if (r.isSpecialMolecule() == 8)
		return "H";
	else if (r.isSpecialMolecule() == 4)
		return "OH";
	else if (r.isSpecialMolecule() == 3)
		return "O";
	else if (r.isSpecialMolecule() == 2)
		return "O2";
	else if (r.isSpecialMolecule() == 5)
		return "HO2";
	else if (r == CH3)
		return "CH3";
	else if (r.isSpecialMolecule() == 12)
		return "CH3O";
	else if (r.isSpecialMolecule() == 14)
		return "CH3O2";

	else if (r == C2H5)
		return "C2H5";
	else
		UTL::error("in HAbs: molecule has not been recognized!");
	return "ERR.";
}

std::string Kinox::nameHAbsRadPlusH(Molecola r)
{
	Molecola CH3;
	CH3.makeCH3();
	Molecola C2H5;
	C2H5.makeC2H5();

	//Molecola CH3O;
	//CH3O.makeCH3O();
	//Molecola CH3O2;
	//CH3O2.makeCH3O2();


	if (r.isSpecialMolecule() == 2)
		return "HO2";
	else if (r.isSpecialMolecule() == 3)
		return "OH";
	else if (r.isSpecialMolecule() == 4)
		return "H2O";
	else if (r.isSpecialMolecule() == 5)
		return "H2O2";
	else if (r.isSpecialMolecule() == 8)
		return "H2";
	else if (r == CH3)
		return "CH4";

	else if (r.isSpecialMolecule() == 12)
		return "CH3OH";
	else if (r.isSpecialMolecule() == 14)
		return "CH3O2H";

	else if (r == C2H5)
		return "C2H6";
	else
		UTL::error("in HAbs: molecule has not been recognized!");
	return "ERR.";

}

