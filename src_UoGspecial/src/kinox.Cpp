
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include "kinox.h"
using namespace std;

//#define FILE_READING_DEBUG
#define R 1.987204


istream& eol(istream& stream);

Kinox::Kinox(void)
{
	leggi((char*)"kinox.kin");
};

void readParameters(std::ifstream& input, int entries, RateRule &RR_)
{
	std::string cell;
	std::string line;
	std::string identifier;
	std::vector<std::string> dimension(RR_.dimensionality());
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
		for(int j = RR_.dimensionality(); j<4; j++)	//discard remaining cells
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

// 0 dimensional arrays
void readParameters(std::ifstream& input, double* A, double* n, double* E)
{
	std::string cell;
	std::string line;

	std::getline(input, line);
	std::stringstream linS(line);
	std::getline(linS, cell, ','); std::getline(linS, cell, ','); std::getline(linS, cell, ','); // discard the first 3 cells

	// get the indexes of the matrix of the values
	std::getline(linS, cell, ','); std::getline(linS, cell, ','); std::getline(linS, cell, ','); std::getline(linS, cell, ',');	// discard 4 cells

	std::getline(linS, cell, ',');		// get A
	std::istringstream Astr(cell);
	Astr >> *A;
	std::getline(linS, cell, ',');		// get n
	std::istringstream nstr(cell);
	nstr >> *n;
	std::getline(linS, cell, ',');		// get E
	std::istringstream Estr(cell);
	Estr >> *E;

#ifdef FILE_READING_DEBUG
	// DEBUG print the matrix
	std::cout << *A;
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << *n;
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << *E;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
#endif
}

// 1 dimensional arrays
template <size_t size>
void readParameters(std::ifstream& input, int entries, std::vector<std::string> names, double(&A_mat)[size], double(&n_mat)[size], double(&E_mat)[size])
{
	std::string cell;
	std::string line;
	int dim = names.size();

	for (int i = 0; i < entries; i++)
	{
		std::getline(input, line);
		std::stringstream linS(line);
		std::getline(linS, cell, ','); std::getline(linS, cell, ','); std::getline(linS, cell, ','); // discard the first 3 cells

		// get the indexes of the matrix of the values
		int ind1 = 0;
		std::getline(linS, cell, ',');
		bool foundInd1 = false;
		for (int j = 0; j < names.size(); j++)
		{
			if (cell == names[j])
			{
				ind1 = j;
				foundInd1 = true;
				break;
			}
		}
		if (foundInd1 == false) 
			std::cout << "ERROR: error in reading parameters, index name '" << cell << "' not valid!" << std::endl;

		std::getline(linS, cell, ','); std::getline(linS, cell, ','); std::getline(linS, cell, ',');	// discard 3 cells
		std::getline(linS, cell, ',');		// get A
		std::istringstream Astr(cell);
		Astr >> A_mat[ind1];
		std::getline(linS, cell, ',');		// get n
		std::istringstream nstr(cell);
		nstr >> n_mat[ind1];
		std::getline(linS, cell, ',');		// get E
		std::istringstream Estr(cell);
		Estr >> E_mat[ind1];
	}
#ifdef FILE_READING_DEBUG
	// DEBUG print the matrix
	for (int j = 0; j < dim; j++)
	{
		std::cout << A_mat[j] << "  ";
	}
	std::cout << std::endl;
	std::cout << std::endl;

	for (int j = 0; j < dim; j++)
	{
		std::cout << n_mat[j] << "  ";
	}
	std::cout << std::endl;

	std::cout << std::endl;

	for (int j = 0; j < dim; j++)
	{
		std::cout << E_mat[j] << "  ";
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
#endif
}

// 2 dimensional arrays
template <size_t size_x, size_t size_y>
void readParameters(std::ifstream& input, int entries, std::vector<std::string> namesRows, std::vector<std::string> namesColumns, double(&A_mat)[size_x][size_y], double(&n_mat)[size_x][size_y], double(&E_mat)[size_x][size_y], bool symmetric)
{
	std::string cell;
	std::string line;
	int rows = namesRows.size();
	int columns = namesColumns.size();

	// matrix to keep track of which values of the matrix are filled
	int filledValues[10][10];
	for (int i = 0; i < rows; i++) for (int j = 0; j < columns; j++) filledValues[i][j] = 0;

	for (int i = 0; i < entries; i++)
	{
		std::getline(input, line);
		std::stringstream linS(line);
		std::getline(linS, cell, ','); std::getline(linS, cell, ','); std::getline(linS, cell, ','); // discard the first 3 cells

		// get the indexes of the matrix of the values
		int ind1 = 0;
		int ind2 = 0;
		std::getline(linS, cell, ',');
		bool foundInd1 = false;
		for (int j = 0; j < namesRows.size(); j++)
		{
			if (cell == namesRows[j])
			{
				ind1 = j;
				foundInd1 = true;
				break;
			}
		}
		if (foundInd1 == false) 
			std::cout << "ERROR: error in reading parameters, index name '" << cell << "' not valid!" << std::endl;
		std::getline(linS, cell, ',');
		bool foundInd2 = false;
		for (int j = 0; j < namesColumns.size(); j++)
		{
			if (cell == namesColumns[j])
			{
				ind2 = j;
				foundInd2 = true;
				break;
			}
		}
		if (foundInd2 == false) 
			std::cout << "ERROR: error in reading kinetic parameters, index name '" << cell << "' not valid!" << std::endl;
		filledValues[ind1][ind2] = 1;

		std::getline(linS, cell, ','); std::getline(linS, cell, ',');	// discard 2 cells
		std::getline(linS, cell, ',');		// get A
		std::istringstream Astr(cell);
		Astr >> A_mat[ind1][ind2];
		std::getline(linS, cell, ',');		// get n
		std::istringstream nstr(cell);
		nstr >> n_mat[ind1][ind2];
		std::getline(linS, cell, ',');		// get E
		std::istringstream Estr(cell);
		Estr >> E_mat[ind1][ind2];
	}
	if (symmetric)
	{
		// fill the matrix
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < columns; j++)
			{
				if (filledValues[i][j] == 1)
				{
					A_mat[j][i] = A_mat[i][j];
					n_mat[j][i] = n_mat[i][j];
					E_mat[j][i] = E_mat[i][j];
				}
				else
				{
					A_mat[i][j] = A_mat[j][i];
					n_mat[i][j] = n_mat[j][i];
					E_mat[i][j] = E_mat[j][i];
				}
			}
		}
	}
#ifdef FILE_READING_DEBUG
	// DEBUG print the matrix
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			std::cout << A_mat[i][j] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			std::cout << n_mat[i][j] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			std::cout << E_mat[i][j] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
#endif
}

// 3 dimensional arrays
template <size_t size_1, size_t size_2, size_t size_3>
void readParameters(std::ifstream& input, int entries, std::vector<std::string> namesDim1, std::vector<std::string> namesDim2, std::vector<std::string> namesDim3, double(&A_mat)[size_1][size_2][size_3], double(&n_mat)[size_1][size_2][size_3], double(&E_mat)[size_1][size_2][size_3])
{
	std::string cell;
	std::string line;
	int dim1 = namesDim1.size();
	int dim2 = namesDim2.size();
	int dim3 = namesDim3.size();

	for (int i = 0; i < entries; i++)
	{
		std::getline(input, line);
		std::stringstream linS(line);
		std::getline(linS, cell, ','); std::getline(linS, cell, ','); std::getline(linS, cell, ','); // discard the first 3 cells

		// get the indexes of the matrix of the values
		int ind1 = 0;
		int ind2 = 0;
		int ind3 = 0;
		std::getline(linS, cell, ',');
		bool foundInd1 = false;
		for (int j = 0; j < namesDim1.size(); j++)
		{
			if (cell == namesDim1[j])
			{
				ind1 = j;
				foundInd1 = true;
				break;
			}
		}
		if (foundInd1 == false) 
			std::cout << "ERROR: error in reading parameters, index name '" << cell << "' not valid!" << std::endl;
		std::getline(linS, cell, ',');
		bool foundInd2 = false;
		for (int j = 0; j < namesDim2.size(); j++)
		{
			if (cell == namesDim2[j])
			{
				ind2 = j;
				foundInd2 = true;
				break;
			}
		}
		if (foundInd2 == false) 
			std::cout << "ERROR: error in reading kinetic parameters, index name '" << cell << "' not valid!" << std::endl;

		std::getline(linS, cell, ',');
		bool foundInd3 = false;
		for (int j = 0; j < namesDim3.size(); j++)
		{
			if (cell == namesDim3[j])
			{
				ind3 = j;
				foundInd3 = true;
				break;
			}
		}
		if (foundInd3 == false) 
			std::cout << "ERROR: error in reading kinetic parameters, index name '" << cell << "' not valid!" << std::endl;

		std::getline(linS, cell, ',');		// discard 1 cell
		std::getline(linS, cell, ',');		// get A
		std::istringstream Astr(cell);
		Astr >> A_mat[ind1][ind2][ind3];
		std::getline(linS, cell, ',');		// get n
		std::istringstream nstr(cell);
		nstr >> n_mat[ind1][ind2][ind3];
		std::getline(linS, cell, ',');		// get E
		std::istringstream Estr(cell);
		Estr >> E_mat[ind1][ind2][ind3];
	}
#ifdef FILE_READING_DEBUG
	// DEBUG print the matrix
	for (int i = 0; i < dim2; i++)
	{
		for (int k = 0; k < dim1; k++)
		{
			for (int j = 0; j < dim3; j++)
			{
				std::cout << setw(7) << A_mat[k][i][j] << "  ";
			}
			std::cout << "          ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	for (int i = 0; i < dim2; i++)
	{
		for (int k = 0; k < dim1; k++)
		{
			for (int j = 0; j < dim3; j++)
			{
				std::cout << setw(7) << n_mat[k][i][j] << "  ";
			}
			std::cout << "          ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	for (int i = 0; i < dim2; i++)
	{
		for (int k = 0; k < dim1; k++)
		{
			for (int j = 0; j < dim3; j++)
			{
				std::cout << setw(7) << E_mat[k][i][j] << "  ";
			}
			std::cout << "          ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
#endif
}

// 3 dimensional arrays
template <size_t size_1, size_t size_2, size_t size_3, size_t size_4>
void readParameters(std::ifstream& input, int entries, std::vector<std::string> namesDim1, std::vector<std::string> namesDim2, std::vector<std::string> namesDim3, std::vector<std::string> namesDim4, double(&A_mat)[size_1][size_2][size_3][size_4], double(&n_mat)[size_1][size_2][size_3][size_4], double(&E_mat)[size_1][size_2][size_3][size_4])
{
	std::string cell;
	std::string line;
	int dim1 = namesDim1.size();
	int dim2 = namesDim2.size();
	int dim3 = namesDim3.size();
	int dim4 = namesDim4.size();

	for (int i = 0; i < entries; i++)
	{
		std::getline(input, line);
		std::stringstream linS(line);
		std::getline(linS, cell, ','); std::getline(linS, cell, ','); std::getline(linS, cell, ','); // discard the first 3 cells

		// get the indexes of the matrix of the values
		int ind1 = 0;
		int ind2 = 0;
		int ind3 = 0;
		int ind4 = 0;
		std::getline(linS, cell, ',');
		bool foundInd1 = false;
		for (int j = 0; j < namesDim1.size(); j++)
		{
			if (cell == namesDim1[j])
			{
				ind1 = j;
				foundInd1 = true;
				break;
			}
		}
		if (foundInd1 == false) 
			std::cout << "ERROR: error in reading parameters, index name '" << cell << "' not valid!" << std::endl;
		std::getline(linS, cell, ',');
		bool foundInd2 = false;
		for (int j = 0; j < namesDim2.size(); j++)
		{
			if (cell == namesDim2[j])
			{
				ind2 = j;
				foundInd2 = true;
				break;
			}
		}
		if (foundInd2 == false) 
			std::cout << "ERROR: error in reading kinetic parameters, index name '" << cell << "' not valid!" << std::endl;

		std::getline(linS, cell, ',');
		bool foundInd3 = false;
		for (int j = 0; j < namesDim3.size(); j++)
		{
			if (cell == namesDim3[j])
			{
				ind3 = j;
				foundInd3 = true;
				break;
			}
		}
		if (foundInd3 == false) 
			std::cout << "ERROR: error in reading kinetic parameters, index name '" << cell << "' not valid!" << std::endl;

		std::getline(linS, cell, ',');
		bool foundInd4 = false;
		for (int j = 0; j < namesDim4.size(); j++)
		{
			if (cell == namesDim4[j])
			{
				ind4 = j;
				foundInd4 = true;
				break;
			}
		}
		if (foundInd4 == false) 
			std::cout << "ERROR: error in reading kinetic parameters, index name '" << cell << "' not valid!" << std::endl;

		//std::getline(linS, cell, ',');		// discard 1 cell
		std::getline(linS, cell, ',');		// get A
		std::istringstream Astr(cell);
		Astr >> A_mat[ind1][ind2][ind3][ind4];
		std::getline(linS, cell, ',');		// get n
		std::istringstream nstr(cell);
		nstr >> n_mat[ind1][ind2][ind3][ind4];
		std::getline(linS, cell, ',');		// get E
		std::istringstream Estr(cell);
		Estr >> E_mat[ind1][ind2][ind3][ind4];
	}
}

Kinox::Kinox(std::string nome)
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
			readParameters(input, 15, initiation);
			break;
		}
		case 2:
		{
			readParameters(input, 39, hAbstraction);
			break;
		}
		case 3:
		{
			readParameters(input, 27, isomerizationR);
			break;
		}
		case 4:
		{
			readParameters(input, 15, betaDecR);
			break;
		}
		case 5:
		{
			readParameters(input, 1, oleFromR);
			break;
		}
		case 6:
		{
			readParameters(input, 3, O2AdditionR);
			break;
		}
		case -6:
		{
			readParameters(input, 12, O2RemovalROO);
			break;
		}
		case 7:
		{
			readParameters(input, 36, isomROO);
			break;
		}
		case -7:
		{
			readParameters(input, 144, isomQOOH);
			break;
		}
		case 8:
		{
			readParameters(input, 3, O2AdditionQOOH);
			break;
		}
		case -8:
		{
			readParameters(input, 12, O2RemovalOOQOOH);
			break;
		}
		case 9:
		{
			readParameters(input, 24, OOQOOHToKHP);
			break;
		}
		case 10:
		{
			readParameters(input, 9, KHPDecomp);
			break;
		}
		case 11:
		{
			readParameters(input, 9, oleFromROO);
			break;
		}
		case 12:
		{
			readParameters(input, 59, etherFromQOOH);
			break;
		}
		case 13:
		{
			readParameters(input, 8, oleFromBetaQOOH);
			std::getline(input, line);
			readParameters(input, 9, oleFromGammaQOOH);
			std::getline(input, line);
			readParameters(input, 9, oleFromDeltaQOOH);
			break;
		}
		case 14:
		{
			readParameters(input, 36,isomOOQOOH);
			break;
		}
		case -14:
		{
			readParameters(input, 36, isomPOOH2);
			break;
		}
		case 15:
		{
			std::vector<std::string> namesRows{ "P", "S", "T" };
			std::vector<std::string> namesColumns{ "P", "S", "T" };
			readParameters(input, 9, oleFromOOQOOH);
			break;
		}
		case 16:
		{
			readParameters(input, 59, etherFromPOOH2);
			break;
		}
		case 17:
		{
			readParameters(input, 8, POOH2Dec1);
			std::getline(input, line);
			readParameters(input, 9, POOHDec2);
			std::getline(input, line);
			readParameters(input, 9, POOHDec3);
			break;
		}
		case 18:
		{
			readParameters(input, 3, oleOOHDec);
			break;
		}
		case 19:
		{
			readParameters(input, 3, etherOOHDec);
			break;
		}
		case 20:
		{
			readParameters(input, 1, cycEthDec);
			break;
		}
		case 21:
		{
			readParameters(input, 3, allylicRadForm);
			break;
		}
		case 22:
		{
			readParameters(input, 1, alkenylROForm);
			break;
		}
		case 23:
		{
			readParameters(input, 4, alkenylRODec);
			break;
		}
		case 24:
		{
			readParameters(input, 1, aldDec);
			break;
		}
		default:
		{
			break;
		}
		}
	}
};

void Kinox::leggi(char nome[80])
{
	ifstream myfile(nome);
	if (!myfile)
	{
		cerr << "Error in opening the kinetic constants file: " << nome;
		exit(1);
	};
	myfile >> *this;
	myfile.close();
};

istream& eol(istream& stream)
/////////////////////////////////////////////////////////////
// It is a stream manipulator:
//     stream >> eol      flushes until the end of the line
/////////////////////////////////////////////////////////////
{
	char ch;
	do { stream.get(ch); } while (ch != '\n');
	return stream;
};

ostream& operator<<(ostream& stream, Kinox& obj)
/////////////////////////////////////////////////////////////
// Write the object into a stream
/////////////////////////////////////////////////////////////
{
	char buff[40];

	return stream;
};

istream& operator>>(istream& stream, Kinox& obj)
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
	case ch3oo:
		return "CH3OO";
		break;
	default:
		std::cout << "Error in HAbs: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

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
		std::cout << "Error in PST: " << c << " is not a valid value!" << std::endl;
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
		std::cout << "Error in PST: " << r << " is not a valid value!" << std::endl;
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
		std::cout << "Error in PST: " << h << " is not a valid value!" << std::endl;
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
	case ao6:
		return "6 members";
		break;
	default:
		std::cout << "Error in ringO_comm: " << a << " is not a valid value!" << std::endl;
		break;
	}
}

reactionComment Kinox::v_initiation(Radicale r1, Radicale r2, int isomers)
{
	std::string RRID;
	std::vector<std::string> param = { MeEtPST(r1) , MeEtPST(r2) };
	std::vector<double> kinPar = initiation.returnRates(param, &RRID);
	A = kinPar[0] * double(isomers);
	n = kinPar[1];
	E = kinPar[2];

	std::string HRComm = "Type radical 1: ";
	HRComm.append(MeEtPST_comm(r1));
	HRComm.append("; Type radical 2: ");
	HRComm.append(MeEtPST_comm(r2));
	
	reactionComment comment(RRID, HRComm, 1, isomers);
	return comment;
}

reactionComment Kinox::v_h_abstraction(HAbsRad r, Carbonio c, int numH, int isomers, std::string corr)
{
	std::string RRID;
	std::vector<std::string> param = { HAbs(r) , PST(c) , corr};
	std::vector<double> kinPar = hAbstraction.returnRates(param, &RRID);
	A = kinPar[0] * double(isomers) * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	
	std::string comm = "Type of attacked carbon: ";
	comm.append(PST_comm(c));
	comm.append("; correction: ");
	comm.append(corr);
	
	reactionComment comment(RRID, comm, numH, isomers);
	return comment;
}

reactionComment Kinox::v_isomerization_r(Radicale r, Idrogeno h, Anello a, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) , PST(h) , ring(a) };
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

reactionComment Kinox::v_beta_dec_r(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), MeEtPST(r2)};
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

reactionComment Kinox::v_o2_add_r(Radicale r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	std::vector<double> kinPar = O2AdditionR.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));
	
	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_o2_rem_roo(Radicale r, int numC)
{
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));
	comm.append("; number of carbons: ");
	std::string numCStr;
	if (numC == 5)
	{
		numCStr = "5";
		comm.append("5");
	}
	if (numC > 5 && numC < 12)
	{
		numCStr = "6 to 11";
		comm.append("6 to 11");
	}
	if (numC > 11 && numC < 15)
	{
		numCStr = "12 to 14";
		comm.append("12 to 14");
	}
	if (numC > 14)
	{
		numCStr = "15+";
		comm.append("15+");
	}
	
	std::string RRID;
	std::vector<std::string> param = { PST(r), numCStr};
	std::vector<double> kinPar = O2RemovalROO.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID, comm);
	return comment;

}

reactionComment Kinox::v_isom_roo(Radicale r, Idrogeno h, Anello a, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r), PST(h), ring(a)};
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

reactionComment Kinox::v_isom_ooqooh(Radicale r, Idrogeno h, Anello a, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r), PST(h), ring(a) };
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

reactionComment Kinox::v_isom_qooh(Radicale r1, Radicale r2, Anello a, int numC)
{
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; C-H type: ");
	comm.append(PST_comm(r2));
	comm.append("; ring size : ");
	comm.append(ring_comm(a));
	comm.append("; number of carbons: ");
	std::string numCStr;
	if (numC == 5)
	{
		comm.append("5");
		numCStr = "5";
	}
	if (numC > 5 && numC < 12)
	{
		comm.append("6 to 11");
		numCStr = "6 to 11";
	}
	if (numC > 11 && numC < 15)
	{
		comm.append("12 to 14");
		numCStr = "12 to 14";
	}
	if (numC > 14)
	{
		comm.append("15+");
		numCStr = "15+";
	}
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2), ring(a), numCStr };
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
	std::vector<std::string> param = { PST(r1), PST(r2), ring(a)};
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

reactionComment Kinox::v_o2_add_qooh(Radicale r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r)};
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

reactionComment Kinox::v_o2_rem_ooqooh(Radicale r, int numC)
{
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));
	comm.append("; number of carbons: ");
	std::string numCStr;
	if (numC == 5)
	{
		numCStr = "5";
		comm.append("5");
	}
	if (numC > 5 && numC < 12)
	{
		numCStr = "6 to 11";
		comm.append("6 to 11");
	}
	if (numC > 11 && numC < 15)
	{
		numCStr = "12 to 14";
		comm.append("12 to 14");
	}
	if (numC > 14)
	{
		numCStr = "15+";
		comm.append("15+");
	}
	std::string RRID;
	std::vector<std::string> param = { PST(r) , numCStr };
	std::vector<double> kinPar = O2RemovalOOQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ooqooh_to_khp(Radicale r1, Radicale r2, Anello a, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) , ring(a)};
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

reactionComment Kinox::v_ole_par_roo(Radicale r1, Carbonio c, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(c)};
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

reactionComment Kinox::v_ole_par_ooqooh(Radicale r1, Carbonio c, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(c) };
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

reactionComment Kinox::v_ole_from_beta_qooh(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = oleFromBetaQOOH.returnRates(param, &RRID);
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

reactionComment Kinox::v_ole_from_gamma_qooh(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = oleFromGammaQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];;
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; radical type: ");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ole_from_delta_qooh(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = oleFromDeltaQOOH.returnRates(param, &RRID);
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

reactionComment Kinox::v_pooh2_dec_1(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
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

reactionComment Kinox::v_pooh2_dec_2(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = POOHDec2.returnRates(param, &RRID);
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

reactionComment Kinox::v_pooh2_dec_3(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = POOHDec3.returnRates(param, &RRID);
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

reactionComment Kinox::v_ether_from_qooh(Radicale r1, Radicale r2, AnelloO a, Carbonio corr)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2), ringO(a), PSTQ(corr)};
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
	comm.append("; carbon in ring: ");
	if (A == 0) // the required value with correction is not available, hence the one withoutis going to be used
	{
		std::cerr << "Warning: in v_ether_from_qooh() the value with the correction is not available" << std::endl;
		param = { PST(r1) , PST(r2), ringO(a), PSTQ(Cs) };
		kinPar = etherFromQOOH.returnRates(param, &RRID);
		A = kinPar[0];
		n = kinPar[1];
		E = kinPar[2];
		comm.append(PSTQ_comm(Cs));
	}
	else
	{
		comm.append(PSTQ_comm(corr));
	}

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ether_from_pooh2(Radicale r1, Radicale r2, AnelloO a, Carbonio corr)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2), ringO(a), PSTQ(corr) };
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
	comm.append(ringO_comm(a));
	if (A == 0) // the required value with correction is not available, hence the one withoutis going to be used
	{
		std::cerr << "Warning: in v_ether_from_qooh() the value with the correction is not available" << std::endl;
		param = { PST(r1) , PST(r2), ringO(a), PSTQ(Cs) };
		kinPar = etherFromPOOH2.returnRates(param, &RRID);
		A = kinPar[0];
		n = kinPar[1];
		E = kinPar[2];
		comm.append(PSTQ_comm(Cs));
	}
	else
	{
		comm.append(PSTQ_comm(corr));
	}
	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_oleooh_dec(Radicale r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r)};
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
	std::vector<double> kinPar = alkenylRODec.returnRates(param,&RRID);
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

//void Kinox::setT(double Temp)
//{
//	T = Temp;
//}

//void wrirea(ofstream &stream,int tiporeaz,double A,double E,Molecola m1 ,...);

void Kinox::wrireaLump(ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E)
{
	// write the lumped kinetic reaction in the form reag -> stoicCoeff[0] prod[1] + stoicCoeff[1] prod[2] + ...      A n E
	stream << reag << " -> ";
	stream << stoicCoeff[0] << " R" << prod[1].numberOfC();
	for (int i = 1; i < stoicCoeff.size(); i++)
	{
		stream << " + " << stoicCoeff[i] << " R" << prod[i + 1].numberOfC();
	}

	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << endl;
}

void Kinox::wrireaLump1(ofstream& stream, Molecola reag, HAbsRad rad, int numC, double A, double n, double E)
{
	// write the lumped kinetic reaction in the form reag + rad -> R_numC      A n E
	stream << reag << " + " << nameHAbsRad(rad) << " -> " << " R" << numC << " + " << nameHAbsRadPlusH(rad);
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
}

void Kinox::wrireaLump3(ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E)
{
	// write the lumped kinetic reaction in the form reag -> stoicCoeff[0] prod[1] + stoicCoeff[1] prod[2] + ...      A n E
	stream << "R" << reag.numberOfC() << " -> ";
	std::string prodName;
	if (prod[1].isole()) prodName = " OLE";
	else if (prod[1].isCRad()) prodName = " R";
	else
	{
		std::cerr << "Error: product of betadecomposition shold be either olefins or C* radical!" << std::endl;
		prodName = " ERR";
	}
	stream << stoicCoeff[0] << prodName << prod[1].numberOfC();

	for (int i = 1; i < stoicCoeff.size(); i++)
	{
		if (prod[i + 1].isole()) prodName = " OLE";
		else if (prod[i + 1].isCRad()) prodName = " R";
		else
		{
			std::cerr << "Error: product of betadecomposition shold be either olefins or C* radical!" << std::endl;
			prodName = " ERR";
		}
		stream << " + " << stoicCoeff[i] << prodName << prod[i + 1].numberOfC();
	}

	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << endl;
}

void Kinox::wrireaLump4(ofstream& stream, int numC, double A, double n, double E)     // RnumC + O2 -> RnumCOO
{
	stream << "R" << numC << " + " << "O2" << " -> " << " R" << numC << "OO";
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}

void Kinox::wrireaLump5(ofstream& stream, int numC, double A, double n, double E)     // RnumCOO -> RnumC + O2
{
	stream << "R" << numC << "OO" << " -> " << "R" << numC << " + " << " O2";
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}

void Kinox::wrireaLump6(ofstream& stream, int numC, double A, double n, double E)     // RnumC + O2 -> OLEnumC + HO2
{
	stream << "R" << numC << " + O2 -> " << "OLE" << numC << " + " << " HO2";
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}

void Kinox::wrireaLump7(ofstream& stream, int numC, double A, double n, double E)     // RnumCOO -> QnumCOOH
{
	stream << "R" << numC << "OO -> " << "Q" << numC << "OOH";
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}

void Kinox::wrireaLump8(ofstream& stream, int numC, double A, double n, double E)     // QnumCOOH -> RnumCOO
{
	stream << "Q" << numC << "OOH -> " << "R" << numC << "OO";
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}

void Kinox::wrireaLump9(ofstream& stream, int numC, double A, double n, double E)     // Q'numC'OOH -> ETER'numC' + OH
{
	stream << "Q" << numC << "OOH -> " << "ETER" << numC << " + OH";
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}

void Kinox::wrireaLump10(ofstream& stream, int numC, double A, double n, double E)     // Q'numC'OOH -> OLE'numC' + HO2
{
	stream << "Q" << numC << "OOH -> " << "OLE" << numC << " + HO2";
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}

void Kinox::wrireaLump11(ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E)
{
	stream << "Q" << reag.numberOfC() << "OOH -> 1.0 OH + ";
	if (prod[1].isole()) {
		stream << stoicCoeff[0] << " OLE" << prod[1].numberOfC();
	}
	else if (prod[1].ischeto())
	{
		stream << stoicCoeff[0] << " C" << prod[1].numberOfC() << "HO";
	}
	else
	{
		std::cerr << "Error: product of beta decomposition shold be either olefins or C* radical!" << std::endl;
		stream << stoicCoeff[0] << " ERR" << prod[1].numberOfC();
	}

	for (int i = 1; i < stoicCoeff.size(); i++)
	{
		if (prod[i + 1].isole())
		{
			stream << " + " << stoicCoeff[i] << " OLE" << prod[i + 1].numberOfC();
		}
		else if (prod[i + 1].ischeto())
		{
			stream << " + " << stoicCoeff[i] << " C" << prod[i + 1].numberOfC() << "HO";
		}
		else
		{
			std::cerr << "Error: product of beta decomposition shold be either olefins or C* radical!" << std::endl;
			stream << " + " << stoicCoeff[i] << " ERR" << prod[i + 1].numberOfC();
		}
	}

	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}

void Kinox::wrireaLump11b(ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff1, double A1, double n1, double E1,
	std::vector<double> stoicCoeff2, double A2, double n2, double E2)
{
	wrireaLump11(stream, reag, prod, stoicCoeff1, A1, n1, E1);
	stream << "DUP" << std::endl;
	wrireaLump11(stream, reag, prod, stoicCoeff2, A2, n2, E2);
	stream << "DUP" << std::endl;
}

void Kinox::wrireaLump12(ofstream& stream, int numC, double A, double n, double E)     // Q'numC'OOH + O2 -> OOQ'numC'OOH
{
	stream << "Q" << numC << "OOH + O2-> " << "OOQ" << numC << "OOH";
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}
void Kinox::wrireaLump13(ofstream& stream, int numC, double A, double n, double E)     // OOQ'numC'OOH -> Q'numC'OOH + O2
{
	stream << "OOQ" << numC << "OOH -> " << "Q" << numC << "OOH + O2";
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}

void Kinox::wrireaLump14(ofstream& stream, int numC, double A, double n, double E)     // OOQ'numC'OOH -> OQ'numC'OOH + OH
{
	stream << "OOQ" << numC << "OOH -> " << "OQ" << numC << "OOH + OH";
	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}

void Kinox::wrireaLump15(ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E)
{
	stream << "OQ" << reag.numberOfC() << "OOH -> 1.0 OH + ";
	if (prod[1].isRO()) {
		stream << stoicCoeff[0] << " R" << prod[1].numberOfC() << "O";
	}
	else if (prod[1].numCrad() == 0 && prod[1].numKeto() == 1)
	{
		stream << stoicCoeff[0] << " C" << prod[1].numberOfC() << "HO";
	}
	else
	{
		std::cerr << "Error: product of ketohydroperoxides decomposition shold be either CHO or RO radical!" << std::endl;
		stream << stoicCoeff[0] << " ERR" << prod[1].numberOfC();
	}

	for (int i = 1; i < stoicCoeff.size(); i++)
	{
		if (prod[i + 1].isRO())
		{
			stream << " + " << stoicCoeff[i] << " R" << prod[i + 1].numberOfC() << "O";
		}
		else if (prod[i + 1].numCrad() == 0 && prod[i + 1].numKeto() == 1)
		{
			stream << " + " << stoicCoeff[i] << " C" << prod[i + 1].numberOfC() << "HO";
		}
		else
		{
			std::cerr << "Error: product of ketohydroperoxides decomposition shold be either CHO or RO radical!" << std::endl;
			stream << " + " << stoicCoeff[i] << " ERR" << prod[i + 1].numberOfC();
		}
	}

	stream << "        A = " << A << ";     n = " << n << ";      E = " << E << endl;
	stream << std::endl;
}



std::string Kinox::nameHAbsRad(HAbsRad rad)
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

std::string Kinox::nameHAbsRadPlusH(HAbsRad rad)
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

// template <class Mol>
void Kinox::wrirea(ofstream& stream, int nreaz, double A, double n, double E, Molecola m1, ...)
/////////////////////////////////////////////////////////////////////////////
//  SCRIVE le reazioni sul file stream
//
//  1  A --> B
//  2  A --> B + C
//  3  A + O2 --> B
//  4  A --> B + O2
//  5  A + O2 --> B + HO2
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
// 21  CH3OO + A --> CH3OOH + B
// 22  HO2 + A --> B + OH
// 23  OH + A --> B + C + H2O
// 24  A --> B + C + D + E + OH 
// 25  OH + A --> B + C + D + H2O
// 26  OH + A --> B + CO + H2O
// 27  A --> B + C + D + E + F + OH 
//
/////////////////////////////////////////////////////////////////////////////
{
	char reaz[200];  // per contenere la reazione
	char cost[121];   // per contenere le costanti
	ostrstream outs(reaz, sizeof(reaz)); // prepara per scrivere la
										// stringa reaz come uno stream
	m1.output = FORMULA;
	//int size = m1.size() * 9;
	//if (sizeReactions > 97) sizeReactions = 97;

	va_list elenco;
	va_start(elenco, m1);  // prepara il puntatore elenco, che servira'
						  // per leggere gli argomenti variabili
	switch (nreaz)
	{
	case 3: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " + O2 -> " << m << '\0';
		break; };
	case 1: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m << '\0';
		break; };
	case 2: { Molecola m[2];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m[0] << " + " << m[1] << '\0';
		break; };
	case 4: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m << " + O2" << '\0';
		break; };
	case 5: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " + O2 -> " << m << " + HO2" << '\0';
		break; };
	case 6: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m << " + HO2" << '\0';
		break; };
	case 7: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m << " + OH" << '\0';
		break; };
	case 8: { Molecola m[2];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m[0] << " + " << m[1] << " + OH" << '\0';
		break; };
	case 9: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " + RH -> " << m << " + R*" << '\0';
		break; };
	case 10: { Molecola m[2];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m[0] << " + " << m[1] << " + 2 OH" << '\0';
		break; };
	case 11: { Molecola m[2];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m[0] << " + " << m[1] << " + HO2" << '\0';
		break; };
	case 12: { Molecola m[3];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		m[2] = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m[0] << " + " << m[1] << " + " << m[2] << " + OH" << '\0';
		break; };
	case 13: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " +  R* -> " << m << '\0';
		break; };
	case 14: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " +   O2 -> " << m << " +  HO2" << '\0';
		break; };
	case 15: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " +   OH -> " << m << " +  H2O" << '\0';
		break; };
	case 16: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " +    H -> " << m << " +   H2" << '\0';
		break; };
	case 17: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " +    O -> " << m << " +   OH" << '\0';
		break; };
	case 18: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " +  HO2 -> " << m << " + H2O2" << '\0';
		break; };
	case 19: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " +  CH3 -> " << m << " +  CH4" << '\0';
		break; };
	case 20: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " + C2H5 -> " << m << " + C2H6" << '\0';
		break; };
	case 21: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " + CH3OO -> " << m << " + CH3OOH" << '\0';
		break; };
	case 22: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " + HO2 -> " << m << " + OH" << '\0';
		break; };
	case 23: { Molecola m[2];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		outs << m1 << " + OH -> " << m[0] << " + " << m[1] << " + H2O" << '\0';
		break; };
	case 24: { Molecola m[4];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		m[2] = va_arg(elenco, Molecola);
		m[3] = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m[0] << " + " << m[1] << " + " << m[2] << " + " << m[3] << " + OH" << '\0';
		break; };
	case 25: { Molecola m[3];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		m[2] = va_arg(elenco, Molecola);
		outs << m1 << " + OH -> " << m[0] << " + " << m[1] << " + " << m[2] << " + H2O" << '\0';
		break; };
	case 26: { Molecola m = va_arg(elenco, Molecola);
		outs << m1 << " + OH -> " << m << " + CO + H2O" << '\0';
		break; };
	case 27: { Molecola m[5];
		m[0] = va_arg(elenco, Molecola);
		m[1] = va_arg(elenco, Molecola);
		m[2] = va_arg(elenco, Molecola);
		m[3] = va_arg(elenco, Molecola);
		m[4] = va_arg(elenco, Molecola);
		outs << m1 << " -> " << m[0] << " + " << m[1] << " + " << m[2] << " + " << m[3] << " + " << m[4] << " + OH" << '\0';
		break; };
	default:
		cerr << "Error in using wrirea! Not such case present in the options.";
		exit(200);
		break;
	};
	sprintf(cost, "/%8.2e/%8.4f/%8.0f/", A, n, E);
	stream << setw(sizeReactions) << reaz << "    " << cost << endl;
};

