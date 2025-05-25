
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
			std::vector<std::string> namesRows{ "me","et","P","S","T" };
			std::vector<std::string> namesColumns{ "me","et","P","S","T" };
			readParameters(input, 15, namesRows, namesColumns, A_initiation, n_initiation, E_initiation, true);
			break;
		}
		case 2:
		{
			std::vector<std::string> namesDim1{ "O2","OH","H","O","HO2","CH3","C2H5", "CH3O2"};
			std::vector<std::string> namesDim2{ "P","S","T" };
			std::vector<std::string> namesDim3{ "none","_1","_2","_3","_01","_02","_03","_11","_11'","_12","_23","_000","_001","_002","_011"};
			readParameters(input, 37, namesDim1, namesDim2, namesDim3, A_h_abstraction, n_h_abstraction, E_h_abstraction);
			break;
		}
		case 3:
		{
			std::vector<std::string> namesDim1{ "P", "S", "T" };
			std::vector<std::string> namesDim2{ "P", "S", "T" };
			std::vector<std::string> namesDim3{ "5", "6", "7" };
			readParameters(input, 27, namesDim1, namesDim2, namesDim3, A_isometization_R, n_isometization_R, E_isometization_R);
			break;
		}
		case 4:
		{
			std::vector<std::string> namesRows{ "P", "S", "T" };
			std::vector<std::string> namesColumns{ "me", "et", "P", "S", "T" };
			readParameters(input, 15, namesRows, namesColumns, A_beta_dec_R, n_beta_dec_R, E_beta_dec_R, false);
			break;
		}
		case 5:
		{
			readParameters(input, &A_ole_from_R, &n_ole_from_R, &E_ole_from_R);
			break;
		}
		case 6:
		{
			std::vector<std::string> names{ "P", "S", "T" };
			readParameters(input, 3, names, A_O2_addition_R, n_O2_addition_R, E_O2_addition_R);
			break;
		}
		case -6:
		{
			std::vector<std::string> names{ "P", "S", "T" };
			std::vector<std::string> names2{ "5", "6 to 11", "12 to 14", "15+" };
			readParameters(input, 12, names, names2, A_O2_removal_ROO, n_O2_removal_ROO, E_O2_removal_ROO, false);
			break;
		}
		case 7:
		{
			std::vector<std::string> namesDim1{ "P", "S", "T" };
			std::vector<std::string> namesDim2{ "P", "S", "T" };
			std::vector<std::string> namesDim3{ "5", "6", "7", "8"};
			readParameters(input, 36, namesDim1, namesDim2, namesDim3, A_isom_ROO, n_isom_ROO, E_isom_ROO);
			break;
		}
		case -7:
		{
			std::vector<std::string> namesDim1{ "P", "S", "T" };
			std::vector<std::string> namesDim2{ "P", "S", "T" };
			std::vector<std::string> namesDim3{ "5", "6", "7" };
			std::vector<std::string> namesDim4{ "5", "6 to 11", "12 to 14", "15+" };
			readParameters(input, 108, namesDim1, namesDim2, namesDim3, namesDim4, A_isom_QOOH, n_isom_QOOH, E_isom_QOOH);
			break;
		}
		case 8:
		{
			std::vector<std::string> names{ "P", "S", "T" };
			readParameters(input, 3, names, A_O2_addition_QOOH, n_O2_addition_QOOH, E_O2_addition_QOOH);
			break;
		}
		case -8:
		{
			std::vector<std::string> names{ "P", "S", "T" };
			std::vector<std::string> names2{ "5", "6 to 11", "12 to 14", "15+" };
			readParameters(input, 12, names, names2, A_O2_removal_OOQOOH, n_O2_removal_OOQOOH, E_O2_removal_OOQOOH, false);
			break;
		}
		case 9:
		{
			std::vector<std::string> namesDim1{ "P", "S", "T" };
			std::vector<std::string> namesDim2{ "P", "S", "T" };
			std::vector<std::string> namesDim3{ "5", "6", "7", "8"};
			readParameters(input, 24, namesDim1, namesDim2, namesDim3, A_OOQOOH_to_KHP, n_OOQOOH_to_KHP, E_OOQOOH_to_KHP);
			break;
		}
		case 10:
		{
			std::vector<std::string> namesDim1{ "P", "S", "T" };
			std::vector<std::string> namesDim2{ "beta", "gamma", "other" };
			readParameters(input, 9, namesDim1, namesDim2, A_KHP_decomp, n_KHP_decomp, E_KHP_decomp, false);
			break;
		}
		case 11:
		{
			std::vector<std::string> namesRows{ "P", "S", "T" };
			std::vector<std::string> namesColumns{ "P", "S", "T" };
			readParameters(input, 9, namesRows, namesColumns, A_ole_from_ROO, n_ole_from_ROO, E_ole_from_ROO, false);
			break;
		}
		case 12:
		{
			std::vector<std::string> namesDim1{ "P", "S", "T" };
			std::vector<std::string> namesDim2{ "P", "S", "T" };
			std::vector<std::string> namesDim3{ "3", "4", "5", "6"};
			std::vector<std::string> namesDim4{ "s", "t", "q" };
			readParameters(input, 59, namesDim1, namesDim2, namesDim3, namesDim4, A_ether_from_QOOH, n_ether_from_QOOH, E_ether_from_QOOH);
			break;
		}
		case 13:
		{
			std::vector<std::string> namesRows{ "P", "S", "T" };
			std::vector<std::string> namesColumns{ "P", "S", "T" };
			readParameters(input, 8, namesRows, namesColumns, A_ole_from_beta_QOOH, n_ole_from_beta_QOOH, E_ole_from_beta_QOOH, false);
			std::getline(input, line);
			readParameters(input, 9, namesRows, namesColumns, A_ole_from_gamma_QOOH, n_ole_from_gamma_QOOH, E_ole_from_gamma_QOOH, false);
			std::getline(input, line);
			readParameters(input, 9, namesRows, namesColumns, A_ole_from_delta_QOOH, n_ole_from_delta_QOOH, E_ole_from_delta_QOOH, false);
			break;
		}
		case 14:
		{
			std::vector<std::string> namesDim1{ "P", "S", "T" };
			std::vector<std::string> namesDim2{ "P", "S", "T" };
			std::vector<std::string> namesDim3{ "5", "6", "7", "8"};
			readParameters(input, 36, namesDim1, namesDim2, namesDim3, A_isom_OOQOOH, n_isom_OOQOOH, E_isom_OOQOOH);
			break;
		}
		case -14:
		{
			std::vector<std::string> namesDim1{ "P", "S", "T" };
			std::vector<std::string> namesDim2{ "P", "S", "T" };
			std::vector<std::string> namesDim3{ "5", "6", "7" };
			readParameters(input, 27, namesDim1, namesDim2, namesDim3, A_isom_POOH2, n_isom_POOH2, E_isom_POOH2);
			break;
		}
		case 15:
		{
			std::vector<std::string> namesRows{ "P", "S", "T" };
			std::vector<std::string> namesColumns{ "P", "S", "T" };
			readParameters(input, 9, namesRows, namesColumns, A_ole_from_OOQOOH, n_ole_from_OOQOOH, E_ole_from_OOQOOH, false);
			break;
		}
		case 16:
		{
			std::vector<std::string> namesDim1{ "P", "S", "T" };
			std::vector<std::string> namesDim2{ "P", "S", "T" };
			std::vector<std::string> namesDim3{ "3", "4", "5", "6"};
			readParameters(input, 36, namesDim1, namesDim2, namesDim3, A_ether_from_POOH2, n_ether_from_POOH2, E_ether_from_POOH2);
			break;
		}
		case 17:
		{
			std::vector<std::string> namesRows{ "P", "S", "T" };
			std::vector<std::string> namesColumns{ "P", "S", "T" };
			readParameters(input, 8, namesRows, namesColumns, A_POOH2_dec_1, n_POOH2_dec_1, E_POOH2_dec_1, false);
			std::getline(input, line);
			readParameters(input, 9, namesRows, namesColumns, A_POOH2_dec_2, n_POOH2_dec_2, E_POOH2_dec_2, false);
			std::getline(input, line);
			readParameters(input, 9, namesRows, namesColumns, A_POOH2_dec_3, n_POOH2_dec_3, E_POOH2_dec_3, false);
			break;
		}
		case 18:
		{
			std::vector<std::string> names{ "P", "S", "T" };
			readParameters(input, 3, names, A_oleOOH_dec, n_oleOOH_dec, E_oleOOH_dec);
			break;
		}
		case 19:
		{
			std::vector<std::string> names{ "P", "S", "T" };
			readParameters(input, 3, names, A_etherOOH_dec, n_etherOOH_dec, E_etherOOH_dec);
			break;
		}
		case 20:
		{
			readParameters(input, &A_cyc_eth_dec, &n_cyc_eth_dec, &E_cyc_eth_dec);
			break;
		}
		case 21:
		{
			std::vector<std::string> names{ "P", "S", "T" };
			readParameters(input, 3, names, A_allylic_rad_form, n_allylic_rad_form, E_allylic_rad_form);
			break;
		}
		case 22:
		{
			readParameters(input, &A_alkenyl_RO_form, &n_alkenyl_RO_form, &E_alkenyl_RO_form);
			break;
		}
		case 23:
		{
			std::vector<std::string> names{ "P", "S", "T", "Q" };
			readParameters(input, 4, names, A_alkenyl_RO_dec, n_alkenyl_RO_dec, E_alkenyl_RO_dec);
			break;
		}
		case 24:
		{
			readParameters(input, &A_ald_dec, &n_ald_dec, &E_ald_dec);
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

int MeEtPST(Radicale r)
{
	switch (r)
	{
	case Rpmet:
		return 0;
		break;
	case Rpet:
		return 1;
		break;
	case Rp:
		return 2;
		break;
	case Rs:
		return 3;
		break;
	case Rt:
		return 4;
		break;
	default:
		std::cout << "Error in MeEtPST: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

int HAbs(HAbsRad r)
{
	switch (r)
	{
	case o2:
		return 0;
		break;
	case oh:
		return 1;
		break;
	case h:
		return 2;
		break;
	case o:
		return 3;
		break;
	case ho2:
		return 4;
		break;
	case ch3:
		return 5;
		break;
	case c2h5:
		return 6;
		break;
	case ch3oo:
		return 7;
		break;
	default:
		std::cout << "Error in HAbs: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

int PST(Carbonio c)
{
	switch (c)
	{
	case Cp:
		return 0;
		break;
	case Cs:
		return 1;
		break;
	case Ct:
		return 2;
		break;
	default:
		std::cout << "Error in PST: " << c << " is not a valid value!" << std::endl;
		break;
	}
}

int PSTQ(Carbonio c)
{
	switch (c)
	{
	case Cp:
		return 0;
		break;
	case Cs:
		return 1;
		break;
	case Ct:
		return 2;
		break;
	case Cq:
		return 3;
		break;
	default:
		std::cout << "Error in PSTQ: " << c << " is not a valid value!" << std::endl;
		break;
	}
}

int PST(Radicale r)
{
	switch (r)
	{
	case Rp:
		return 0;
		break;
	case Rs:
		return 1;
		break;
	case Rt:
		return 2;
		break;
	default:
		std::cout << "Error in PST: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

int PST(Idrogeno h)
{
	switch (h)
	{
	case Hp:
		return 0;
		break;
	case Hs:
		return 1;
		break;
	case Ht:
		return 2;
		break;
	default:
		std::cout << "Error in PST: " << h << " is not a valid value!" << std::endl;
		break;
	}
}

int ring(Anello a)
{
	switch (a)
	{
	case a5:
		return 0;
		break;
	case a6:
		return 1;
		break;
	case a7:
		return 2;
		break;
	case a8:
		return 3;
		break;
	default:
		std::cout << "Error in ring: " << a << " is not a valid value!" << std::endl;
		break;
	}
}

int ringO(AnelloO a)
{
	switch (a)
	{
	case ao3:
		return 0;
		break;
	case ao4:
		return 1;
		break;
	case ao5:
		return 2;
		break;
	case ao6:
		return 3;
		break;
	default:
		std::cout << "Error in ringO: " << a << " is not a valid value!" << std::endl;
		break;
	}
}

double Kinox::v_initiation(Radicale r1, Radicale r2, int isomers, double Temp)
{
	//enum Radicale { Rpmet, Rpet, Rp, Rs, Rt };
	int i = MeEtPST(r1);
	int j = MeEtPST(r2);
	A = A_initiation[i][j] * double(isomers);
	n = n_initiation[i][j];
	E = E_initiation[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_h_abstraction(HAbsRad r, Carbonio c, int numH, int isomers, std::string corr, double Temp)
{
	int k = 0;
	if (corr == "none")
		k = 0;
	else if (corr == "_1")
		k = 1;
	else if (corr == "_2")
		k = 2;
	else if (corr == "_3")
		k = 3;
	else if (corr == "_01")
		k = 4;
	else if (corr == "_02")
		k = 5;
	else if (corr == "_03")
		k = 6;
	else if (corr == "_11")
		k = 7;
	else if (corr == "_11'")
		k = 8;
	else if (corr == "_12")
		k = 9;
	else if (corr == "_23")
		k = 10;
	else if (corr == "_000")
		k = 11;
	else if (corr == "_001")
		k = 12;
	else if (corr == "_002")
		k = 13;
	else if (corr == "_011")
		k = 14;

	int i = HAbs(r);
	int j = PST(c);
	A = A_h_abstraction[i][j][k] * double(isomers) * double(numH);
	n = n_h_abstraction[i][j][k];
	E = E_h_abstraction[i][j][k];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_isomerization_r(Radicale r, Idrogeno h, Anello a, int numH, double Temp)
{
	int i = PST(r);
	int j = PST(h);
	int k = ring(a);
	A = A_isometization_R[i][j][k] * double(numH);
	n = n_isometization_R[i][j][k];
	E = E_isometization_R[i][j][k];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_beta_dec_r(Radicale r1, Radicale r2, double Temp)
{
	int i = PST(r1);
	int j = MeEtPST(r2);
	A = A_beta_dec_R[i][j];
	n = n_beta_dec_R[i][j];
	E = E_beta_dec_R[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ole_par_r(int numH, double Temp)
{
	A = A_ole_from_R * double(numH);
	n = n_ole_from_R;
	E = E_ole_from_R;
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_o2_add_r(Radicale r, double Temp)
{
	int i = PST(r);
	A = A_O2_addition_R[i];
	n = n_O2_addition_R[i];
	E = E_O2_addition_R[i];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_o2_rem_roo(Radicale r, int numC, double Temp)
{
	int i = PST(r);
	int j = 0;
	if (numC == 5)
		j = 0;
	if (numC > 5 && numC < 12)
		j = 1;
	if (numC > 11 && numC < 15)
		j = 2;
	if (numC > 14)
		j = 3;
	A = A_O2_removal_ROO[i][j];
	n = n_O2_removal_ROO[i][j];
	E = E_O2_removal_ROO[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_isom_roo(Radicale r, Idrogeno h, Anello a, int numH, double Temp)
{
	int i = PST(r);
	int j = PST(h);
	int k = ring(a);
	A = A_isom_ROO[i][j][k] * double(numH);
	n = n_isom_ROO[i][j][k];
	E = E_isom_ROO[i][j][k];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_isom_ooqooh(Radicale r, Idrogeno h, Anello a, int numH, double Temp)
{
	int i = PST(r);
	int j = PST(h);
	int k = ring(a);
	A = A_isom_OOQOOH[i][j][k] * double(numH);
	n = n_isom_OOQOOH[i][j][k];
	E = E_isom_OOQOOH[i][j][k];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_isom_qooh(Radicale r1, Radicale r2, Anello a, int numC, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	int k = ring(a);
	int h = 0;
	if (numC == 5)
		h = 0;
	if (numC > 5 && numC < 12)
		h = 1;
	if (numC > 11 && numC < 15)
		h = 2;
	if (numC > 14)
		h = 3;
	A = A_isom_QOOH[i][j][k][h];
	n = n_isom_QOOH[i][j][k][h];
	E = E_isom_QOOH[i][j][k][h];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_isom_pooh2(Radicale r1, Radicale r2, Anello a, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	int k = ring(a);
	A = A_isom_POOH2[i][j][k];
	n = n_isom_POOH2[i][j][k];
	E = E_isom_POOH2[i][j][k];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_o2_add_qooh(Radicale r, double Temp)
{
	int i = PST(r);
	A = A_O2_addition_QOOH[i];
	n = n_O2_addition_QOOH[i];
	E = E_O2_addition_QOOH[i];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_o2_rem_ooqooh(Radicale r, int numC, double Temp)
{
	int i = PST(r);
	int j = 0;
	if (numC == 5)
		j = 0;
	if (numC > 5 && numC < 12)
		j = 1;
	if (numC > 11 && numC < 15)
		j = 2;
	if (numC > 14)
		j = 3;
	A = A_O2_removal_OOQOOH[i][j];
	n = n_O2_removal_OOQOOH[i][j];
	E = E_O2_removal_OOQOOH[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ooqooh_to_khp(Radicale r1, Radicale r2, Anello a, int numH, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	int k = ring(a);
	A = A_OOQOOH_to_KHP[i][j][k] * double(numH);
	n = n_OOQOOH_to_KHP[i][j][k];
	E = E_OOQOOH_to_KHP[i][j][k];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_khp_decomp(Radicale r, int dist, double Temp)
{
	int i = PST(r);
	int j = dist - 1;
	if (j > 2)
		j = 2;
	A = A_KHP_decomp[i][j];
	n = n_KHP_decomp[i][j];
	E = E_KHP_decomp[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ole_par_roo(Radicale r1, Carbonio c, int numH, double Temp)
{
	int i = PST(r1);
	int j = PST(c);
	A = A_ole_from_ROO[i][j] * double(numH);
	n = n_ole_from_ROO[i][j];
	E = E_ole_from_ROO[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ole_par_ooqooh(Radicale r1, Carbonio c, int numH, double Temp)
{
	int i = PST(r1);
	int j = PST(c);
	A = A_ole_from_OOQOOH[i][j] * double(numH);
	n = n_ole_from_OOQOOH[i][j];
	E = E_ole_from_OOQOOH[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ole_from_beta_qooh(Radicale r1, Radicale r2, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	A = A_ole_from_beta_QOOH[i][j];
	n = n_ole_from_beta_QOOH[i][j];
	E = E_ole_from_beta_QOOH[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ole_from_gamma_qooh(Radicale r1, Radicale r2, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	A = A_ole_from_gamma_QOOH[i][j];
	n = n_ole_from_gamma_QOOH[i][j];
	E = E_ole_from_gamma_QOOH[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ole_from_delta_qooh(Radicale r1, Radicale r2, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	A = A_ole_from_delta_QOOH[i][j];
	n = n_ole_from_delta_QOOH[i][j];
	E = E_ole_from_delta_QOOH[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_pooh2_dec_1(Radicale r1, Radicale r2, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	A = A_POOH2_dec_1[i][j];
	n = n_POOH2_dec_1[i][j];
	E = E_POOH2_dec_1[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_pooh2_dec_2(Radicale r1, Radicale r2, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	A = A_POOH2_dec_2[i][j];
	n = n_POOH2_dec_2[i][j];
	E = E_POOH2_dec_2[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_pooh2_dec_3(Radicale r1, Radicale r2, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	A = A_POOH2_dec_3[i][j];
	n = n_POOH2_dec_3[i][j];
	E = E_POOH2_dec_3[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ether_from_qooh(Radicale r1, Radicale r2, AnelloO a, Carbonio corr, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	int k = ringO(a);
	int z = PSTQ(corr)-1;
	A = A_ether_from_QOOH[i][j][k][z];
	n = n_ether_from_QOOH[i][j][k][z];
	E = E_ether_from_QOOH[i][j][k][z];
	if (A == 0) // the required value with correction is not available, hence the one withoutis going to be used
	{
		std::cerr << "Warning: in v_ether_from_qooh() the value with the correction is not available" << std::endl;
		z = PSTQ(Cs);
		A = A_ether_from_QOOH[i][j][k][z];
		n = n_ether_from_QOOH[i][j][k][z];
		E = E_ether_from_QOOH[i][j][k][z];
	}
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ether_from_pooh2(Radicale r1, Radicale r2, AnelloO a, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	int k = ringO(a);
	A = A_ether_from_POOH2[i][j][k];
	n = n_ether_from_POOH2[i][j][k];
	E = E_ether_from_POOH2[i][j][k];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_beta_dec_qooh(Radicale r1, Radicale r2, double Temp)
{
	int i = PST(r1);
	int j = PST(r2);
	A = A_QOOH_beta_dec[i][j];
	n = n_QOOH_beta_dec[i][j];
	E = E_QOOH_beta_dec[i][j];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_OH_abs(Carbonio c, int numH, int isom, double Temp)
{
	int i = PST(c);
	A = A_OH_abs[i] * double(numH) * double(isom);
	n = n_OH_abs[i];
	E = E_OH_abs[i];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ooqooh_dec_tert(double Temp)
{
	A = A_OOQOOH_dec_tert;
	n = n_OOQOOH_dec_tert;
	E = E_OOQOOH_dec_tert;
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ooqooh_dec_quat(double Temp)
{
	A = A_OOQOOH_dec_quat;
	n = n_OOQOOH_dec_quat;
	E = E_OOQOOH_dec_quat;
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_roo_h_abs(double Temp)
{
	A = A_ROO_h_abs;
	n = n_ROO_h_abs;
	E = E_ROO_h_abs;
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_r_h_abs(double Temp)
{
	A = A_R_h_abs;
	n = n_R_h_abs;
	E = E_R_h_abs;
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_oleooh_dec(Radicale r, double Temp)
{
	int i = PST(r);
	A = A_oleOOH_dec[i];
	n = n_oleOOH_dec[i];
	E = E_oleOOH_dec[i];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_etherooh_dec(Radicale r, double Temp)
{
	int i = PST(r);
	A = A_etherOOH_dec[i];
	n = n_etherOOH_dec[i];
	E = E_etherOOH_dec[i];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_cyc_eth_dec(double Temp)
{
	A = A_cyc_eth_dec;
	n = n_cyc_eth_dec;
	E = E_cyc_eth_dec;
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_allylic_rad_form(Carbonio c, int numH, double Temp)
{
	int i = PST(c);
	A = A_allylic_rad_form[i] * double(numH);
	n = n_allylic_rad_form[i];
	E = E_allylic_rad_form[i];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_alkenyl_ro_form(double Temp)
{
	A = A_alkenyl_RO_form;
	n = n_alkenyl_RO_form;
	E = E_alkenyl_RO_form;
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_alkenyl_ro_dec(Carbonio c, double Temp)
{
	int i = PSTQ(c);
	A = A_alkenyl_RO_dec[i];
	n = n_alkenyl_RO_dec[i];
	E = E_alkenyl_RO_dec[i];
	return A * pow(Temp, n) * exp(-E / R / Temp);
}

double Kinox::v_ald_dec(double Temp)
{
	A = A_ald_dec;
	n = n_ald_dec;
	E = E_ald_dec;
	return A * pow(Temp, n) * exp(-E / R / Temp);
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

