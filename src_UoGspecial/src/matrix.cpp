#include "matrix.h"
#include <fstream>
#include <iostream>

Matrix::Matrix()
	: width(0),
	height(0),
	elements(0)
{
}

Matrix::Matrix(int _height, int _width)
	: width(_width),
	height(_height),
	elements(_width* _height)
{
}

double& Matrix::operator()(int row, int column)
{
	if (row < 0 || row > height - 1)
		std::cerr << "Row out of index!" << std::endl;
	if (column < 0 || column > width - 1)
		std::cerr << "Column out of index!" << std::endl;
	return elements[row * width + column];
}

double Matrix::operator()(int row, int column) const
{
	if (row < 0 || row > height - 1)
		std::cerr << "Row out of index!" << std::endl;
	if (column < 0 || column > width - 1)
		std::cerr << "Column out of index!" << std::endl;
	return elements[row * width + column];
}

void Matrix::resize(int _height, int _width)
{
	elements.resize(_width * _height);
	width = _width;
	height = _height;
}

void Matrix::print()
{
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			std::cout << elements[j + i * width] << " ";
		}
		std::cout << std::endl;
	}
}

