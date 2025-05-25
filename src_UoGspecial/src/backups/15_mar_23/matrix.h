#pragma once
#include <vector>
#include <stdexcept>

class Matrix
{
public:
	Matrix(int _height, int _width);
	Matrix();
	
	double& operator()(int row, int column);
	double operator()(int row, int column) const;
	
	void resize(int _height, int _width);

	void print();

private:
	std::vector<double> elements;
	int width;
	int height;
};


