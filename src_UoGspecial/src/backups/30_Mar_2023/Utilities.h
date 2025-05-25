#pragma once

#include <iostream>
#include <fstream>
#include <string>


inline void deleteLineNumber(std::string path, int lineNumber) {
	std::string line;
	std::ifstream fin;

	fin.open(path);
	// contents of path must be copied to a temp file then
	// renamed back to the path file
	std::ofstream temp;
	temp.open("temp.txt");

	int counter = 1;
	while (getline(fin, line)) {
		// write all lines to temp other than the line marked for erasing
		if (counter != lineNumber)
			temp << line << std::endl;
		else
			std::cout << "Line removed from " << path << std::endl << "       " << line << std::endl;
		counter++;
	}

	temp.close();
	fin.close();

	// required conversion for remove and rename functions
	const char* p = path.c_str();
	remove(p);
	rename("temp.txt", p);
}