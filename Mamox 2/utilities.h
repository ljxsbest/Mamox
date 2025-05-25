// HEADER ONLY LIBRARY
#pragma once
#include <iostream>
#include <string>
#include <stdlib.h>
#define NOMINMAX
#include <Windows.h>
#include <WinBase.h>
#include <direct.h>


namespace UTL {

	/**
	* Print an error message on the terminal in red and terminate the program
	* @param message the message to print
	* @return void
	*/
	inline void fatalError(std::string message);

	/**
	* Print an error message in red on the terminal wihtout terminating the program
	* @param message the message to print
	* @return void
	*/
	inline void error(std::string message);

	/**
	* Print a warning message in yellow on the terminal.
	* @param message the message to print
	* @return void
	*/
	inline void warning(std::string message);

	/**
	* return true if s represents an integer (and can be converted into integer
	* with std::stoi).
	* @param s string to check
	* @return bool the results of the check
	*/
	inline bool is_int(const std::string& s);

	/**
	* Add the element to the vector only if it is not already present. If the elemnt is present
	* return the index of the element in the vector otherwise return -1.
	* @param vec pointer to the vector to which the element has to be addedd
	* @param elem element to add
	* @return integer containing the index of the element if it is already present, otherwise
	*         return -1
	*/
	template<typename T> inline int addUnique(std::vector<T>* vec, T elem);

	/**
	* Concatenate toAdd vector t host vector
	* @param host pointer to the vector that is going to be extended by concatenating
	*			  the toAdd vector
	* @param toAdd the vector to concatenate to host vector
	* @return void
	*/
	template<typename T> inline void concatenate(std::vector<T>* host, std::vector<T>* toAdd);
	
	/**
	* Concatenate to host all the vector of toAdd that are not already present in host
	* @param host pointer to the vector that is going to be extended
	* @param toAdd pointer to the vector to concatenate to the host vector
	* @return void
	*/
	template<typename T> inline void concatenateUnique(std::vector<T>* host,
		std::vector<T>* toAdd);

	/**
	* Check if the two vectors have the same number and tye of elements not 
	* considering the order (i.e. return true for {1, 4, 7, 4} and {4, 4, 7, 1})
	* 
	* @param vec1 pointer to the first vector to check
	* @param vec2 pointer to the second vector to check
	* @return bool true if the two vectors are equivalent
	*/
	template<typename T> inline bool areEquivalent(std::vector<T>* vec1,
		std::vector<T>* vec2);

	/**
	* Check if the element is present in the vector
	* 
	* @param vec1 pointer to the vector
	* @param elem element to check
	* @return bool true if the element is present in the vector
	*/
	template<typename T> inline bool isPresent(std::vector<T>* vec1,
		T elem);

	/**
	* Return the y vector of the linear interpolation of the known function on the
	* provided X vector
	* 
	* @param knownX vector of x values of the function to interpolate (must be sorted
	*			minumun to maximum)
	* @param knownY vector of y values of the function to interpolate
	* @param X vector of x values to get the interpolated values for
	* @return a vector with the interpolated values
	*/
	inline std::vector<double> linInterpolation(std::vector<double> knownX,
		std::vector<double> knownY, std::vector<double> X);

	/**
	* Set the cursor of the console terminal to the beginning of the line
	* 
	* @return void
	*/
	inline void setCursorToLineStart();

	// ~~ is_number
	//  - returns true if s represent a number, either integer or double (and can be converted into
	//    integer or double using std::stoi and std::stod)
	inline bool is_number(const std::string& s);

	// ~~ printEmbeddedString
	//  - print on the terminal a line fillin the whole width of the terminal filled with chars c
	//    containing title in the middle. If title is an empty string it prints just a solid line
	//    of chars c
	inline void printEmbeddedString(char c, std::string title);
	
	// ~~ printCharLine
	//  - print on the terminal a line of chars c filling the whole width of the terminal window
	inline void printCharLine(char c);

	// ~~ printAshtagLine
	//  - print on the terminal a line of ashtags filling the whole width of the terminal window
	inline void printAshtagLine();

	// ~~ printTitle
	//  - print a title of the terminal, composed by three ashtags lines embeddint the string title
	//    in the middle
	inline void printTitle(std::string title);
	
	// ~~ removeSpaces
	//  - remove spaces, tabs and new lines from str and return the updated string. Spaces contained
	//    between " and ' chars are kept
	inline std::string removeSpaces(std::string str);
	
	// ~~ checkFileExistence
	//  - return true if the file exists
	inline bool checkFileExistence(std::string path);

	// ~~ addUniqueString
	//  - add str to vec only if it is not contained already. If str is present in vec return its 
	//    index otherwise return -1
	inline int addUniqueString(std::vector<std::string> vec, std::string str);
	
	// ~~ dirExists
	//  - return true if dirName_in is a directory
	inline bool dirExists(const std::string& dirName_in);

	// ~~ createDirectory
	//	- create directory with path, throw fatal error if it is not possible.
	inline void createDirectory(std::string path);





	inline void fatalError(std::string message)
	{
		//std::cout << "\033[1; 31m" << "ERROR: " << message << std::endl;
		std::cout << "\x1B[31mERROR: " << message << std::endl <<
			"       Program will be closed.\033[0m" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	inline void error(std::string message)
	{
		//std::cout << "\033[1; 31m" << "ERROR: " << message << std::endl;
		std::cout << "\x1B[31mERROR: " << message << std::endl <<
			"       Program will continue to execute.\033[0m" << std::endl;
	}


	inline void warning(std::string message)
	{
		std::cout << "\x1B[33mWARNING: " << message << "\033[0m" << std::endl;
	}

	inline bool is_int(const std::string& s)
	{
		return !s.empty() && std::find_if(s.begin(),
			s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
	}
	inline bool is_number(const std::string& s)
	{
		return !s.empty() && std::find_if(s.begin(),
			s.end(), [](unsigned char c) { return !(std::isdigit(c) || c == '.'); }) == s.end();
	}


	inline void printEmbeddedString(char c, std::string title)
	{
		CONSOLE_SCREEN_BUFFER_INFO csbi;
		int columns, rows;

		GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
		columns = csbi.srWindow.Right - csbi.srWindow.Left + 1;

		if (title.size() > 0)
		{
			std::string temp = " ";
			temp.append(title);
			temp.append(" ");
			title = temp;
		}

		int titleLen = title.size();
		int numChars = std::floor((columns - titleLen) / 2);
		if (numChars < 0)
			numChars = 0;

		for (int i = 0; i < numChars; i++)
			std::cout << c;
		std::cout << title;
		for (int i = numChars + titleLen; i < columns; i++)
			std::cout << c;
		std::cout << std::endl;
	}

	inline void printCharLine(char c)
	{
		UTL::printEmbeddedString(c, "");
	}

	inline void printAshtagLine()
	{
		UTL::printCharLine('#');
	}

	inline void printTitle(std::string title)
	{
		UTL::printAshtagLine();
		UTL::printEmbeddedString('#', title);
		UTL::printAshtagLine();
	}

	inline std::string removeSpaces(std::string str)
	{
		std::string newStr = "";
		std::vector<bool> safeSpaces(str.size(), false);
		bool safezone = false;
		for (int i = 0; i < str.size(); i++)	
		{
			if ((str[i] == '\"' || str[i] == '\'') && safezone == false)
			{
				safezone = true;
				continue;
			}
			if ((str[i] == '\"' || str[i] == '\'') && safezone == true)
			{
				safezone = false;
				continue;
			}
			if (safezone)
			{
				safeSpaces[i] = true;
			}
		}
		for (int i = 0; i < str.size(); i++)
		{
			if ((str[i] == ' ' || str[i] == '\t' || str[i] == '\n') && safeSpaces[i] == false)
				continue;
			if (str[i] == '\"' || str[i] == '\'')
				continue;
			newStr.push_back(str[i]);
		}
		return newStr;
	}

	inline bool checkFileExistence(std::string path)
	{
		std::wstring stemp = std::wstring(path.begin(), path.end());
		LPCWSTR pa = stemp.c_str();
		GetFileAttributes(pa); // from winbase.h
		if (INVALID_FILE_ATTRIBUTES == GetFileAttributes(pa) && GetLastError() == ERROR_FILE_NOT_FOUND)
		{
			return false;
		}
		return true;
	}

	inline int addUniqueString(std::vector<std::string>* vec, std::string str)
	{
		return addUnique(vec, str);
	}

	inline bool dirExists(const std::string& dirName_in)
	{
		DWORD ftyp = GetFileAttributesA(dirName_in.c_str());
		if (ftyp == INVALID_FILE_ATTRIBUTES)
			return false;  //something is wrong with your path!

		if (ftyp & FILE_ATTRIBUTE_DIRECTORY)
			return true;   // this is a directory!

		return false;    // this is not a directory!
	}

	inline void createDirectory(std::string path)
	{
		if (_mkdir(path.c_str()) == -1)
			if (errno == ENOENT)
				UTL::fatalError("Unable to create directory \"" + path + "\"");
	}

	template<typename T> inline int addUnique(std::vector<T>* vec, T elem)
	{
		for (int i = 0; i < vec->size(); i++)
		{
			if ((*vec)[i] == elem)
				return i;
		}
		vec->push_back(elem);
		return -1;
	}

	template<typename T> inline void concatenate(std::vector<T>* host, std::vector<T>* toAdd)
	{
		host->insert(host->end(), toAdd->begin(), toAdd->end());
	}

	template<typename T> inline void concatenateUnique(std::vector<T>* host,
		std::vector<T>* toAdd)
	{
		for (auto& elem : *toAdd)
		{
			addUnique(host, elem);
		}
	}

	template<typename T> inline bool areEquivalent(std::vector<T>* vec1,
		std::vector<T>* vec2)
	{
		if (vec1->size() != vec2->size())
			return false;

		std::vector<bool> alreadyFound(vec2->size(), false);
		for (int i = 0; i < vec1->size(); i++)
		{
			bool found = false;
			for (int j = 0; j < vec2->size(); j++)
			{
				if (alreadyFound[j] == true)
					continue;
				if ((*vec1)[i] == (*vec2)[j])
				{
					found = true;
					alreadyFound[j] = true;
					break;
				}
			}
			if (found == false)
				return false;
		}
		return true;
	}

	template<typename T> inline bool isPresent(std::vector<T>* vec1,
		T elem)
	{
		for (auto& el : *vec1)
			if (el == elem)
				return true;
		return false;
	}

	inline std::vector<double> linInterpolation(std::vector<double> knownX,
		std::vector<double> knownY, std::vector<double> X)
	{
		std::vector<double> Y(X.size());
		for (int i = 0; i < X.size(); i++)
		{
			if (X[i] < knownX[0])
			{
				warning("In UTL::linInterpolation a value is outside the range of provided values, extrapolation has been used");
				Y[i] = knownY[0] + (knownY[1] - knownY[0]) * 
					(X[i] - knownX[0]) / (knownX[1] - knownX[0]);
			}
			else if (X[i] > knownX[knownX.size()-1])
			{
				warning("In UTL::linInterpolation a value is outside the range of provided values, extrapolation has been used");
				Y[i] = knownY[knownX.size() - 2] + (knownY[knownX.size() - 1] - 
					knownY[knownX.size() - 2]) * (X[i] - knownX[knownX.size() - 2]) 
					/ (knownX[knownX.size() - 1] - knownX[knownX.size() - 2]);
			}
			else
			{
				int ind = 0;
				for (int j = 0; j < knownX.size(); j++)
				{
					if (X[i] > knownX[j])
					{
						ind = j;
						break;
					}
				}
				Y[i] = knownY[ind] + (knownY[ind+1] - knownY[ind]) *
					(X[i] - knownX[ind]) / (knownX[ind+1] - knownX[ind]);
			}
		}
		return Y;
	}

	inline void setCursorToLineStart()
	{
		CONSOLE_SCREEN_BUFFER_INFO conInf;
		GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &conInf);
		COORD coord = conInf.dwCursorPosition;
		coord.X = 0;
		SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), coord);
	}
}

