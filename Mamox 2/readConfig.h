// HEADER ONLY LIBRARY
#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "utilities.h"
#include <sstream>


namespace rdcfg {


	inline std::vector<std::string> divideChunks(std::ifstream* infile)
	{
		infile->clear();
		infile->seekg(0);
		std::stringstream buffSS;
		std::string line;
		// remove comments
		while (!infile->eof())
		{
			std::getline(*infile, line, '!');
			buffSS << line;
			std::getline(*infile, line);
		}
		std::string buffS = buffSS.str();
		// remove whitespaces
		//buffS.erase(std::remove_if(buffS.begin(), buffS.end(), ::isspace), buffS.end());
		buffS = UTL::removeSpaces(buffS);
		buffSS.str(buffS);
		std::vector<std::string> chunks;
		std::string ln;
		while (!buffSS.eof())
		{
			std::getline(buffSS, ln, ';');
			chunks.push_back(ln);
		}
		return chunks;
	}

	inline bool isKeywordPresent(std::ifstream* infile, std::string keyword)
	{
		std::vector<std::string> chunks = divideChunks(infile);
		std::string word;
		std::stringstream chunkSS;
		for (int i = 0; i < chunks.size(); i++)
		{
			chunkSS.str(chunks[i]);
			std::getline(chunkSS, word, '=');
			if (word == keyword)
				return true;
		}
		return false;
	}

	inline std::string getRawValue(std::ifstream* infile, std::string keyword)
	{
		infile->clear();
		infile->seekg(0);
		std::vector<std::string> chunks = divideChunks(infile);
		std::string kw;
		std::string value;
		for (int i = 0; i < chunks.size(); i++)
		{
			std::stringstream chunkSS(chunks[i]);
			std::getline(chunkSS, kw, '=');
			std::getline(chunkSS, value, '=');
			if (kw == keyword)
				return value;
		}
		return "";
	}

	inline int readInt(std::ifstream* infile, std::string keyword)
	{
		std::string value = getRawValue(infile, keyword);
		if (UTL::is_int(value))
			return std::stoi(value);
		else
			UTL::fatalError("in readCongif.h, readInt: string is not an integer!");
		
		return -1;
	}
	
	inline double readDouble(std::ifstream* infile, std::string keyword)
	{
		std::string value = getRawValue(infile, keyword);
		if (UTL::is_number(value))
			return std::stof(value);
		else
			UTL::fatalError("in readCongif.h, readDouble: string is not an double!");

		return -1;
	}

	inline std::string readString(std::ifstream* infile, std::string keyword)
	{
		std::string value = getRawValue(infile, keyword);
		return value;
	}

	inline std::vector<double> readDoubleVec(std::ifstream* infile, std::string keyword)
	{
		std::string value = getRawValue(infile, keyword);
		std::stringstream valueSS(value);
		std::vector<double> vec;
		while (!valueSS.eof())
		{
			std::getline(valueSS, value, ',');
			if (UTL::is_number(value))
				vec.push_back(std::stof(value));
			else
				UTL::fatalError("in readCongif.h, readDoubleVec: string properly formatted!");
		}
		return vec;
	}

	inline std::vector<std::string> readStringVec(std::ifstream* infile, std::string keyword)
	{
		std::string value = getRawValue(infile, keyword);
		std::stringstream valueSS(value);
		std::vector<std::string> vec;
		while (!valueSS.eof())
		{
			std::getline(valueSS, value, ',');
			vec.push_back(value);
		}
		return vec;
	}

	//inline std::string readString(std::ifstream* infile, std::string keyword)
	//{
	//
	//}
}