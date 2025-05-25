#include "XML.h"
#include <sstream>

XML::XML(std::string fileName_, char mode_)
{
	fileName = fileName_;
	mode = mode_;

	if (mode == 'n')
	{
		file.open(fileName, ios_base::out);
		file << "<?xml version=\"1.0\" encoding=\"UTF - 8\"?>" << std::endl;
		file << "<list>" << std::endl;
	}
	else if (mode == 'r')
	{
		file.open(fileName, ios_base::in);
	}

}

void XML::addMolecule(Molecola* mol, std::string name)
{
	//file << "\t<molecule>" << std::endl;
	//
	//file << "\t\t<name>" << name << "</name>" << std::endl;
	//
	//file << "\t\t<size>" << mol->size() << "</size>" << std::endl;
	//file << "\t\t<structure>";
	//for (int i = 0; i < mol->size(); i++)
	//	for (int j = 0; j < 5; j++)
	//		file << mol->struttura[i + 1][j + 1] << " ";
	//file << "</structure>" << std::endl;
	//
	//file << "\t\t<numberOlefins>" << mol->numberOle() << "</numberOlefins>" << std::endl;
	//file << "\t\t<olefinsBonds>";
	//for (int i = 0; i < mol->numberOle(); i++)
	//	file << mol->ole[i + 1][0] << " " << mol->ole[i + 1][1] << " ";
	//file << "</olefinsBonds>" << std::endl;
	//
	//file << "\t\t<numberEthers>" << mol->numberEthero() << "</numberEthers>" << std::endl;
	//file << "\t\t<ethersBonds>";
	//for (int i = 0; i < mol->numberEthero(); i++)
	//	file << mol->ete[i + 1][0] << " " << mol->ete[i + 1][1] << " ";
	//file << "</ethersBonds>" << std::endl;
	//
	//file << "\t\t<numberKetos>" << mol->numberEthero() << "</numberKetos>" << std::endl;
	//file << "\t\t<ketoPositions>";
	//for (int i = 0; i < mol->numberEthero(); i++)
	//	file << mol->doppioO[i + 1] << " ";
	//file << "</ketoPositions>" << std::endl;
	//
	//file << "\t</molecule>" << std::endl;
}

void XML::close()
{
	if (mode == 'n')
	{
		file << "</list>" << std::endl;
		file.close();
	}
	else if (mode == 'r')
	{
		file.close();
	}
}

void XML::getMolecules(std::vector<Molecola>* vec, std::vector<std::string>* names)
{
	bool fileEnded = false;
	std::string keyword;
	std::string word;

	while (!fileEnded)
	{
		keyword = getNextKeyword();
		if (keyword == "molecule")
		{
			std::string molName;
			goToKeyword("name");
			std::getline(file, molName, '<');


			goToKeyword("size");
			std::getline(file, word, '<');
			int size = stoi(word);
			std::vector<int> structure(0);
			if (size > 0)
			{
				goToKeyword("structure");
				for (int i = 0; i < size; i++)
				{
					for (int j = 0; j < 5; j++)
					{
						std::getline(file, word, ' ');
						structure.push_back(stoi(word));
					}
				}
			}
			else
				continue;

			goToKeyword("numberOlefins");
			std::getline(file, word, '<');
			int numberOlefins = stoi(word);
			std::vector<int> olefinsBonds(0);
			goToKeyword("olefinsBonds");
			for (int i = 0; i < numberOlefins; i++)
			{
				std::getline(file, word, ' ');
				olefinsBonds.push_back(stoi(word));
				std::getline(file, word, ' ');
				olefinsBonds.push_back(stoi(word));
			}

			goToKeyword("numberEthers");
			std::getline(file, word, '<');
			int numberEthers = stoi(word);
			std::vector<int> ethersBonds(0);
			goToKeyword("ethersBonds");
			for (int i = 0; i < numberEthers; i++)
			{
				std::getline(file, word, ' ');
				ethersBonds.push_back(stoi(word));
				std::getline(file, word, ' ');
				ethersBonds.push_back(stoi(word));
			}

			goToKeyword("numberKetos");
			std::getline(file, word, '<');
			int numberKetos = stoi(word);
			std::vector<int> ketosPositions(0);
			goToKeyword("ketoPositions");
			for (int i = 0; i < numberKetos; i++)
			{
				std::getline(file, word, ' ');
				ketosPositions.push_back(stoi(word));
			}

			Molecola mol(size, structure, numberOlefins, olefinsBonds, numberEthers,
				ethersBonds, numberKetos, ketosPositions);
			//std::cout << mol << std::endl;
			vec->push_back(mol);
			names->push_back(molName);
			/*
			std::cout << "NAME = " << molName << std::endl;

			std::cout << "	size		= " << size << std::endl;
			std::cout << "	structure	= ";
			for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < 5; j++)
					std::cout << structure[i * 5 + j] << " ";
				std::cout << std::endl;
				if(i != size-1)
					std::cout << "			  ";
			}

			std::cout << "	numOlefins	= " << numberOlefins << std::endl;
			std::cout << "	oleBonds	= ";
			for (int i = 0; i < numberOlefins; i++)
			{
				std::cout << olefinsBonds[i * 2] << " " << olefinsBonds[i * 2 + 1] << " ";
				std::cout << std::endl;
				if(i != numberOlefins-1)
					std::cout << "			  ";
			}
			if (numberOlefins == 0)
				std::cout << std::endl;

			std::cout << "	numEthers	= " << numberEthers << std::endl;
			std::cout << "	ethersBonds	= ";
			for (int i = 0; i < numberEthers; i++)
			{
				std::cout << ethersBonds[i * 2] << " " << ethersBonds[i * 2 + 1] << " ";
				std::cout << std::endl;
				if (i != numberEthers - 1)
					std::cout << "			  ";
			}
			if (numberEthers == 0)
				std::cout << std::endl;

			std::cout << "	numKetos	= " << numberKetos << std::endl;
			std::cout << "	ketosPos	= ";
			for (int i = 0; i < numberKetos; i++)
			{
				std::cout << ketosPositions[i] << " ";
				std::cout << std::endl;
				if (i != numberKetos - 1)
					std::cout << "			  ";
			}
			if (numberKetos == 0)
				std::cout << std::endl;
			std::cout << std::endl;
			*/
		}
		if (keyword == "/list")
			fileEnded = true;
	}
}

std::string XML::getNextKeyword()
{
	std::string word;
	std::getline(file, word, '>');
	std::stringstream ssword(word);
	std::getline(ssword, word, '<');
	std::getline(ssword, word);
	return word;
}

void XML::goToKeyword(std::string keyword)
{
	while (true)
	{
		if (getNextKeyword() == keyword)
			break;
	}
}