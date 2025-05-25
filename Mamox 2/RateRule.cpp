#include "RateRule.h"

RateRule::RateRule(int numDim, bool isSymm)
{
	isSymmetric = isSymm;
	numDimensions = numDim;
	if (numDim == 0)
		zeroDim = true;
}

void RateRule::addRate(std::string ID, std::vector<std::string> param, double A, double n, double E)
{
	Rate newRate(ID, param, A, n, E);
	rates.push_back(newRate);
}

std::vector<double> RateRule::returnRates(std::vector<std::string> param, std::string *ident)
{
	if (zeroDim)
	{
		if (rates.size() == 1)
		{
			*ident = rates[0].ID;
			return  std::vector<double> {rates[0].A, rates[0].n, rates[0].E};
		}
		else
		{
			std::cerr << "ERROR in returnRates: rate rule not found" << std::endl;
			return std::vector<double> {0.0, 0.0, 0.0};
		}
	}
	int index = -1;
	for (int i = 0; i < rates.size(); i++)
	{
		bool match = true;
		if (isSymmetric == false)
		{
			for (int j = 0; j < numDimensions; j++)
			{
				if (param[j] != rates[i].parameters[j])
				{
					match = false;
					continue;
				}
			}
		}
		else
		{
			std::vector<bool> alreadyMatched(param.size(), false);
			for (int j = 0; j < param.size(); j++)
			{
				for (int k = 0; k < param.size(); k++)
				{
					if (param[j] == rates[i].parameters[k] && alreadyMatched[k] == false)
					{
						alreadyMatched[k] = true;
						break;
					}
				}
			}
			for (int j = 0; j < alreadyMatched.size(); j++)
			{
				if (alreadyMatched[j] == false)
				{
					match = false;
					break;
				}
			}
		}

		if (match)
		{
			index = i;
			break;
		}
	}

	if (index == -1)
	{
		std::cerr << "ERROR in returnRates: rate rule not found" << std::endl;
		return std::vector<double> {0.0, 0.0, 0.0};
	}
	*ident = rates[index].ID;
	return std::vector<double> {rates[index].A, rates[index].n, rates[index].E};
}


std::vector<double> RateRule::returnRates(std::vector<std::string> param)
{
	std::string dummy;
	return returnRates(param, &dummy);
}