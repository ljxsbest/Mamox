#pragma once
#include <string>

struct reactionComment
{
	std::string rateRule = "";
	std::string humanReadable = "";
	int HMultiplier = 1;
	int multiPathMultiplier = 1;
	int isomerMultiplier = 1;

	reactionComment()
	{
	}

	reactionComment(std::string rateRule_)
	{
		rateRule = rateRule_;
	}

	reactionComment(std::string rateRule_, int Hmul)
	{
		rateRule = rateRule_;
		HMultiplier = Hmul;
	}

	reactionComment(std::string rateRule_, int Hmul, int isomMul)
	{
		rateRule = rateRule_;
		HMultiplier = Hmul;
		isomerMultiplier = isomMul;
	}

	reactionComment(std::string rateRule_, std::string humRead)
	{
		rateRule = rateRule_;
		humanReadable = humRead;
	}

	reactionComment(std::string rateRule_, std::string humRead, int Hmul)
	{
		rateRule = rateRule_;
		humanReadable = humRead;
		HMultiplier = Hmul;
	}

	reactionComment(std::string rateRule_, std::string humRead, int Hmul, int isomMul)
	{
		rateRule = rateRule_;
		humanReadable = humRead;
		HMultiplier = Hmul;
		isomerMultiplier = isomMul;
	}

	std::string comment()
	{
		std::string comm;
		comm.append("RR ");
		comm.append(rateRule);
		comm.append(" M ");
		comm.append(std::to_string(HMultiplier * multiPathMultiplier * isomerMultiplier));
		
		//if (humanReadable.size() > 0)
		//{	
		//	comm.append(" ; ");
		//	comm.append(humanReadable);
		//}
		//if (HMultiplier > 1)
		//{
		//	comm.append(" ; A is multiplied by ");
		//	comm.append(std::to_string(HMultiplier));
		//	comm.append(" because of equivalent Hs");
		//}
		//if (multiPathMultiplier > 1)
		//{
		//	comm.append(" ; A is multiplied by ");
		//	comm.append(std::to_string(multiPathMultiplier));
		//	comm.append(" because of multiple pathways");
		//}
		//if (isomerMultiplier > 1)
		//{
		//	comm.append(" ; A is multiplied by ");
		//	comm.append(std::to_string(isomerMultiplier));
		//	comm.append(" because of multiple isomers");
		//}
		return comm;
	}
};