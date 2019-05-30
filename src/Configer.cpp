//#####################################################################
// Copyright 2017. Haotian Zhang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Configer.h"
#include "Logger.h"
#include <fstream>
#include <cstring>

using namespace std;

#define STRING_LEN 80
const string Configer::m_configFile = "../config.ini";
Configer* Configer::m_configer = NULL;


Configer::~Configer()
{
	for(auto& sec : m_sections)
		delete sec;
	m_sections.clear();
	if(m_configer)
		delete m_configer;
}

Configer* Configer::getConfiger()
{
	if(!m_configer)
		m_configer = new Configer();
	return m_configer;
}

void Configer::init(const char* configFile)
{
	const char* realPath = configFile? configFile : m_configFile.c_str();
	ifstream fin(realPath);
	assert(!fin.fail() && "Can not open config.ini!");

	string buffer;
	while(!fin.eof())
	{
		getline(fin, buffer);
		parseIntoSections(buffer);
	}
	fin.close();
}

void Configer::parseIntoSections(string& line)
{
	const char* tmpBuff = strchr(line.c_str(), '[');

	if (tmpBuff)
	{
		const char* tmpEnd = strchr(tmpBuff, ']');
		unsigned startPos = static_cast<unsigned>(tmpBuff - line.c_str()) + 1;
		unsigned numChar = static_cast<unsigned>(tmpEnd - line.c_str()) - 1;
		string secName = line.substr(startPos, numChar);
		ConfigSection* newSec = new ConfigSection(secName);
		m_sections.push_back(newSec);
	}
	else
	{
		const char* tmpPairBuff = strchr(line.c_str(), '=');
		if (tmpPairBuff)
		{	
			ConfigKeyPair pair;
			unsigned numChar = static_cast<unsigned>(tmpPairBuff - line.c_str());
			pair.key = line.substr(0, numChar);

			const char* valueField = tmpPairBuff + 1;
			unsigned startPos = static_cast<unsigned>(valueField - line.c_str());
			
			// numChar = static_cast<unsigned>(line.length() - startPos);
			numChar = static_cast<unsigned>(line.length() - startPos - 1);

			pair.value = line.substr(startPos, numChar);
			if (m_sections.size())
			{
				ConfigSection* currentSection = m_sections.back();
				currentSection->addKeyPair(pair);
			}

		}
	}
}

void ConfigSection::addKeyPair(const ConfigKeyPair& pair)
{
	m_pairs.push_back(pair);
}

bool Configer::getString(const char* section, const char* key, char* text) const
{
	const ConfigSection* sec = findSection(section);
	if(!sec)
		return false;
	const string* stringValue = sec->findValue(key);
	if(!stringValue)
		return false;
	strncpy(text, stringValue->c_str(), STRING_LEN);
	return true;
}

const ConfigSection* Configer::findSection(const char* section) const
{
	for(auto& sec: m_sections)
	{
		if(!strcmp(sec->getName().c_str(), section))
		{
			return sec;
		}
	}
	return NULL;
}

const string* ConfigSection::findValue(const char* key) const
{
	for(auto& pair : m_pairs)
	{
		if(!strcmp(pair.key.c_str(), key))
		{
			return &pair.value;
		}
	}
	return NULL;
}

bool Configer::getInt(const char* section, const char* key, int& value) const
{
	char text[STRING_LEN];
	if(getString(section, key, text))
	{
		value = atoi(text);
		LOG << "[Configer]: " << key << " = " << value <<"\n";
		return true;
	}
	return false;
}

bool Configer::getBool(const char* section, const char* key, bool& value) const
{
	char text[STRING_LEN];
	if(getString(section, key, text))
	{
		value = ( atoi(text) == 1);
		LOG << "[Configer]: " << key << " = " << value <<"\n";
		return true;
	}
	return false;
}

bool Configer::getFloat(const char* section, const char* key, float& value) const
{
	char text[STRING_LEN];
	if(getString(section, key, text))
	{
		value = atof(text);
		LOG << "[Configer]: " << key << " = " << value <<"\n";
		return true;
	}
	return false;
}

bool Configer::getDouble(const char* section, const char* key, double& value) const
{
	char text[STRING_LEN];
	if(getString(section, key, text))
	{
		value = atof(text);
		LOG << "[Configer]: " << key << " = " << value <<"\n";
		return true;
	}
	return false;
}

bool Configer::getString(const char* section, const char* key, string& value) const
{
	char text[STRING_LEN];
	if(getString(section, key, text))
	{
		value = string(text);
		LOG << "[Configer]: " << key << " = " << value.c_str() <<"\n";
		return true;
	}
	return false;
}

bool Configer::getVec3i(const char* section, const char* key, Vector3i& vec) const
{
	vector<double> vec_d;
	if( getVecXd(section, key, vec_d) && vec_d.size() == 3)
	{
		vec = Vector3i((int)vec_d[0], (int)vec_d[1], (int)vec_d[2]);
		return true;
	}
	return false;
}

bool Configer::getVec3d(const char* section, const char* key, Vector3d& vec) const
{
	vector<double> vec_d;
	if( getVecXd(section, key, vec_d) && vec_d.size() == 3)
	{
		vec = Vector3d(vec_d[0], vec_d[1], vec_d[2]);
		return true;
	}
	return false;
}


bool Configer::getVec2i(const char* section, const char* key, Vector2i& vec) const
{
	vector<double> vec_d;
	if( getVecXd(section, key, vec_d) && vec_d.size() == 2)
	{
		vec = Vector2i((int)vec_d[0], (int)vec_d[1]);
		return true;
	}
	return false;
}

bool Configer::getVec2d(const char* section, const char* key, Vector2d& vec) const
{
	vector<double> vec_d;
	if( getVecXd(section, key, vec_d) && vec_d.size() == 2)
	{
		vec = Vector2d(vec_d[0], vec_d[1]);
		return true;
	}
	return false;
}

bool Configer::getVecXd(const char* section, const char* key, vector<double>& vec) const
{
	char text[STRING_LEN];
	if(getString(section, key, text))
	{
		int i = 0, s = 1, l = strlen(text);
		while(i < l)
		{
			if(text[i] == ',')
			{
				vec.push_back(atof(text + s));
				s = i + 1;
			}
			else if(text[i] == ')')
			{
				vec.push_back(atof(text + s));
				LOG << "[Configer]: " << key << " = " << '('; 
				for(int j = 0; j < vec.size(); j++)
				{
					LOG << vec[j];
					if(j != vec.size()-1)
						LOG << ", ";
				}
				LOG << ')' <<"\n";
 				return true;	
			}
			i ++;
		}
	}
	return false;
}

bool Configer::getVecXi(const char* section, const char* key, vector<int>& vec) const
{
	vector<double> vec_d;
	if( getVecXd(section, key, vec_d) )
	{
		for(auto& v : vec_d)
			vec.push_back((int)v);
		return true;
	}
	return false;
}

bool Configer::getVecXb(const char* section, const char* key, vector<bool>& vec) const
{
	vector<double> vec_d;
	if( getVecXd(section, key, vec_d) )
	{
		for(auto& v : vec_d)
			vec.push_back( v == 1);
		return true;
	}
	return false;
}
