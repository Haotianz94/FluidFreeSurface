#ifndef _CONFIGER_H_
#define _CONFIGER_H_

#include <string>
#include <vector>

struct Vector2i
{
	Vector2i(int xx, int yy): x(xx), y(yy){}
	Vector2i(){}

	union{
		int val[2];
	};
	struct{
		int x;
		int y;
	};
};

struct Vector2d
{
	Vector2d(double xx, double yy): x(xx), y(yy){}
	Vector2d(){}

	union{
		double val[2];
	};
	struct{
		double x;
		double y;
	};
};


struct Vector3i
{
	Vector3i(int xx, int yy, int zz): x(xx), y(yy), z(zz){}
	Vector3i(){}

	union{
		int val[3];
	};
	struct{
		int x;
		int y;
		int z;
	};
};

struct Vector3d
{
	Vector3d(double xx, double yy, double zz): x(xx), y(yy), z(zz){}
	Vector3d(){}

	union{
		double val[3];
	};
	struct{
		double x;
		double y;
		double z;
	};
};

struct ConfigKeyPair
{
	std::string key;
	std::string value;
};

class ConfigSection
{
private:
	std::string m_secName;
	std::vector<ConfigKeyPair> m_pairs;

public:
	ConfigSection(std::string name): m_secName(name){}
	const std::string& getName() 
	{ 
		return m_secName;
	}
	void addKeyPair(const ConfigKeyPair& pair);
	const std::string* findValue(const char* key) const;
};

class Configer
{
private:
	Configer() 
	{
		init();
	}
	~Configer();

	static Configer* m_configer;
	const static std::string m_configFile;
	std::vector<ConfigSection*> m_sections;

	const ConfigSection* findSection(const char* section) const;
	bool getString(const char* section, const char* key, char* text) const;
	void parseIntoSections(std::string& buffer);

public:
	static Configer* getConfiger();
	void init(const char* configFile = NULL);

	bool getInt(const char* section, const char* key, int& value) const;
	bool getBool(const char* section, const char* key, bool& value) const;
	bool getFloat(const char* section, const char* key, float& value) const;
	bool getDouble(const char* section, const char* key, double& value) const;
	bool getString(const char* section, const char* key, std::string& value) const;
	bool getVec2i(const char* section, const char* key, Vector2i& value) const;
	bool getVec2d(const char* section, const char* key, Vector2d& value) const;
	bool getVec3i(const char* section, const char* key, Vector3i& vec) const;
	bool getVec3d(const char* section, const char* key, Vector3d& vec) const;
	bool getVecXi(const char* section, const char* key, std::vector<int>& vec) const;
	bool getVecXd(const char* section, const char* key, std::vector<double>& vec) const;
	bool getVecXb(const char* section, const char* key, std::vector<bool>& vec) const;
};


#endif