#include "stdafx.h"
#include "fluidCube.h"


FluidCube::FluidCube()
{
	bool Simulation2D;
	std::string renderTypeStr;
	std::string simulationStr;
	assert(Configer::getConfiger()->getBool("Base", "DebugPrint", DEBUGPRINT)); 
	assert(Configer::getConfiger()->getBool("Base", "CreateBlobby", CREATEBLOBBY)); 
	assert(Configer::getConfiger()->getBool("Base", "Simulation2D", Simulation2D));
	if(Simulation2D)
		simulationStr = std::string("Simulation2D");
	else
		simulationStr = std::string("Simulation3D");

	assert(Configer::getConfiger()->getInt(simulationStr.c_str(), "GridSize", GRIDSIZE));
	assert(Configer::getConfiger()->getInt(simulationStr.c_str(), "ParticlePerGrid", PARTICLEPERGRID));
	assert(Configer::getConfiger()->getFloat(simulationStr.c_str(), "Gravity", GRAVITY)); 
	assert(Configer::getConfiger()->getFloat(simulationStr.c_str(), "Viscosity", VISCOSITY)); 
	assert(Configer::getConfiger()->getFloat(simulationStr.c_str(), "FrameRate", FRAMERATE)); 
	assert(Configer::getConfiger()->getInt(simulationStr.c_str(), "MaxIteration", MAXITERATION));
	assert(Configer::getConfiger()->getString(simulationStr.c_str(), "RenderType", renderTypeStr));
	assert(Configer::getConfiger()->getString(simulationStr.c_str(), "SceneType", SCENETYPE)); 
	assert(Configer::getConfiger()->getString(simulationStr.c_str(), "ObstacleType", OBSTACLETYPE)); 
	assert(Configer::getConfiger()->getString(simulationStr.c_str(), "FlowInType", FLOWINTYPE)); 

	// RenderType
	if(renderTypeStr.compare("PARTICLE") == 0)
		RENDERTYPE = PARTICLE;
	else if(renderTypeStr.compare("PRESSURE") == 0)
		RENDERTYPE = PRESSURE;
	else if(renderTypeStr.compare("VELOSITYY") == 0)
		RENDERTYPE = VELOSITYY;
	else if(renderTypeStr.compare("VELOSITYX") == 0)
		RENDERTYPE = VELOSITYX;
	else if(renderTypeStr.compare("DIVERGENCE") == 0)
		RENDERTYPE = DIVERGENCE;
	else if(renderTypeStr.compare("FLUIDGRID") == 0)
		RENDERTYPE = FLUIDGRID;
	else if(renderTypeStr.compare("BLOBBY") == 0)
		RENDERTYPE = BLOBBY;
	else
	{
		PRINT("RenderType not known!");
		exit(0);
	}
}


FluidCube::~FluidCube()
{
	delete [] div;
	delete [] type;
	delete [] type0;
	for(int i = 0; i < NUMGRID; i++)
		delete[] invertedList[i];
	delete [] invertedList;

	delete [] pos2index;
	for(int i = 0; i < NUMGRID; i++)
		delete[] neighbor[i];
	delete [] neighbor;
	delete [] neighNoneSolid;
	delete [] neighAir;

	delete [] fai_b;
	delete [] fai_f;
}