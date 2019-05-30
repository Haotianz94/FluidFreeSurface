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
	assert(Configer::getConfiger()->getString(simulationStr.c_str(), "SceneType", SCENETYPE)); 
	assert(Configer::getConfiger()->getString(simulationStr.c_str(), "RenderType", renderTypeStr));
	assert(Configer::getConfiger()->getBool(simulationStr.c_str(), "Obstacle", OBSTACLE)); 
	assert(Configer::getConfiger()->getBool(simulationStr.c_str(), "FlowIn", FLOWIN)); 

	// RenderType
	if(renderTypeStr == std::string("PARTICLE"))
		RENDERTYPE = PARTICLE;
	else if(renderTypeStr == std::string("PRESSURE"))
		RENDERTYPE = PRESSURE;
	else if(renderTypeStr == std::string("VELOSITYY"))
		RENDERTYPE = VELOSITYY;
	else if(renderTypeStr == std::string("VELOSITYX"))
		RENDERTYPE = VELOSITYX;
	else if(renderTypeStr == std::string("DIVERGENCE"))
		RENDERTYPE = DIVERGENCE;
	else if(renderTypeStr == std::string("FLUIDGRID"))
		RENDERTYPE = FLUIDGRID;
	else if(renderTypeStr == std::string("BLOBBY"))
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