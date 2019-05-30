#include "stdafx.h"
#include "fluidCube3D.h"
#include "display.h"
#include <Eigen/Core>
#include <omp.h>

int main(int argc, char* argv[])
{
	omp_set_num_threads(16);
	Eigen::setNbThreads(16);
	Eigen::initParallel();

	bool CREATEBLOBBY;
	assert(Configer::getConfiger()->getBool("Base", "CreateBlobby", CREATEBLOBBY));

	if(!CREATEBLOBBY)
	{
		PRINT("Program Starting");

		glutInit(&argc, argv);
		//initialize the window 
		initialize();
		
		PRINT("Entering Main Loop");
		glutMainLoop(); //this starts the infinite loop
		PRINT("Exiting Program");
	}
	else
	{
		FluidCube3D *cube = new FluidCube3D();
		cube->createBlobby(1000);
	}

	return 0;
}