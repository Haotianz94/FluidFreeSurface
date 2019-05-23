#include "stdafx.h"
#include "fluidCube3D.h"
#include "display.h"
#include <fstream>
#include <Eigen/Core>
#include <omp.h>

int main(int argc, char* argv[])
{
	omp_set_num_threads(16);
	Eigen::setNbThreads(16);
	Eigen::initParallel();

#ifndef CREATEBLOBBY
	PRINT("Program Starting");

	glutInit(&argc, argv);
	//initialize the window 
	initialize();
	
	PRINT("Entering Main Loop");
	glutMainLoop(); //this starts the infinite loop
	PRINT("Exiting Program");
#else
	FluidCube3D *cube = new FluidCube3D(VISCOSITY, FRAMERATE, MYSCENE, MYRENDER);
	cube->createBlobby(1000);

#endif
	return 0;
}