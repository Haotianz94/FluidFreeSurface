#include "stdafx.h"
#include "fluidCube3D.h"
#include "display.h"
#include <fstream>

int main(int argc, char* argv[])
{

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
//1.  h2;  1/ dt
//2.  Store A