#include "stdafx.h"
#include "fluidCube3D.h"
#include "display.h"
#include "trianglemesh.h"
#include "quadmesh.h"
#include <Eigen/Core>
#include <omp.h>

int main(int argc, char* argv[])
{
	 omp_set_num_threads(16);
	 Eigen::setNbThreads(16);
	 Eigen::initParallel();

	 // std::string obj_path;
	 // assert(Configer::getConfiger()->getString("Volcano", "ObjectPath", obj_path));
	 // QuadMesh mesh(obj_path);
	// mesh.refine();
	// mesh.dumpObj("../obj/lava_right_refine.obj");
	 // return 0;

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
		cube->createBlobby();
	}

	return 0;
}
