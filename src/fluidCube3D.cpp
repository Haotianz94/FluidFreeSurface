#include "stdafx.h"
#include "fluidCube3D.h"
#include "CIsoSurface.h"
#include "quadmesh.h"

#include <GL/freeglut.h>
#include <memory.h>
#include <math.h>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <omp.h>
#include <sys/stat.h>


#define IX(x, y, z) ((x) + (y)*(NUMGRIDX+2) + (z)*(NUMGRIDX+2)*(NUMGRIDY+2) )
#define DISTANCE(x1, y1, z1, x2, y2, z2) ( sqrtf((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) )
#define DISTANCE2(x1, y1, z1, x2, y2, z2) ( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) )
#define BOUNDED(x, y, z) ( (x >= 1 && x <= NUMGRIDX && y >= 1 && y <= NUMGRIDY && z >= 1 && z <= NUMGRIDZ)? true : false)

extern float px;
extern float py;
extern float pz;


FluidCube3D::FluidCube3D()
{
	float CUBELENGTH;
	assert(Configer::getConfiger()->getInt("Simulation3D", "NumGridX", NUMGRIDX)); 
	assert(Configer::getConfiger()->getInt("Simulation3D", "NumGridY", NUMGRIDY)); 
	assert(Configer::getConfiger()->getInt("Simulation3D", "NumGridZ", NUMGRIDZ)); 
	assert(Configer::getConfiger()->getFloat("Simulation2D", "CubeLength", CUBELENGTH));

	NUMGRID = (NUMGRIDX+2) * (NUMGRIDY+2) * (NUMGRIDZ+2);
	h = CUBELENGTH / NUMGRIDX;
	h2 = h * h;
	hi = 1 / h;
	frameTime = 1.0 / FRAMERATE;
	ctime = 0;
	totalTime = 0;
	iteration = 0;

	max_vx = 0;
	max_vy = 0;
	max_vz = 0;
	max_v = 0;
	max_p = 0;

	Vx = new float [NUMGRID]; 
	Vy = new float [NUMGRID]; 
	Vz = new float [NUMGRID]; 
	Vx0 = new float [NUMGRID]; 
	Vy0 = new float [NUMGRID]; 
	Vz0 = new float [NUMGRID]; 
	div = new float [NUMGRID];
	memset(Vx, 0, sizeof(float) * NUMGRID);
	memset(Vy, 0, sizeof(float) * NUMGRID);
	memset(Vz, 0, sizeof(float) * NUMGRID);
	memset(Vx0, 0, sizeof(float) * NUMGRID);
	memset(Vy0, 0, sizeof(float) * NUMGRID);
	memset(Vz0, 0, sizeof(float) * NUMGRID);
	memset(div, 0, sizeof(float) * NUMGRID);
	type = new GridType [NUMGRID];
	type0 = new GridType [NUMGRID]; 
	invertedList = new std::vector<int>* [NUMGRID];

	//Advection using BFECC
	fai_b = new float [NUMGRID];
	fai_f = new float [NUMGRID];
	memset(fai_b, 0, sizeof(float) * NUMGRID);
	memset(fai_f, 0, sizeof(float) * NUMGRID);

	//Projection using Conjugate Gradient
	dir[0] = Eigen::Vector3i(0, 0, -1);
	dir[1] = Eigen::Vector3i(0, -1, 0);
	dir[2] = Eigen::Vector3i(-1, 0, 0);
	dir[3] = Eigen::Vector3i(1, 0, 0);
	dir[4] = Eigen::Vector3i(0, 1, 0);
	dir[5] = Eigen::Vector3i(0, 0, 1);
	pos2index = new int [NUMGRID]; 
	neighNoneSolid = new int [NUMGRID];
	neighAir = new int [NUMGRID];
	neighbor = new int* [NUMGRID];
	for(int i = 0; i < NUMGRID; i++)
	{
		neighbor[i] = new int[6];
		pos2index[i] = -1;
		type[i] = SOLID;
		type0[i] = SOLID;
		invertedList[i] = new std::vector<int>();
	}
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
				type[IX(x,y,z)] = AIR;

	//Blobby
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++)
				dir2[i + j*3 + k*9] = Eigen::Vector3i(i-1, j-1, k-1);

	//Extrapolate
	layer = new int [NUMGRID];
	for(int i = 0; i < NUMGRID; i++)
		layer[i] == -1;

	initFluid();
	initSolid();
}

FluidCube3D::~FluidCube3D()
{
	delete [] Vx;
	delete [] Vy;
	delete [] Vz;
	delete [] Vx0;
	delete [] Vy0;
	delete [] Vz0;
}

void FluidCube3D::initFluid()
{
	srand(time(0));
	originFluid = 0;
	fluidNum = 0;
	if(SCENETYPE.compare("CUBEFALL") == 0)
	{
		for(int z = NUMGRIDZ/3.0; z <= NUMGRIDZ/3.0*2; z++)
			for(int y = NUMGRIDY/4.0; y <= NUMGRIDY/2.0; y++)
				for(int x = NUMGRIDX/3.0; x <= NUMGRIDX/3.0*2; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
	}
	else if(SCENETYPE.compare("SPHEREFALL") == 0)
	{
		int cx = NUMGRIDX/2;
		int cy = NUMGRIDY/2;
		int cz = NUMGRIDZ/2;
		float R = NUMGRIDX/3;
		for(int z = 1; z <= NUMGRIDZ; z++)
			for(int y = 1; y <= NUMGRIDY; y++)
				for(int x = 1; x <= NUMGRIDX; x++)
					if(DISTANCE(x, y, z, cx, cy, cz) <= R)
					{
						originFluid ++;
						fillParticleInGrid(x, y, z);
					}
	}
	else if(SCENETYPE.compare("CONTAINER") == 0)
	{
		
		for(int z = NUMGRIDZ/2.0-2; z <= NUMGRIDZ/2.0+2; z++)
			for(int y = NUMGRIDY/4.0+5; y <= NUMGRIDY/4.0+10; y++)
				for(int x = NUMGRIDX/2.0-2; x <= NUMGRIDX/2.0+2; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
		

		for(int z = 1; z <= NUMGRIDZ; z++)
			for(int y = 1; y <= NUMGRIDY/4.0; y++)
				for(int x = 1; x <= NUMGRIDX; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
	}
	else if(SCENETYPE.compare("DAMBREAK") == 0)
	{
		for(int z = 1; z <= NUMGRIDZ/4.0; z++)
			for(int y = 1; y <= NUMGRIDY/3.0*2; y++)
				for(int x = 1; x <= NUMGRIDX; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
	}
	else if(SCENETYPE.compare("DOUBLEDAM") == 0)
	{
		for(int z = 1; z <= NUMGRIDZ/4.0; z++)
			for(int y = 1; y <= NUMGRIDY/3.0*2; y++)
				for(int x = 1; x <= NUMGRIDX; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
		
		for(int z = NUMGRIDZ/4.0*3; z <= NUMGRIDZ; z++)
			for(int y = 1; y <= NUMGRIDY/3.0*2; y++)
				for(int x = 1; x <= NUMGRIDX; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
	}
	else if(SCENETYPE.compare("NONE") == 0)
	{

	}
	// else
	// {
	// 	PRINT("SceneType not known!");
	// 	exit(0);
	// }
}

void FluidCube3D::initSolid()
{	
	if(SCENETYPE.compare("VOLCANO") == 0 || SCENETYPE.compare("VOLCANOBURST") == 0)
	{
		std::string obj_path;
		assert(Configer::getConfiger()->getString("Volcano", "VolcanoPath", obj_path));
		volcano = new QuadMesh(obj_path);
		volcano->normalize(NUMGRIDX);
		volcano->dumpObj("../obj/volcano_simulate.obj");
		volcano->calculateHeightMap(NUMGRIDX);
		for(int x = 1; x <= NUMGRIDX; x++)
			for(int z = 1; z <= NUMGRIDZ; z++)
				for(int y = 1; y <= volcano->getHeight(x-1, z-1)+1; y++)
					type[IX(x, y, z)] = type0[IX(x, y, z)] = SOLID;
	}
	else if(SCENETYPE.compare("GROUND") == 0)
	{
		std::string obj_path;
		assert(Configer::getConfiger()->getString("Volcano", "GroundPath", obj_path));
		ground = new QuadMesh(obj_path);
		ground->normalize(NUMGRIDX);
		ground->dumpObj("../obj/ground_simulate.obj");
		ground->calculateHeightMap(NUMGRIDX);
		for(int x = 1; x <= NUMGRIDX; x++)
			for(int z = 1; z <= NUMGRIDZ; z++)
				for(int y = 1; y <= ground->getHeight(x-1, z-1)+1; y++)
					type[IX(x, y, z)] = type0[IX(x, y, z)] = SOLID;
	}
	else if(OBSTACLETYPE.compare("CENTERWALL") == 0)
	{
		for(int z = NUMGRIDZ/16.0*7; z <= NUMGRIDZ/16.0*9; z++)
			for(int y = 1; y <= NUMGRIDY/3.0*2; y++)
				for(int x = NUMGRIDX/4.0; x <= NUMGRIDX/4.0*3; x++)
					type[IX(x, y, z)] = type0[IX(x, y, z)] = SOLID;
	}
}

void FluidCube3D::addFlowIn()
{
	if(SCENETYPE.compare("VOLCANO") == 0)
	{
		for(int z = NUMGRIDZ/2.0 - 30; z <= NUMGRIDZ/2.0 + 30; z+=1)
			for(int x = NUMGRIDX/2.0 - 30; x <= NUMGRIDX/2.0 + 30; x+=1)
			{
				int height = volcano->getHeight(x-1, z-1) + 2;
				type[IX(x, height, z)] = type0[IX(x, height, z)] = FLUIDIN;
				fillParticleInGrid(x, height, z);
				Vy[IX(x, height, z)] = Vy[IX(x, height+1, z)] = 4;
			}

		// for(int z = NUMGRIDZ/2.0 - 80; z <= NUMGRIDZ/2.0 + 20; z+=1)
		// 	for(int x = NUMGRIDX/2.0 - 60; x <= NUMGRIDX/2.0 + 40; x+=1)
		// 	{
		// 		int height = volcano->getHeight(x-1, z-1) + 2;
		// 		type[IX(x, height, z)] = type0[IX(x, height, z)] = FLUIDIN;
		// 		fillParticleInGrid(x, height, z);
		// 		Vy[IX(x, height, z)] = Vy[IX(x, height+1, z)] = 2;
		// }
	}
	if(SCENETYPE.compare("VOLCANOBURST") == 0)
	{	
		// srand(0);
		// for(int z = NUMGRIDZ/2.0 - 5; z <= NUMGRIDZ/2.0 + 5; z+=5)
		// 	for(int x = NUMGRIDX/2.0 - 5; x <= NUMGRIDX/2.0 + 5; x+=5)
		// 	{
		// 		int height = volcano->getHeight(x-1, z-1) + 2;
		// 		type[IX(x, height, z)] = type0[IX(x, height, z)] = FLUIDIN;
		// 		fillParticleInGrid(x, height, z);
		// 		Vy[IX(x, height, z)] = Vy[IX(x, height+1, z)] = rand() % 10 + 15;
		// 		Vx[IX(x, height, z)] = Vx[IX(x, height+1, z)] = rand() % 6 - 3;
		// 		Vz[IX(x, height, z)] = Vz[IX(x, height+1, z)] = rand() % 6 - 3;
		//  }

		// srand(0);
		for(int z = NUMGRIDZ/2.0 - 25; z <= NUMGRIDZ/2.0 + 25; z+=5)
			for(int x = NUMGRIDX/2.0 - 35; x <= NUMGRIDX/2.0 + 15; x+=5)
			{
				int height = volcano->getHeight(x-1, z-1) + 2;
				type[IX(x, height, z)] = type0[IX(x, height, z)] = FLUIDIN;
				fillParticleInGrid(x, height, z);
				Vy[IX(x, height, z)] = Vy[IX(x, height+1, z)] = rand() % 30 + 30;
				Vx[IX(x, height, z)] = Vx[IX(x, height+1, z)] = rand() % 10 - 5;
				Vz[IX(x, height, z)] = Vz[IX(x, height+1, z)] = rand() % 10 - 5;
			}
	}
	else if(SCENETYPE.compare("GROUND") == 0)
	{
		srand(0);
		for(int z = 1; z <= NUMGRIDZ; z+=20)
			for(int x = 1; x <= NUMGRIDX; x+=20)
			{
				int xx = x + (rand() % 30 - 15);
				int zz = z + (rand() % 30 - 15);
				xx = std::max(1, std::min(xx, NUMGRIDX));
				zz = std::max(1, std::min(zz, NUMGRIDZ));
				// int height = NUMGRIDY/2.0-10;
				int height = ground->getHeight(xx-1, zz-1) + 2;
				type[IX(xx, height, zz)] = type0[IX(xx, height, zz)] = FLUIDIN;
				fillParticleInGrid(xx, height, zz);
				Vy[IX(xx, height, zz)] = Vy[IX(xx, height+1, zz)] = (rand() % 3 + 1) / 10.0;
			}
	}

	if(FLOWINTYPE.compare("TOP") == 0)
	{
		// fulid come from top 
		for(int z = NUMGRIDZ/2.0-1; z <= NUMGRIDZ/2.0+1; z++)
			for(int x = NUMGRIDX/2.0-1; x <= NUMGRIDX/2.0+1; x++)
			{
				type[IX(x, NUMGRIDY+1, z)] = type0[IX(x, NUMGRIDY, z)] = FLUIDIN;
				fillParticleInGrid(x, NUMGRIDY, z);
				Vy[IX(x, NUMGRIDY, z)] = Vy[IX(x, NUMGRIDY+1, z)] = -2;
			}
	}
	else if(FLOWINTYPE.compare("SIDE") == 0)
	{ 	
		// fluid come from side
		for(int y = NUMGRIDY/2.0 - 5; y <= NUMGRIDY/2.0 + 5; y++)
			for(int x = NUMGRIDX/2.0 - 5; x <= NUMGRIDX/2.0 + 5; x++)
			{
				type[IX(x, y, 0)] = type0[IX(x, y, 0)] = FLUIDIN;
				fillParticleInGrid(x, y, 1);
				Vz[IX(x, y, 0)] = Vz[IX(x, y, 1)] = 2;
			}
	}
	else if(FLOWINTYPE.compare("BOTTOM") == 0)
	{ 
		// fluid come from bottom
		for(int z = NUMGRIDZ/2.0 - 1; z <= NUMGRIDZ/2.0 + 1; z++)
			for(int x = NUMGRIDX/2.0 - 1; x <= NUMGRIDX/2.0 + 1; x++)
			{
				type[IX(x, 0, z)] = type0[IX(x, 1, z)] = FLUIDIN;
				fillParticleInGrid(x, 1, z);
				Vy[IX(x, 0, z)] = Vy[IX(x, 1, z)] = 5;
			}
	}
}

void FluidCube3D::simulate()
{
	clock_t start, current = clock();

	bool draw = calculateTimeStep();
	
	addFlowIn();
	LOG << "addFlowIn "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	updateParticles();
	LOG << "updateParticles "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	updateGrid();
	LOG << "updateGrid "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	set_bnd();
	LOG << "set_bnd "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	vel_step();
	LOG << "vel_step "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	clock_t simTime = clock() - start;
	report(simTime);

	if(draw)
	{	
		if(!CREATEBLOBBY)
			render();
		else if(iteration % blobbyFrameStride == 0)
			createBlobbySurface();
	}
}

void FluidCube3D::vel_step()
{
	if(fluidNum == 0)
		return; 
	
	clock_t start, current = clock();

	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	SWAP(Vz0, Vz);
	advectVelosity();
	set_bnd();
	LOG << "advectVelosity "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	addForce();
	set_bnd();
	LOG << "addForce "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	projectVelosity();
	LOG << "projectVelosity "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	extrapolate();
	set_bnd();
	LOG << "extrapolate "<< (clock() - current) / 1e6 << "s\n";

	//errorRemove();
}

void FluidCube3D::addForce()
{
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
				if(type[IX(x, y, z)] == FLUID)
					Vy[IX(x, y, z)] -= dt * GRAVITY;
}

void FluidCube3D::diffuseVelosity()
{
	diffuse(1, Vx0, Vx, VISCOSITY);
	diffuse(2, Vy0, Vy, VISCOSITY);
	diffuse(3, Vz0, Vz, VISCOSITY);
}

void FluidCube3D::advectVelosity()
{
	//BFECC
	/*
	advect(1, Vx0, fai_b, true);
	advect(1, fai_b, fai_f, false);
	for(int i = 0; i < NUMGRID; i++)
		fai_b[i] = Vx0[i] + (Vx0[i] - fai_f[i]) * 0.5;
	advect(1, fai_b, Vx, true);

	advect(2, Vy0, fai_b, true);
	advect(2, fai_b, fai_f, false);
	for(int i = 0; i < NUMGRID; i++)
		fai_b[i] = Vy0[i] + (Vy0[i] - fai_f[i]) * 0.5;
	advect(2, fai_b, Vy, true);
	*/

	advect(1, Vx0, Vx, true);
	advect(2, Vy0, Vy, true);
	advect(3, Vz0, Vz, true);
}

void FluidCube3D::projectVelosity()
{

	max_vx = 0;
	max_vy = 0;
	max_vz = 0;

#ifdef GAUSS_SEIDEL
	//Gauss_Seidel
	float *p = fai_f;

	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
			{
				if(type[IX(x, y, z)] != FLUID)
					continue;

				div[IX(x, y, z)] = -h * (Vx[IX(x+1,y,z)]-Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y,z)]
										+ Vz[IX(x,y,z+1)]-Vz[IX(x,y,z)]);
				p[IX(x, y, z)] = 0;
			}
	
	for(int k = 0; k < ITERATION; k++)
	{
		for(int z = 1; z <= NUMGRIDY; z++)
			for(int y = 1; y <= NUMGRIDY; y++)
				for(int x = 1; x <= NUMGRIDX; x++)
					if(type[IX(x, y, z)] == FLUID)
						p[IX(x, y, z)] = (div[IX(x,y,z)] + p[IX(x-1,y,z)] + p[IX(x+1,y,z)] + p[IX(x,y-1,z)] + p[IX(x,y+1,z)]
										  + p[IX(x,y,z-1)] + p[IX(x,y,z+1)]) / 6;
	}

#pragma omp parallel for
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
			{
				if(type[IX(x, y, z)] != FLUID)
					continue;

				double p1, p2;
				p2 = p[IX(x, y, z)];
				//Vx
				if(type[IX(x-1, y, z)] == AIR)
					p1 = 0;
				else if(type[IX(x-1, y, z)] == FLUID)
					p1 = p[IX(x-1, y, z)];
				else
					p1 = p[IX(x, y, z)];
				Vx[IX(x, y, z)] -= (p2 - p1) * hi;
			
				//Vy
				if(type[IX(x, y-1, z)] == AIR)
					p1 = 0;
				else if(type[IX(x, y-1, z)] == FLUID)
					p1 = p[IX(x, y-1, z)];
				else
					p1 = p[IX(x, y, z)];
				Vy[IX(x, y, z)] -= (p2 - p1) * hi;

				//Vz
				if(type[IX(x, y, z-1)] == AIR)
					p1 = 0;
				else if(type[IX(x, y, z-1)] == FLUID)
					p1 = p[IX(x, y, z-1)];
				else
					p1 = p[IX(x, y, z)];
				Vz[IX(x, y, z)] -= (p2 - p1) * hi;

				if(fabsf(Vx[IX(x, y, z)]) > max_vx)
					max_vx = fabsf(Vx[IX(x, y, z)]);
				if(fabsf(Vy[IX(x, y, z)]) > max_vy)
					max_vy = fabsf(Vy[IX(x, y, z)]);
				if(fabsf(Vz[IX(x, y, z)]) > max_vz)
					max_vz = fabsf(Vz[IX(x, y, z)]);
			}
#else
	//Conjugate Gradient
	/*
	//check div before project
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
			{
				if(type[IX(x, y, z)] != FLUID)
					continue;
				div[IX(x, y, z)] = (Vx[IX(x+1,y,z)]-Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y,z)]
									+ Vz[IX(x,y,z+1)]-Vz[IX(x,y,z)]);
			}
	output(div);
	*/
	Eigen::VectorXd b(fluidNum);
	int index = 0;

// #pragma omp parallel for (cannot be parallel)
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
			{
				if(type[IX(x, y, z)] != FLUID)
					continue;

				b[index++] =  -h * (Vx[IX(x+1,y,z)]-Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y,z)]
											+ Vz[IX(x,y,z+1)]-Vz[IX(x,y,z)]);
			}
	
	p.resize(fluidNum);
	p = solver.solve(b);

	max_p = -9999;
	for(int i = 0; i < fluidNum; i++)
		if(p[i] > max_p)
			max_p = p[i];

#pragma omp parallel for
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
			{
				if(type[IX(x, y, z)] != FLUID)
					continue;

				float p1, p2;
				p2 = p[pos2index[IX(x, y, z)]];
				//Vx
				if(type[IX(x-1, y, z)] == AIR)
					p1 = 0;
				else if(type[IX(x-1, y, z)] == FLUID)
					p1 = p[pos2index[IX(x-1, y, z)]];
				else
					p1 = p2;
				Vx[IX(x, y, z)] -= (p2 - p1) * hi;
			
				//Vy
				if(type[IX(x, y-1, z)] == AIR)
					p1 = 0;
				else if(type[IX(x, y-1, z)] == FLUID)
					p1 = p[pos2index[IX(x, y-1, z)]];
				else
					p1 = p2;
				Vy[IX(x, y, z)] -= (p2 - p1) * hi;

				//Vz
				if(type[IX(x, y, z-1)] == AIR)
					p1 = 0;
				else if(type[IX(x, y, z-1)] == FLUID)
					p1 = p[pos2index[IX(x, y, z-1)]];
				else
					p1 = p2;
				Vz[IX(x, y, z)] -= (p2 - p1) * hi;

				if(fabsf(Vx[IX(x, y, z)]) > max_vx)
					max_vx = fabsf(Vx[IX(x, y, z)]);
				if(fabsf(Vy[IX(x, y, z)]) > max_vy)
					max_vy = fabsf(Vy[IX(x, y, z)]);
				if(fabsf(Vz[IX(x, y, z)]) > max_vz)
					max_vz = fabsf(Vz[IX(x, y, z)]);
			}
	
#endif

}

void FluidCube3D::diffuse(int b, float *u0, float *u, float diffusion)
{
	float a = dt * diffusion / h2;

	//Gauss Seidel relexation
	//in this way, the initcial value for u may be important 
	/*
	for(int k = 0; k < ITERATION; k++)
	{
		for(int z = 1; z <= NUMGRIDZ; z++)
			for(int y = 1; y <= NUMGRIDY; y++)
				for(int x = 1; x <= NUMGRIDX; x++)
					if(type[IX(x, y, z)] == FLUID)
					{	
						int fnum = 0;
						float v = 0;
						for(int i = 0; i < 6; i++)
						{
							int xx = x + dir[i][0];
							int yy = y + dir[i][1];
							int zz = z + dir[i][2];
							if(type[IX(xx, yy, zz)] == FLUID)
							{
								v += u[IX(xx, yy, zz)];
								fnum ++;
							}
						}
						u[IX(x, y, z)] = (u0[IX(x, y, z)] + a * v) / (1+fnum*a);
					}
	}
	*/

	//can also try unstable way
#pragma omp parallel for
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
				if(type[IX(x, y, z)] == FLUID)
				{
					u[IX(x, y, z)] = u0[IX(x, y, z)];
					for(int i = 0; i < 6; i++)
					{
						int xx = x + dir[i][0];
						int yy = y + dir[i][1];
						int zz = z + dir[i][2];
						if(type[IX(xx, yy, zz)] == FLUID)
						{
							u[IX(x, y, z)] += a * (u0[IX(xx, yy, zz)] - u0[IX(x, y, z)]);
						}	
					}
				}
}

void FluidCube3D::advect(int b, float *u0, float *u,  bool backward)
{
#pragma omp parallel for
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
				if(type[IX(x, y, z)] == FLUID)
				{
					Pos3D pos = traceParticle(b, x, y, z, backward);
					u[IX(x, y, z)] = getVelosity(b, pos.x, pos.y, pos.z, u0);
				}
}

void FluidCube3D::set_bnd()
{
	//try to use free-slip condition
	/*
	for(int y = 1; y <= NUMGRIDY; y++)
	{
		Vx[IX(1, y)] = -Vx[IX(2, y)];
		Vx[IX(NUMGRIDX+1, y)] = -Vx[IX(NUMGRIDX, y)];
	}
	for(int x = 1; x <= NUMGRIDX; x++)
	{
		Vy[IX(x, 1)] = -Vy[IX(x, 2)];
		Vy[IX(x, NUMGRIDY+1)] = -Vy[IX(x, NUMGRIDY)];
	}
	*/
	if(fluidNum == 0)
		return; 

#pragma omp parallel for
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
			{
				if(type[IX(x, y, z)] != FLUID)
					continue;

				if(type[IX(x-1, y, z)] == SOLID)
				{
					Vx[IX(x, y, z)] = 0;
					Vy[IX(x-1, y, z)] = Vy[IX(x, y, z)];
					Vz[IX(x-1, y, z)] = Vz[IX(x, y, z)];
				}
				if(type[IX(x, y-1, z)] == SOLID)
				{
					Vy[IX(x, y, z)] = 0;
					Vx[IX(x, y-1, z)] = Vx[IX(x, y, z)];
					Vz[IX(x, y-1, z)] = Vz[IX(x, y, z)];
				}
				if(type[IX(x, y, z-1)] == SOLID)
				{
					Vz[IX(x, y, z)] = 0;
					Vx[IX(x, y, z-1)] = Vx[IX(x, y, z)];
					Vy[IX(x, y, z-1)] = Vy[IX(x, y, z)];
				}
				if(type[IX(x+1, y, z)] == SOLID)
				{
					Vx[IX(x+1, y, z)] = 0;
					Vy[IX(x+1, y, z)] = Vy[IX(x, y, z)];
					Vz[IX(x+1, y, z)] = Vz[IX(x, y, z)];
				}
				if(type[IX(x, y+1, z)] == SOLID)
				{
					Vy[IX(x, y+1, z)] = 0;
					Vx[IX(x, y+1, z)] = Vx[IX(x, y, z)];
					Vz[IX(x, y+1, z)] = Vz[IX(x, y, z)];
				}
				if(type[IX(x, y, z+1)] == SOLID)
				{
					Vz[IX(x, y, z+1)] = 0;
					Vx[IX(x, y, z+1)] = Vx[IX(x, y, z)];
					Vy[IX(x, y, z+1)] = Vy[IX(x, y, z)];
				}

				switch(neighAir[IX(x, y, z)])
				{
				case 0:
					//full fluid cell
					break;
				case 1: 
					//when only 1 surface-face, solve the div = 0 directly 
					{
					if(type[IX(x-1, y, z)] == AIR)
						Vx[IX(x,y,z)] = Vx[IX(x+1,y,z)] + Vy[IX(x,y+1,z)] - Vy[IX(x,y,z)] + Vz[IX(x,y,z+1)] - Vz[IX(x,y,z)];
					else if(type[IX(x+1, y, z)] == AIR)
						Vx[IX(x+1,y,z)] = Vx[IX(x,y,z)] - Vy[IX(x,y+1,z)] + Vy[IX(x,y,z)] - Vz[IX(x,y,z+1)] + Vz[IX(x,y,z)];
					else if(type[IX(x, y-1, z)] == AIR)
						Vy[IX(x,y,z)] = Vy[IX(x,y+1,z)] + Vx[IX(x+1,y,z)] - Vx[IX(x,y,z)] + Vz[IX(x,y,z+1)] - Vz[IX(x,y,z)];
					else if(type[IX(x, y+1, z)] == AIR)
						Vy[IX(x,y+1,z)] = Vy[IX(x,y,z)] - Vx[IX(x+1,y,z)] + Vx[IX(x,y,z)] - Vz[IX(x,y,z+1)] + Vz[IX(x,y,z)];
					else if(type[IX(x, y, z-1)] == AIR)
						Vz[IX(x,y,z)] = Vz[IX(x,y,z+1)] + Vx[IX(x+1,y,z)] - Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)] - Vy[IX(x,y,z)];
					else if(type[IX(x, y, z+1)] == AIR)
						Vz[IX(x,y,z+1)] = Vz[IX(x,y,z)] - Vx[IX(x+1,y,z)] + Vx[IX(x,y,z)] - Vy[IX(x,y+1,z)] + Vy[IX(x,y,z)];
					break;
					}
				case 2:
					//when 2 surface-face, two situations:
					//a. 2 surface-face oppsite, do nothing
					//b. for surface-face, copy the value from the opposite, and add half of the left 2 difference to 2 surface-face 
					//totally 12 branches
					{
					if(type[IX(x-1, y, z)] == AIR)
					{
						if(type[IX(x, y-1, z)] == AIR)
						{
							float half = 0.5 * (Vz[IX(x, y, z)] - Vz[IX(x, y, z+1)]);
							Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)] + half;
							Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)] + half;
						}
						else if(type[IX(x, y+1, z)] == AIR)
						{
							float half = 0.5 * (Vz[IX(x, y, z)] - Vz[IX(x, y, z+1)]);
							Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)] + half;
							Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)] + half;
						}
						else if(type[IX(x, y, z-1)] == AIR)
						{
							float half = 0.5 * (Vy[IX(x, y, z)] - Vy[IX(x, y+1, z)]);
							Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)] + half;
							Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)] + half;
						}
						else if(type[IX(x, y, z+1)] == AIR)
						{
							float half = 0.5 * (Vy[IX(x, y, z)] - Vy[IX(x, y+1, z)]);
							Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)] + half;
							Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)] + half;
						}
					}
					else if(type[IX(x+1, y, z)] == AIR)
					{
						if(type[IX(x, y-1, z)] == AIR)
						{
							float half = 0.5 * (Vz[IX(x, y, z)] - Vz[IX(x, y, z+1)]);
							Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)] + half;
							Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)] + half;
						}
						else if(type[IX(x, y+1, z)] == AIR)
						{
							float half = 0.5 * (Vz[IX(x, y, z)] - Vz[IX(x, y, z+1)]);
							Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)] + half;
							Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)] + half;
						}
						else if(type[IX(x, y, z-1)] == AIR)
						{
							float half = 0.5 * (Vy[IX(x, y, z)] - Vy[IX(x, y+1, z)]);
							Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)] + half;
							Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)] + half;
						}
						else if(type[IX(x, y, z+1)] == AIR)
						{
							float half = 0.5 * (Vy[IX(x, y, z)] - Vy[IX(x, y+1, z)]);
							Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)] + half;
							Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)] + half;
						}
					}
					else if(type[IX(x, y-1, z)] == AIR)
					{
						if(type[IX(x, y, z-1)] == AIR)
						{
							float half = 0.5 * (Vx[IX(x, y, z)] - Vx[IX(x+1, y, z)]);
							Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)] + half;
							Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)] + half;
						}
						else if(type[IX(x, y, z+1)] == AIR)
						{
							float half = 0.5 * (Vx[IX(x, y, z)] - Vx[IX(x+1, y, z)]);
							Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)] + half;
							Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)] + half;
						}
					}
					else if(type[IX(x, y+1, z)] == AIR)
					{
						if(type[IX(x, y, z-1)] == AIR)
						{
							float half = 0.5 * (Vx[IX(x, y, z)] - Vx[IX(x+1, y, z)]);
							Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)] + half;
							Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)] + half;
						}
						else if(type[IX(x, y, z+1)] == AIR)
						{
							float half = 0.5 * (Vx[IX(x, y, z)] - Vx[IX(x+1, y, z)]);
							Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)] + half;
							Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)] + half;
						}
					}
					break;
					}
				case 3:
					//when 3 surface-face, two situations:
					//a. 3 not oppsite each other, copy the value from the oppsite
					//b. for 1 surface-face not oppsite, solve it as the case 1
					{
					bool tx0 = (type[IX(x-1, y, z)] == AIR && type[IX(x+1, y, z)] != AIR);
					bool tx1 = (type[IX(x+1, y, z)] == AIR && type[IX(x-1, y, z)] != AIR);
					bool ty0 = (type[IX(x, y-1, z)] == AIR && type[IX(x, y+1, z)] != AIR);
					bool ty1 = (type[IX(x, y+1, z)] == AIR && type[IX(x, y-1, z)] != AIR);
					bool tz0 = (type[IX(x, y, z-1)] == AIR && type[IX(x, y, z+1)] != AIR);
					bool tz1 = (type[IX(x, y, z+1)] == AIR && type[IX(x, y, z-1)] != AIR);
					if(tx0 && ty0 && tz0)
					{
						Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)];
						Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)];
						Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)];
					}
					else if(tx1 && ty0 && tz0)
					{
						Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)];
						Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)];
						Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)];
					}
					else if(tx0 && ty1 && tz0)
					{
						Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)];
						Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)];
						Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)];
					}
					else if(tx1 && ty1 && tz0)
					{
						Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)];
						Vy[IX(x+1, y, z)] = Vy[IX(x, y, z)];
						Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)];
					}
					else if(tx0 && ty0 && tz1)
					{
						Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)];
						Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)];
						Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)];
					}
					else if(tx1 && ty0 && tz1)
					{
						Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)];
						Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)];
						Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)];
					}
					else if(tx0 && ty1 && tz1)
					{
						Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)];
						Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)];
						Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)];
					}
					else if(tx1 && ty1 && tz1)
					{
						Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)];
						Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)];
						Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)];
					}
					else if(tx0 && !(ty0||ty1) )
						Vx[IX(x,y,z)] = Vx[IX(x+1,y,z)] + Vy[IX(x,y+1,z)] - Vy[IX(x,y,z)] + Vz[IX(x,y,z+1)] - Vz[IX(x,y,z)];
					else if(tx1 && !(ty0||ty1) )
						Vx[IX(x+1,y,z)] = Vx[IX(x,y,z)] - Vy[IX(x,y+1,z)] + Vy[IX(x,y,z)] - Vz[IX(x,y,z+1)] + Vz[IX(x,y,z)];
					else if(ty0 && !(tx0||tx1) )
						Vy[IX(x,y,z)] = Vy[IX(x,y+1,z)] + Vx[IX(x+1,y,z)] - Vx[IX(x,y,z)] + Vz[IX(x,y,z+1)] - Vz[IX(x,y,z)];
					else if(ty1 && !(tx0||tx1) )
						Vy[IX(x,y+1,z)] = Vy[IX(x,y,z)] - Vx[IX(x+1,y,z)] + Vx[IX(x,y,z)] - Vz[IX(x,y,z+1)] + Vz[IX(x,y,z)];
					else if(tz0 && !(tx0||tx1) )
						Vz[IX(x,y,z)] = Vz[IX(x,y,z+1)] + Vx[IX(x+1,y,z)] - Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)] - Vy[IX(x,y,z)];
					else if(tz1 && !(tx0||tx1) )
						Vz[IX(x,y,z+1)] = Vz[IX(x,y,z)] - Vx[IX(x+1,y,z)] + Vx[IX(x,y,z)] - Vy[IX(x,y+1,z)] + Vy[IX(x,y,z)];
					break;
					}
				case 4:
					//when 4 surface-face, two situations:
					//a. each 2 opposite, add the quater of the div to each surface-face
					//b. for 2 opposite non-surface face, solve them as the case 2
					{
					if(type[IX(x+1, y, z)] != AIR)
					{
						if(type[IX(x, y+1, z)] != AIR)
						{
							float half = 0.5 * (Vz[IX(x, y, z)] - Vz[IX(x, y, z+1)]);
							Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)] + half;
							Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)] + half;
						}
						else if(type[IX(x, y-1, z)] != AIR)
						{
							float half = 0.5 * (Vz[IX(x, y, z)] - Vz[IX(x, y, z+1)]);
							Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)] + half;
							Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)] + half;
						}
						else if(type[IX(x, y, z+1)] != AIR)
						{
							float half = 0.5 * (Vy[IX(x, y, z)] - Vy[IX(x, y+1, z)]);
							Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)] + half;
							Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)] + half;
						}
						else if(type[IX(x, y, z-1)] != AIR)
						{
							float half = 0.5 * (Vy[IX(x, y, z)] - Vy[IX(x, y+1, z)]);
							Vx[IX(x, y, z)] = Vx[IX(x+1, y, z)] + half;
							Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)] + half;
						}
					}
					else if(type[IX(x-1, y, z)] != AIR)
					{
						if(type[IX(x, y+1, z)] != AIR)
						{
							float half = 0.5 * (Vz[IX(x, y, z)] - Vz[IX(x, y, z+1)]);
							Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)] + half;
							Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)] + half;
						}
						else if(type[IX(x, y-1, z)] != AIR)
						{
							float half = 0.5 * (Vz[IX(x, y, z)] - Vz[IX(x, y, z+1)]);
							Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)] + half;
							Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)] + half;
						}
						else if(type[IX(x, y, z+1)] != AIR)
						{
							float half = 0.5 * (Vy[IX(x, y, z)] - Vy[IX(x, y+1, z)]);
							Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)] + half;
							Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)] + half;
						}
						else if(type[IX(x, y, z-1)] != AIR)
						{
							float half = 0.5 * (Vy[IX(x, y, z)] - Vy[IX(x, y+1, z)]);
							Vx[IX(x+1, y, z)] = Vx[IX(x, y, z)] + half;
							Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)] + half;
						}
					}
					else if(type[IX(x, y+1, z)] != AIR)
					{
						if(type[IX(x, y, z+1)] != AIR)
						{
							float half = 0.5 * (Vx[IX(x, y, z)] - Vx[IX(x+1, y, z)]);
							Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)] + half;
							Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)] + half;
						}
						else if(type[IX(x, y, z-1)] != AIR)
						{
							float half = 0.5 * (Vx[IX(x, y, z)] - Vx[IX(x+1, y, z)]);
							Vy[IX(x, y, z)] = Vy[IX(x, y+1, z)] + half;
							Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)] + half;
						}
					}
					else if(type[IX(x, y-1, z)] != AIR)
					{
						if(type[IX(x, y, z+1)] != AIR)
						{
							float half = 0.5 * (Vx[IX(x, y, z)] - Vx[IX(x+1, y, z)]);
							Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)] + half;
							Vz[IX(x, y, z)] = Vz[IX(x, y, z+1)] + half;
						}
						else if(type[IX(x, y, z-1)] != AIR)
						{
							float half = 0.5 * (Vx[IX(x, y, z)] - Vx[IX(x+1, y, z)]);
							Vy[IX(x, y+1, z)] = Vy[IX(x, y, z)] + half;
							Vz[IX(x, y, z+1)] = Vz[IX(x, y, z)] + half;
						}
					}
					else if(type[IX(x-1, y, z)] != AIR && type[IX(x+1, y, z)] != AIR)
					{
						float div4 = -(Vx[IX(x+1,y,z)]-Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y,z)]
									 + Vz[IX(x,y,z+1)]-Vz[IX(x,y,z)]) / 4;
						Vy[IX(x, y, z)] += div4;
						Vz[IX(x, y, z)] += div4;
						Vy[IX(x, y+1, z)] += div4;
						Vz[IX(x, y, z+1)] += div4;
					}
					else if(type[IX(x, y-1, z)] != AIR && type[IX(x, y+1, z)] != AIR)
					{
						float div4 = -(Vx[IX(x+1,y,z)]-Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y,z)]
									 + Vz[IX(x,y,z+1)]-Vz[IX(x,y,z)]) / 4;
						Vx[IX(x, y, z)] += div4;
						Vz[IX(x, y, z)] += div4;
						Vx[IX(x+1, y, z)] += div4;
						Vz[IX(x, y, z+1)] += div4;
					}
					else if(type[IX(x, y, z-1)] != AIR && type[IX(x, y, z+1)] != AIR)
					{
						float div4 = -(Vx[IX(x+1,y,z)]-Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y,z)]
									 + Vz[IX(x,y,z+1)]-Vz[IX(x,y,z)]) / 4;
						Vx[IX(x, y, z)] += div4;
						Vy[IX(x, y, z)] += div4;
						Vx[IX(x+1, y, z)] += div4;
						Vy[IX(x, y+1, z)] += div4;
					}
					break;
					}
				case 5:
					//when 5 surface-face, solve the one opposite non-surface face as the case 1
					{
					if(type[IX(x+1, y, z)] != AIR)
						Vx[IX(x,y,z)] = Vx[IX(x+1,y,z)] + Vy[IX(x,y+1,z)] - Vy[IX(x,y,z)] + Vz[IX(x,y,z+1)] - Vz[IX(x,y,z)];
					else if(type[IX(x-1, y, z)] != AIR)
						Vx[IX(x+1,y,z)] = Vx[IX(x,y,z)] - Vy[IX(x,y+1,z)] + Vy[IX(x,y,z)] - Vz[IX(x,y,z+1)] + Vz[IX(x,y,z)];
					else if(type[IX(x, y+1, z)] != AIR)
						Vy[IX(x,y,z)] = Vy[IX(x,y+1,z)] + Vx[IX(x+1,y,z)] - Vx[IX(x,y,z)] + Vz[IX(x,y,z+1)] - Vz[IX(x,y,z)];
					else if(type[IX(x, y-1, z)] != AIR)
						Vy[IX(x,y+1,z)] = Vy[IX(x,y,z)] - Vx[IX(x+1,y,z)] + Vx[IX(x,y,z)] - Vz[IX(x,y,z+1)] + Vz[IX(x,y,z)];
					else if(type[IX(x, y, z+1)] != AIR)
						Vz[IX(x,y,z)] = Vz[IX(x,y,z+1)] + Vx[IX(x+1,y,z)] - Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)] - Vy[IX(x,y,z)];
					else if(type[IX(x, y, z-1)] != AIR)
						Vz[IX(x,y,z+1)] = Vz[IX(x,y,z)] - Vx[IX(x+1,y,z)] + Vx[IX(x,y,z)] - Vy[IX(x,y+1,z)] + Vy[IX(x,y,z)];
					break;
					}
				case 6:
					//when 6 surface-face, do nothing
					break;
				}
			}
}

bool FluidCube3D::calculateTimeStep()
{
	iteration ++;

	if(max_v == 0)
		dt = frameTime;
	else
		dt = h / max_v;
	if(dt > frameTime)
		dt = frameTime;

	dt = frameTime;
	totalTime += dt;
	return true;

	// if(ctime + dt >= frameTime)
	// {
	// 	dt = frameTime - ctime;
	// 	totalTime += dt;
	// 	ctime = 0;
	// 	return true;
	// }
	// else
	// {
	// 	ctime += dt;
	// 	totalTime += dt;
	// 	return false;
	// }
}

void FluidCube3D::updateParticles()
{
	std::vector<Pos3D> particles_new;
	std::vector<Velo3D> velosities_new;
// #pragma omp parallel for
	for(unsigned i = 0; i < particles.size(); i++)
	{
		float x0 = particles[i].x;
		float y0 = particles[i].y;
		float z0 = particles[i].z;
		Velo3D v0 = getVelosity(x0, y0, z0, Vx, Vy, Vz);

		float x1 = x0 + dt * v0.x * hi;
		float y1 = y0 + dt * v0.y * hi;
		float z1 = z0 + dt * v0.z * hi;

		if (x1 < 1 || x1 > NUMGRIDX || y1 < 1 || y1 > NUMGRIDY || z1 < 1 || z1 > NUMGRIDZ)
			continue;

		if(type[IX(int(x1), int(y1), int(z1))] == AIR)
		{
			particles_new.push_back(Pos3D(x1, y1, z1));
			velosities_new.push_back(v0);
			continue;
		}
		else if(type[IX(int(x1), int(y1), int(z1))] == SOLID)
		{
			if(int(x1) > int(x0))
				x1 = int(x0)+0.99;
			else if(int(x1) < int(x0))
				x1 = int(x0)+0.01;
			if(int(y1) > int(y0))
				y1 = int(y0)+0.99;
			else if(int(y1) < int(y0))
				y1 = int(y0)+0.01;
			if(int(z1) > int(z0))
				z1 = int(z0)+0.99;
			else if(int(z1) < int(z0))
				z1 = int(z0)+0.01;
			particles_new.push_back(Pos3D(x1, y1, z1));
			velosities_new.push_back(v0);
			continue;
		}

		// if FLUID
		Velo3D v1 = getVelosity(x1, y1, z1, Vx, Vy, Vz);
		x1 = x0 + dt * 0.5 * (v0.x + v1.x) * hi;
		y1 = y0 + dt * 0.5 * (v0.y + v1.y) * hi;
		z1 = z0 + dt * 0.5 * (v0.z + v1.z) * hi;

		if (x1 < 1 || x1 > NUMGRIDX || y1 < 1 || y1 > NUMGRIDY || z1 < 1 || z1 > NUMGRIDZ)
			continue;

		if(type[IX(int(x1), int(y1), int(z1))] == SOLID)
		{
			if(int(x1) > int(x0))
				x1 = int(x0)+0.99;
			else if(int(x1) < int(x0))
				x1 = int(x0)+0.01;
			if(int(y1) > int(y0))
				y1 = int(y0)+0.99;
			else if(int(y1) < int(y0))
				y1 = int(y0)+0.01;
			if(int(z1) > int(z0))
				z1 = int(z0)+0.99;
			else if(int(z1) < int(z0))
				z1 = int(z0)+0.01;
		}
		
		particles_new.push_back(Pos3D(x1, y1, z1));
		velosities_new.push_back(Velo3D((v0.x+v1.x)*0.5, (v0.y+v1.y)*0.5, (v0.z+v1.z)*0.5));
	}
	particles = particles_new;
	velosities = velosities_new;
}

void FluidCube3D::updateGrid()
{
	clock_t start, current = clock();

	//swap 
	GridType *tmp = type;
	type = type0;
	type0 = tmp;

	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
			{
				if(type0[IX(x, y, z)] == AIR || type0[IX(x, y, z)] == FLUID)
				{
					type[IX(x, y, z)] = AIR;
					invertedList[IX(x, y, z)]->clear();
				}
				layer[IX(x, y, z)] = -1;
			}

	LOG << "seg 0 "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	for(unsigned i = 0; i < particles.size(); i++)
	{
		int x = int(particles[i].x);
		int y = int(particles[i].y);
		int z = int(particles[i].z);
		type[IX(x, y, z)] = FLUID;
		layer[IX(x, y, z)] = 0;
		invertedList[IX(x, y, z)]->push_back(i);
	}

	LOG << "seg 1 "<< (clock() - current) / 1e6 << "s\n";
	current = clock();
	
	fluidNum = 0;
	for(int i = 0; i < NUMGRID; i++)
		pos2index[i] = -1;

	LOG << "seg 2 "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

// #pragma omp parallel for (this loop cannot be paralleled)
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
				if(type[IX(x, y, z)] == FLUID)
					pos2index[IX(x, y, z)] = fluidNum++;

	LOG << "seg 3 "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	if(fluidNum == 0)
		return; 

#pragma omp parallel for
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
			{
				if(type[IX(x, y, z)] == FLUID)
				{
					neighNoneSolid[IX(x, y, z)] = 0;
					neighAir[IX(x, y, z)] = 0;
					for(int i = 0; i < 6; i++)
					{
						int xx = x+dir[i][0];
						int yy = y+dir[i][1];
						int zz = z+dir[i][2];
						neighbor[IX(x, y, z)][i] = pos2index[IX(xx, yy, zz)];
						//if(type[IX(xx, yy, zz)] != SOLID)
						if(type[IX(xx, yy, zz)] == AIR || type[IX(xx, yy, zz)] == FLUID)
							neighNoneSolid[IX(x, y, z)] ++;
						if(type[IX(xx, yy, zz)] == AIR)
							neighAir[IX(x, y, z)] ++;
					}
				}

			//set velosity for new fluid cell
			float downRate = 0.7;
			if(type0[IX(x, y, z)] == AIR && type[IX(x, y, z)] == FLUID)
			{
				float dist = 100;
				int pid = -1;
				//Vx
				std::vector<int> *list = invertedList[IX(x, y, z)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					float x0 = particles[list->at(i)].x;
					//float y0 = particles[list->at(i)].y;
					//float z0 = particles[list->at(i)].z;
					/*
					if(DISTANCE(x0, y0, x, y+0.5)< dist)
					{
						dist = DISTANCE(x0, y0, x, y+0.5);
						pid = list->at(i);
					}
					*/
					if(x0 - x < dist)
					{
						dist = x0 - x;
						pid = list->at(i);
					}
				}
				list = invertedList[IX(x-1, y, z)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					float x0 = particles[list->at(i)].x;
					//float y0 = particles[list->at(i)].y;
					/*
					if(DISTANCE(x0, y0, x, y+0.5) < dist)
					{
						dist = DISTANCE(x0, y0, x, y+0.5);
						pid = list->at(i);
					}
					*/
					if(x - x0 < dist)
					{
						dist = x - x0;
						pid = list->at(i);
					}
				}
				Vx[IX(x, y, z)] = velosities[pid].x * downRate;
				//Vx[IX(x, y, z)] = getVelosity(1, particles[pid].x, particles[pid].y, particles[pid].z, Vx);

				//Vy
				dist = 100;
				list = invertedList[IX(x, y, z)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					//float x0 = particles[list->at(i)].x;
					float y0 = particles[list->at(i)].y;
					/*
					if(DISTANCE(x0, y0, x+0.5, y) < dist)
					{
						dist = DISTANCE(x0, y0, x+0.5, y);
						pid = list->at(i);
					}
					*/
					if(y0 - y < dist)
					{
						dist = y0 - y;
						pid = list->at(i);
					}
				}
				list = invertedList[IX(x, y-1, z)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					//float x0 = particles[list->at(i)].x;
					float y0 = particles[list->at(i)].y;
					/*
					if(DISTANCE(x0, y0, x+0.5, y) < dist)
					{
						dist = DISTANCE(x0, y0, x+0.5, y);
						pid = list->at(i);
					}
					*/
					if(y - y0 < dist)
					{
						dist = y - y0;
						pid = list->at(i);
					}
				}
				Vy[IX(x, y, z)] = velosities[pid].y * downRate;
				//Vy[IX(x, y, z)] = getVelosity(2, particles[pid].x, particles[pid].y, particles[pid].z, Vy);

				//Vz
				dist = 100;
				list = invertedList[IX(x, y, z)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					//float x0 = particles[list->at(i)].x;
					//float y0 = particles[list->at(i)].y;
					float z0 = particles[list->at(i)].z;
					/*
					if(DISTANCE(x0, y0, x+0.5, y) < dist)
					{
						dist = DISTANCE(x0, y0, x+0.5, y);
						pid = list->at(i);
					}
					*/
					if(z0 - z < dist)
					{
						dist = z0 - z;
						pid = list->at(i);
					}
				}
				list = invertedList[IX(x, y, z-1)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					//float x0 = particles[list->at(i)].x;
					//float y0 = particles[list->at(i)].y;
					float z0 = particles[list->at(i)].z;
					/*
					if(DISTANCE(x0, y0, x+0.5, y) < dist)
					{
						dist = DISTANCE(x0, y0, x+0.5, y);
						pid = list->at(i);
					}
					*/
					if(z - z0 < dist)
					{
						dist = z - z0;
						pid = list->at(i);
					}
				}
				Vz[IX(x, y, z)] = velosities[pid].z * downRate;
				//Vz[IX(x, y, z)] = getVelosity(3, particles[pid].x, particles[pid].y, particles[pid].z, Vz);
			}

		}

	LOG << "seg 4 "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	//init Matrix
	A = Eigen::SparseMatrix<double>(fluidNum, fluidNum);         // default is column major
	A.reserve(Eigen::VectorXi::Constant(fluidNum, 7));
	int index = 0;
// #pragma omp parallel for (cannot be parallel)
	for(int z = 1; z <= NUMGRIDZ; z++)
		for(int y = 1; y <= NUMGRIDY; y++)
			for(int x = 1; x <= NUMGRIDX; x++)
			{
				if(type[IX(x, y, z)] != FLUID)
					continue;

				for(int i = 0; i < 3; i++)
				{
					int neighid = neighbor[IX(x, y, z)][i];
					if(neighid != -1)
					{
						A.insert(index, neighid) = -1;
					}
				}
				A.insert(index, index) = neighNoneSolid[IX(x, y, z)];
				for(int i = 3; i < 6; i++)
				{
					int neighid = neighbor[IX(x, y, z)][i];
					if(neighid != -1)
					{
						A.insert(index, neighid) = -1;
					}
				}
				index ++;
			}
	LOG << "seg 5 "<< (clock() - current) / 1e6 << "s\n";
	current = clock();

	A.makeCompressed();
	solver.compute(A);
	LOG << "seg 6 "<< (clock() - current) / 1e6 << "s\n";
	current = clock();
}

void FluidCube3D::extrapolate()
{
	if(max_vx > max_vy)
		if(max_vx > max_vz)
			max_v = max_vx;
		else
			max_v = max_vz;
	else
		if(max_vy > max_vz)
			max_v = max_vy;
		else
			max_v = max_vz;

	//make sure all the fluid cell with layer = 0, others = -1
	int iter = int(max_v * dt * hi) + 2;
	REPORT(iter);

#pragma omp parallel for
	for(int k = 1; k <= 3; k++)
	{
		int start[3], end[3], du[3];
		int maxu[3] = {NUMGRIDX, NUMGRIDY, NUMGRIDZ};
		for(int i = 0; i < 3; i++)
		{
			du[i] = (rand() % 2 == 0)? 1 : -1;
			if(du[i] > 0)
			{
				start[i] = 0;
				end[i] = maxu[i];
			}
			else
			{
				start[i] = maxu[i]-1;
				end[i] = -1;
			}
		}
		for(int z = start[2]; z != end[2]; z += du[2])
			for(int y = start[1]; y != end[1]; y += du[1])
				for(int x = start[0]; x != end[0]; x += du[0])
				{
					if(type[IX(x, y, z)] == AIR && layer[IX(x, y, z)] == -1)
					{
						int nei = 0;
						float velo[] = {0, 0, 0};
						for(int j = 0; j < 6; j++)
						{
							int xx = x + dir[j][0];
							int yy = y + dir[j][1];
							int zz = z + dir[j][2];
							if(layer[IX(xx, yy, zz)] == k-1)
							{
								nei ++;
								velo[0] += Vx[IX(xx, yy, zz)];
								velo[1] += Vy[IX(xx, yy, zz)];
								velo[2] += Vz[IX(xx, yy, zz)];
							}
						}
						if(nei > 0)
						{
							Vx[IX(x, y, z)] = velo[0] / nei;
							Vy[IX(x, y, z)] = velo[1] / nei;
							Vz[IX(x, y, z)] = velo[2] / nei;
							layer[IX(x, y, z)] = k;
						}
					}
				}
	}
}

void FluidCube3D::errorRemove()
{
	double eps = 1e-12;

	for(int i = 0; i < NUMGRID; i++)
	{
		if(fabs(Vx[i]) < eps)
			Vx[i] = 0;
		if(fabs(Vy[i]) < eps)
			Vy[i] = 0;
		if(fabs(Vz[i]) < eps)
			Vz[i] = 0;
	}
}

float FluidCube3D::getVelosity(int index, float x, float y, float z, float *u)
{
	switch(index)
	{
	case 1:
		y -= 0.5;
		z -= 0.5;
		break;
	case 2:
		x -= 0.5;
		z -= 0.5;
		break;
	case 3:
		x -= 0.5;
		y -= 0.5;
		break;
	}

	// if(x < 0 || x >= NUMGRIDX+1 || y < 0 || y >= NUMGRIDY+1 || z < 0 || z >= NUMGRIDZ+1)
	// {
	// 	std::cout<<"Get velosity out of bound"<<std::endl;
	// 	REPORT(x);
	// 	REPORT(y);
	// 	REPORT(z);
	// }
	if(x < 0)
		x = 0;
	if(x >= NUMGRIDX+1)
		x = NUMGRIDX+0.999;
	if(y < 0)
		y = 0;
	if(y >= NUMGRIDY+1)
		y = NUMGRIDY+0.999;
	if(z < 0)
		z = 0;
	if(z >= NUMGRIDZ+1)
		z = NUMGRIDZ+0.999;


	int i0 = int(x), i1 = i0 + 1;
	int j0 = int(y), j1 = j0 + 1;
	int k0 = int(z), k1 = k0 + 1;
	float s1 = x - i0, s0 = 1 - s1;
	float t1 = y - j0, t0 = 1 - t1;
	float r1 = z - k0, r0 = 1 - r1;

	float a1 = r0 * u[IX(i0, j0, k0)] + r1 * u[IX(i0, j0, k1)];
	float a2 = r0 * u[IX(i0, j1, k0)] + r1 * u[IX(i0, j1, k1)];
	float b1 = r0 * u[IX(i1, j0, k0)] + r1 * u[IX(i1, j0, k1)];
	float b2 = r0 * u[IX(i1, j1, k0)] + r1 * u[IX(i1, j1, k1)];

	float c1 = t0 * a1 + t1 * a2;
	float c2 = t0 * b1 + t1 * b2;

	return s0 * c1 + s1 * c2;
}

Velo3D FluidCube3D::getVelosity(float x, float y, float z, float *vx, float *vy, float *vz)
{
	return Velo3D(getVelosity(1,x,y,z,vx), getVelosity(2,x,y,z,vy), getVelosity(3,x,y,z,vz));
}

Pos3D FluidCube3D::traceParticle(int index, int x, int y, int z, bool backward)
{
	float x0, y0, z0;
	switch(index)
	{
	case 1:
		x0 = x;
		y0 = y + 0.5;
		z0 = z + 0.5;
		break;
	case 2:
		x0 = x + 0.5;
		y0 = y;
		z0 = z + 0.5;
		break;
	case 3:
		x0 = x + 0.5;
		y0 = y + 0.5;
		z0 = z;
		break;
	}

	Velo3D v0 = getVelosity(x0, y0, z0, Vx0, Vy0, Vz0);
	float t = (backward)? -dt : dt;
	//
	//Pos3D p = Pos3D(x0 + v0.x*t*hi, y0 + v0.y*t*hi);
	Velo3D v1 = getVelosity(x0 + v0.x*t*hi, y0 + v0.y*t*hi, z0 + v0.z*t*hi, Vx0, Vy0, Vz0);
	return Pos3D(x0 + 0.5*t*(v0.x+v1.x)*hi, y0 + 0.5*t*(v0.y+v1.y)*hi, z0 + 0.5*t*(v0.z+v1.z)*hi);
}

void FluidCube3D::fillParticleInGrid(int x, int y, int z)
{
	int nump = PARTICLEPERGRID;
	int sample = 10000;
	float subSize = 1.0 / nump;
	for(int i = 0; i < nump; i++)
		for(int j = 0; j < nump; j++)
			for(int k = 0; k < nump; k++)
			{
				float x0 = subSize * (rand()%sample) / sample + i * subSize;
				float y0 = subSize * (rand()%sample) / sample + j * subSize;
				float z0 = subSize * (rand()%sample) / sample + k * subSize;
				//float x0 = i * step;
				//float y0 = j * step;
				//float z0 = k * step;
				particles.push_back(Pos3D(x+x0, y+y0, z+z0) );
				velosities.push_back(Velo3D(0, 0, 0));
			}
}

void FluidCube3D::report(clock_t simTime)
{
	REPORT(iteration);
	REPORT(max_vx);
	REPORT(max_vy);
	REPORT(max_vz);
	REPORT(max_p);
	float fluidShrink = 1.0 * fluidNum / originFluid;
	REPORT(fluidShrink);
	REPORT(totalTime);
	std::cout<<"Simulation finished in "<<simTime<<"ms"<<std::endl;
	LOGSEG;
}

void FluidCube3D::output(float *u)
{
	if(!DEBUGPRINT)
		return;  
	for(int z = 10; z <= 15; z++)
		for(int y = 10; y <= 15; y++)
			for(int x = 5; x <= 10; x++)
			{
				std::cout<<u[IX(x, y, z)]<<' ';
				if(x == 10)
					std::cout<<std::endl;
			}
	LOGSEG;
}

void FluidCube3D::render()
{
	glMatrixMode(GL_MODELVIEW);
	// Reset transformations
	glLoadIdentity();
	glScalef(0.02f, 0.02f, 0.02f);

	// Set the camera
	gluLookAt(	px, py, pz,
				0, 0, 0,
				0.0f, 1.0f, 0.0f);

	//clear the screen to a desired color in range [0-1] RGBA
	glClearColor(0.5, 0.5, 0.5, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//enable blending for translucency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	int LENGTH = NUMGRIDX * GRIDSIZE;
	//calculate divergence
	if(RENDERTYPE == DIVERGENCE)
	{
#pragma omp parallel for
		for(int z = 1; z <= NUMGRIDZ; z++)
			for(int y = 1; y <= NUMGRIDY; y++)
				for(int x = 1; x <= NUMGRIDX; x++)
				{
					if(type[IX(x, y, z)] != FLUID)
						continue;
					div[IX(x, y, z)] = (Vx[IX(x+1,y,z)]-Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y,z)]
										+ Vz[IX(x,y,z+1)]-Vz[IX(x,y,z)]);
				}
	}

	glTranslatef(-LENGTH/2, -LENGTH/2, -LENGTH/2);

// #pragma omp parallel for (cannot make rendering parallel)
	for(int k = 0; k < NUMGRIDZ; k++)	
		for(int j = 0; j < NUMGRIDY; j++)
			for(int i = 0; i < NUMGRIDX; i++)
			{
				int x = i+1;
				int y = j+1;
				int z = k+1;
				float color;
				if(type[IX(x, y, z)] == SOLID)
					glColor4f(0, 0.25, 0, 1);
				else if(type[IX(x, y, z)] == FLUID)
				{
					switch(RENDERTYPE)
					{
					case FLUIDGRID:
						glColor3f(0, 0, 0.7);
						break;
					case PRESSURE:
						glColor3f(p[pos2index[IX(x,y,z)]]/max_p, 0, 0);
						break;
					case VELOSITYY:
						if(Vy[IX(x, y, z)] >= 0)
							glColor3f(Vy[IX(x, y, z)]/max_vy, 0, 0);
						else
							glColor3f(0, -Vy[IX(x, y, z)]/max_vy, 0);
						break;
					case VELOSITYX:
						if(Vx[IX(x, y, z)] >= 0)
							glColor3f(Vx[IX(x, y, z)]/max_vx, 0, 0);
						else
							glColor3f(0, -Vx[IX(x, y, z)]/max_vx, 0);
						break;
					case DIVERGENCE:
						if(fabs(div[IX(x, y, z)]) > 1e-5)
							glColor3f(1, 0, 0);
						else
							glColor3f(0, 0, 0.7);
						break;
					case PARTICLE:
						continue;
						break;
					default:
						glColor3f(0, 0, 0.7);
						break;
					//vorticity
					//float w = 0.5 * (Vy[IX(x+1, y)] - Vy[IX(x-1, y)]);
					//		  - 0.5 * (Vx[IX(x, y+1)] - Vx[IX(x, y-1)]);
					}
				}
				else //AIR
					continue;

				//draw cube 
				glBegin(GL_QUADS);
				//hold k
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				//hold k+1
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				//hold j
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				//hold j+1
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				//hold i
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				//hold i+1
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				glEnd();
				//if(GRIDSIZE >= 10 && type[IX(x, y)] == FLUID)
				//	draw_velo(i, j, Vx[IX(x, y)], Vy[IX(x, y)]);
			}

	//draw particles
	if(RENDERTYPE == PARTICLE)
	{
		glColor3f(0.3, 0, 0);
		glBegin(GL_POINTS);
		for(unsigned i = 0; i < particles.size(); i++)
		{
			glVertex3f((particles[i].x-1)*GRIDSIZE, (particles[i].y-1)*GRIDSIZE, (particles[i].z-1)*GRIDSIZE);
		}
		glEnd();
	}

	
	//draw box
	glColor3f(1, 1, 1);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0, LENGTH, 0);
	glVertex3f(0, LENGTH, 0);
	glVertex3f(LENGTH, LENGTH, 0);
	glVertex3f(LENGTH, LENGTH, 0);
	glVertex3f(LENGTH, 0, 0);
	glVertex3f(LENGTH, 0, 0);
	glVertex3f(0, 0, 0);

	glVertex3f(0, 0, LENGTH);
	glVertex3f(0, LENGTH, LENGTH);
	glVertex3f(0, LENGTH, LENGTH);
	glVertex3f(LENGTH, LENGTH, LENGTH);
	glVertex3f(LENGTH, LENGTH, LENGTH);
	glVertex3f(LENGTH, 0, LENGTH);
	glVertex3f(LENGTH, 0, LENGTH);
	glVertex3f(0, 0, LENGTH);

	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, LENGTH);
	glVertex3f(0, LENGTH, LENGTH);
	glVertex3f(0, LENGTH, 0);
	glVertex3f(LENGTH, LENGTH, LENGTH);
	glVertex3f(LENGTH, LENGTH, 0);
	glVertex3f(LENGTH, 0, LENGTH);
	glVertex3f(LENGTH, 0, 0);
	glEnd();

	glutSwapBuffers();
}

void FluidCube3D::createBlobby()
{
	char obj_dir_base[] = "../surface3D";
	mkdir(obj_dir_base, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	clock_t timeStamp = clock();
	sprintf(obj_dir, "%s/%ld", obj_dir_base, timeStamp);
	mkdir(obj_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	int frameSum = 1000;
	Configer::getConfiger()->getInt("Blobby", "FrameSum", frameSum);
	Configer::getConfiger()->getInt("Blobby", "FrameStride", blobbyFrameStride);
	while(iteration < frameSum)
	{
		simulate();
	}
}

void FluidCube3D::createBlobbySurface()
{
	int blobbyGridSize = 1;
	assert(Configer::getConfiger()->getInt("Blobby", "GridSize", blobbyGridSize));

	double blobbyRadius = 1.0 * blobbyGridSize / PARTICLEPERGRID;
	// assert(Configer::getConfiger()->getDouble("Blobby", "Radius", blobbyRadius));

	double blobbyH = 3 * blobbyRadius;
	double blobbyH2i = 1 / (blobbyH * blobbyH);
	double thresh = blobbyKernel( (blobbyRadius * blobbyRadius) / (blobbyH * blobbyH) );
	
	int border = 6;
	float *scalarField = new float [(NUMGRIDX+border) * (NUMGRIDY+border) * (NUMGRIDZ+border) * blobbyGridSize * blobbyGridSize * blobbyGridSize]; 

	PRINT("Start to sample the implict function on the grid...");
	clock_t start = clock();

#pragma omp parallel for
	for(int z = 0; z < (NUMGRIDZ+border) * blobbyGridSize; z++)
		for(int y = 0; y < (NUMGRIDY+border) * blobbyGridSize; y++)
			for(int x = 0; x < (NUMGRIDX+border) * blobbyGridSize; x++)
			{
				int index = x + y * (NUMGRIDX+border) * blobbyGridSize + z * (NUMGRIDX+border) * blobbyGridSize * (NUMGRIDY+border) * blobbyGridSize;

				int X = x / blobbyGridSize - border/2;
				int Y = y / blobbyGridSize - border/2;
				int Z = z / blobbyGridSize - border/2;

				/*
				if(BOUNDED(X,Y,Z) && type[IX(X, Y, Z)] == FLUID && neighNoneSolid[IX(X, Y, Z)] - neighAir[IX(X, Y, Z)] == 6)
				{
					scalarField[index] = 2 * thresh;
					//scalarField[index] = -1;
					continue;
				}
				*/
				
				//origin blobby sited by Bridson
				double F = 0;
				//improved version by Zhu and Bridson
				//double Fai = -1;
				//Eigen::Vector3d Xu(0, 0, 0);
				//double ru = 0;
				std::vector<int> *list;
				for(int i = 0; i < 27; i++)
				{
					int xx = X + dir2[i][0];
					int yy = Y + dir2[i][1];
					int zz = Z + dir2[i][2];
					if(BOUNDED(xx, yy, zz))
						list = invertedList[IX(xx, yy, zz)];
					else
						continue;
					for(unsigned j = 0; j < list->size(); j++)
					{
						float x0 = (particles[list->at(j)].x + border/2) * blobbyGridSize;
						float y0 = (particles[list->at(j)].y + border/2) * blobbyGridSize;
						float z0 = (particles[list->at(j)].z + border/2) * blobbyGridSize;
						F += blobbyKernel(DISTANCE2(x0, y0, z0, x, y, z) * blobbyH2i);
						//double ker = blobbyKernel(DISTANCE2(x0, y0, z0, x, y, z) * blobbyH2i);
						//Xu += ker * Eigen::Vector3d(x0, y0, z0);
						//ru += ker * r; //for simple use the average r, but the radii to the closest particle is better
						//F += ker;
					}
				}
				scalarField[index] = F;

				//if(F > 0)
				//	Fai = DISTANCE(x, y, z, Xu[0]/F, Xu[1]/F, Xu[2]/F) - ru/F;
				//scalarField[index] = Fai;
			}
	
	clock_t timeCost = clock() - start;
	std::cout<<"Finish sampling in "<<timeCost<<" ms"<<std::endl;
	PRINT("Start to calculate the surface mesh using Marching Cubes...");
	start = clock();


	//thresh = 0;
	//create surface and write to obj
	CIsoSurface<float> builder;
	builder.GenerateSurface(scalarField, thresh, (NUMGRIDX+border)*blobbyGridSize-1, (NUMGRIDY+border)*blobbyGridSize-1, (NUMGRIDZ+border)*blobbyGridSize-1, 1, 1, 1);

	if(!builder.IsSurfaceValid())
		system("pause");

	char ext[] = ".obj";
	char file_path[100];
	sprintf(file_path, "%s/%03d%s", obj_dir, iteration, ext);
	std::ofstream fout(file_path);
	unsigned Nver = builder.m_nVertices;
	unsigned Ntri = builder.m_nTriangles;
	
	fout<<"#  v "<<Nver<<" f "<<Ntri<<std::endl;
	for(unsigned i = 0; i < Nver; i++)
	{
		fout << "v " <<builder.m_ppt3dVertices[i][0] - border/2 <<' ' \
					 <<builder.m_ppt3dVertices[i][1] - border/2 <<' ' \
					 <<builder.m_ppt3dVertices[i][2] - border/2 << std::endl;
	}

	for(unsigned i = 0; i < Ntri; i++)
	{
		fout<< "f " <<builder.m_piTriangleIndices[i*3]+1<<' '<<builder.m_piTriangleIndices[i*3+2]+1<<' '<<builder.m_piTriangleIndices[i*3+1]+1<<std::endl;
	}
	fout.close();

	timeCost = clock() - start;
	std::cout<<"Finish Marching Cubes in "<<timeCost<<" ms"<<std::endl<<std::endl;
	delete[] scalarField;
}

double FluidCube3D::blobbyKernel(double s2)
{
	if(s2 < 1)
		return (1-s2)*(1-s2)*(1-s2);
	else
		return 0;
}