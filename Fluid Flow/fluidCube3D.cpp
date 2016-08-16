#include "stdafx.h"
#ifndef SIMULATION_2D

#include "fluidCube3D.h"
#include <memory.h>
#include <math.h>
#include <ctime>
#include <cstdio>
#include <Eigen\Eigen>

extern float px;
extern float py;
extern float pz;

FluidCube3D::FluidCube3D(float viscosity, float fr, SCENETYPE sc, RENDERTYPE rt)
{
	size = (_X+2) * (_Y+2) * (_Z+2);
	h = _L / _X;
	h2 = h * h;
	hi = 1 / h;
	visc = viscosity;
	frameTime = 1.0 / fr;
	scene = sc;
	renderType = rt;
	ctime = 0;
	totalTime = 0;
	iteration = 0;

	max_vx = 0;
	max_vy = 0;
	max_vz = 0;
	max_p = 0;

	Vx = new float [size]; 
	Vy = new float [size]; 
	Vz = new float [size]; 
	Vx0 = new float [size]; 
	Vy0 = new float [size]; 
	Vz0 = new float [size]; 
	div = new float [size];
	type = new GRIDTYPE [size];

	//Advection using BFECC
	fai_b = new float [size];
	fai_f = new float [size];

	//Projection using Conjugate Gradient
	dir[0] = Eigen::Vector3i(0, 0, -1);
	dir[1] = Eigen::Vector3i(0, -1, 0);
	dir[2] = Eigen::Vector3i(-1, 0, 0);
	dir[3] = Eigen::Vector3i(1, 0, 0);
	dir[4] = Eigen::Vector3i(0, 1, 0);
	dir[5] = Eigen::Vector3i(0, 0, 1);
	pos2index = new int [size]; 
	neighNoneSolid = new int [size];
	neighAir = new int [size];
	neighbor = new int* [size];
	for(int i = 0; i < size; i ++)
		neighbor[i] = new int[6];

	//MAC
	type0 = new GRIDTYPE [size]; 
	invertedList = new std::vector<int>* [size];

	for(int i = 0; i < size; i++)
	{
		pos2index[i] = -1;
		type[i] = SOLID;
		type0[i] = SOLID;
		invertedList[i] = new std::vector<int>();
	}
	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
				type[IX(x,y,z)] = AIR;

	//init fluid
	//cube fall
	originFluid = 0;
	fluidNum = 0;
	switch(scene)
	{
	case CUBEFALL:
	{
		for(int z = _Z/3.0; z <= _Z/3.0*2; z++)
			for(int y = _Y/4.0; y <= _Y/2.0; y++)
				for(int x = _X/3.0; x <= _X/3.0*2; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
		break;
	}
	//sphere fall
	case SPHEREFALL:
	{
		int cx = _X/2;
		int cy = _Y/2;
		int cz = _Z/2;
		float R = _X/6;
		for(int z = 1; z <= _Z; z++)
			for(int y = 1; y <= _Y; y++)
				for(int x = 1; x <= _X; x++)
					if(DISTANCE(x, y, z, cx, cy, cz) <= R)
					{
						originFluid ++;
						fillParticleInGrid(x, y, z);
					}
		break;
	}
	//contain bottom
	case CONTAINER:
	{
		/*
		for(int z = _Z*3.0; z <= _Z*3.0/2; z++)
			for(int y = _Y/2.0; y <= _Y/3.0; y++)
				for(int x = _X/3.0; x <= _X/3.0*2; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
		*/

		for(int z = 1; z <= _Z; z++)
			for(int y = 1; y <= _Y/4.0; y++)
				for(int x = 1; x <= _X; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
		break;
	}
	//dam break
	case DAMBREAK:
	{
		for(int z = 1; z <= _Z/4.0; z++)
			for(int y = 1; y <= _Y/3.0*2; y++)
				for(int x = 1; x <= _X; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
		break;
	}
	case DOUBLEDAM:
	{
		for(int z = 1; z <= _Z/4.0; z++)
			for(int y = 1; y <= _Y/3.0*2; y++)
				for(int x = 1; x <= _X; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
		
		for(int z = _Z/4.0*3; z <= _Z; z++)
			for(int y = 1; y <= _Y/3.0*2; y++)
				for(int x = 1; x <= _X; x++)
				{
					originFluid ++;
					fillParticleInGrid(x, y, z);
				}
		break;
	}
	case EMPTY:
		break;
	}

#ifdef OBSTACLE
	for(int z = _Z/16.0*7; z <= _Z/16.0*9; z++)
		for(int y = 1; y <= _Y/3.0*2; y++)
			for(int x = _X/4.0; x <= _X/4.0*3; x++)
			{
				type[IX(x, y, z)] = type0[IX(x, y, z)] = SOLID;
			}
		/*for(int y = 1; y <= _Y; y++)
		for(int x = 1; x <= _X; x++)
			if(type[IX(x, y)] == SOLID)
			{
				for(int i = 0; i < 4; i++)
					if(type[IX(x+dir[i].x, y+dir[i].y)] == FLUID)
					{
						obstacle.push_back(Pos(x, y));
						break;
					}
			}
			else //FLUID
			{
				neighNum[IX(x, y)] = 0;
				for(int i = 0; i < 4; i++)
				{
					int xx = x+dir[i].x;
					int yy = y+dir[i].y;
					neighbor[IX(x,y)][i] = pos2index[IX(xx, yy)];
					if(type[IX(xx,yy)] != SOLID)
						neighNum[IX(x, y)] ++;
				}
			}
	*/
#endif

	memset(Vx, 0, sizeof(float) * size);
	memset(Vy, 0, sizeof(float) * size);
	memset(Vz, 0, sizeof(float) * size);
	memset(Vx0, 0, sizeof(float) * size);
	memset(Vy0, 0, sizeof(float) * size);
	memset(Vz0, 0, sizeof(float) * size);
	memset(div, 0, sizeof(float) * size);
	memset(fai_b, 0, sizeof(float) * size);
	memset(fai_f, 0, sizeof(float) * size);
}

FluidCube3D::~FluidCube3D()
{
	delete [] Vx;
	delete [] Vy;
	delete [] Vz;
	delete [] Vx0;
	delete [] Vy0;
	delete [] Vz0;
	delete [] div;
	delete [] type;

	delete [] pos2index;
	for(int i = 0; i < size; i++)
		delete[] neighbor[i];
	delete [] neighbor;
	delete [] neighNoneSolid;
	delete [] neighAir;

	delete [] fai_b;
	delete [] fai_f;

	delete [] type0;
	for(int i = 0; i < size; i++)
		delete[] invertedList[i];
	delete [] invertedList;

}

void FluidCube3D::vel_step()
{
	/*
	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	SWAP(Vz0, Vz);
	diffuseVelosity();
	set_bnd();
	*/

	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	SWAP(Vz0, Vz);
	advectVelosity();
	set_bnd();

	addForce();
	set_bnd();

	projectVelosity();
	set_bnd();

	//errorRemove();
}

void FluidCube3D::addForce()
{
	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
				if(type[IX(x, y, z)] == FLUID)
					Vy[IX(x, y, z)] -= dt * GRAVITY;
}

void FluidCube3D::diffuseVelosity()
{
	diffuse(1, Vx0, Vx, visc);
	diffuse(2, Vy0, Vy, visc);
	diffuse(3, Vz0, Vz, visc);
}

void FluidCube3D::advectVelosity()
{
	//BFECC
	/*
	advect(1, Vx0, fai_b, true);
	advect(1, fai_b, fai_f, false);
	for(int i = 0; i < size; i++)
		fai_b[i] = Vx0[i] + (Vx0[i] - fai_f[i]) * 0.5;
	advect(1, fai_b, Vx, true);

	advect(2, Vy0, fai_b, true);
	advect(2, fai_b, fai_f, false);
	for(int i = 0; i < size; i++)
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

	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
			{
				if(type[IX(x, y, z)] != FLUID)
					continue;

				div[IX(x, y, z)] = -h * (Vx[IX(x+1,y,z)]-Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y,z)]
										+ Vz[IX(x,y,z+1)]-Vz[IX(x,y,z)]);
				p[IX(x, y, z)] = 0;
			}
	
	for(int k = 0; k < ITERATION; k++)
	{
		for(int z = 1; z <= _Y; z++)
			for(int y = 1; y <= _Y; y++)
				for(int x = 1; x <= _X; x++)
					if(type[IX(x, y, z)] == FLUID)
						p[IX(x, y, z)] = (div[IX(x,y,z)] + p[IX(x-1,y,z)] + p[IX(x+1,y,z)] + p[IX(x,y-1,z)] + p[IX(x,y+1,z)]
										  + p[IX(x,y,z-1)] + p[IX(x,y,z+1)]) / 6;
	}

	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
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
	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
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
	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
			{
				if(type[IX(x, y, z)] != FLUID)
					continue;

				b[index++] =  -h * (Vx[IX(x+1,y,z)]-Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y,z)]
											+ Vz[IX(x,y,z+1)]-Vz[IX(x,y,z)]);
			}
	
	//std::cout<<b<<std::endl;
	p.resize(fluidNum);
	p = solver.solve(b);

	//for(int i = 100; i < 110; i++)
	//	cout<<p[i]<<' ';
	//cout<<endl;
	/*
	std::cout<<"bbbbbbb"<<std::endl;
	for(int y = 1; y <= _Y; y++)
		for(int x = 1; x <= _X; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			std::cout<<b[pos2index[IX(x,y)]]<<std::endl;
			if(x == _X)
				std::cout<<std::endl;
		}
	
	std::cout<<"ppppppp"<<std::endl;
	for(int y = 1; y <= _Y; y++)
	{
		float lastp = -1;
		for(int x = 1; x <= _X; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			//std::cout<<p[pos2index[IX(x,10)]]<<std::endl;
			//if(x == _X)
				//std::cout<<std::endl;
			float curp = p[pos2index[IX(x,y)]];
			if(x == 1)
				lastp = curp;
			else
			{
				if(curp != lastp)
				{
					std::cout<<x<<' '<<y<<' '<<lastp<<' '<<curp<<std::endl;
					lastp = curp;
				}
			}
		}
	}
	*/

	max_p = -9999;
	for(int i = 0; i < fluidNum; i++)
		if(p[i] > max_p)
			max_p = p[i];

	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
			{
				if(type[IX(x, y, z)] != FLUID)
					continue;

				double p1, p2;
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
		for(int z = 1; z <= _Z; z++)
			for(int y = 1; y <= _Y; y++)
				for(int x = 1; x <= _X; x++)
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
	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
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
	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
				if(type[IX(x, y, z)] == FLUID)
				{
					Pos pos = traceParticle(b, x, y, z, backward);
					u[IX(x, y, z)] = getVelosity(b, pos.x, pos.y, pos.z, u0);
				}
}

void FluidCube3D::set_bnd()
{
	//try to use free-slip condition
	/*
	for(int y = 1; y <= _Y; y++)
	{
		Vx[IX(1, y)] = -Vx[IX(2, y)];
		Vx[IX(_X+1, y)] = -Vx[IX(_X, y)];
	}
	for(int x = 1; x <= _X; x++)
	{
		Vy[IX(x, 1)] = -Vy[IX(x, 2)];
		Vy[IX(x, _Y+1)] = -Vy[IX(x, _Y)];
	}
	*/
	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
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

void FluidCube3D::simulate()
{	
	//while(true)
	//{
		bool draw = calculateTimeStep();
	
		updateParticles();
	
		updateGrid();
	
		set_bnd();
	
		vel_step();

		if(draw)
			render();

		report();
	//}
}

void FluidCube3D::output(float *u)
{
#ifdef OUTPUT 
	for(int z = 10; z <= 15; z++)
		for(int y = 10; y <= 15; y++)
			for(int x = 5; x <= 10; x++)
			{
				std::cout<<u[IX(x, y, z)]<<' ';
				if(x == 10)
					std::cout<<std::endl;
			}
	PRINT("=================================\n");
#endif
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

	//calculate divergence
	if(renderType == DIVERGENCE)
	{
		for(int z = 1; z <= _Z; z++)
			for(int y = 1; y <= _Y; y++)
				for(int x = 1; x <= _X; x++)
				{
					if(type[IX(x, y, z)] != FLUID)
						continue;
					div[IX(x, y, z)] = (Vx[IX(x+1,y,z)]-Vx[IX(x,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y,z)]
										+ Vz[IX(x,y,z+1)]-Vz[IX(x,y,z)]);
				}
	}

	glTranslatef(-LENGTH/2, -LENGTH/2, -LENGTH/2);

	for(int k = 0; k < _Z; k++)	
		for(int j = 0; j < _Y; j++)
			for(int i = 0; i < _X; i++)
			{
				int x = i+1;
				int y = j+1;
				int z = k+1;
				float color;
				if(type[IX(x, y, z)] == SOLID)
					glColor4f(0, 0.5, 0, 1);
				else if(type[IX(x, y, z)] == FLUID)
				{
					switch(renderType)
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
				else
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
	if(renderType == PARTICLE)
	{
		glColor3f(0, 0, 0.7);
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

bool FluidCube3D::calculateTimeStep()
{
	iteration ++;

	float max_v;
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

	if(max_v == 0)
		dt = frameTime;
	else
		dt = h / max_v;
	if(dt > frameTime)
		dt = frameTime;

	totalTime += dt;
	return true;

	if(ctime + dt >= frameTime)
	{
		dt = frameTime - ctime;
		ctime = 0;
		return true;
	}
	else
	{
		ctime += dt;
		return false;
	}
	//if(dt > h2 /(6*visc))
}

void FluidCube3D::updateParticles()
{
	for(unsigned i = 0; i < particles.size(); i++)
	{
		float x0 = particles[i].x;
		float y0 = particles[i].y;
		float z0 = particles[i].z;
		Velo v0 = getVelosity(x0, y0, z0, Vx, Vy, Vz);

		float x1 = x0 + dt * v0.x * hi;
		float y1 = y0 + dt * v0.y * hi;
		float z1 = z0 + dt * v0.z * hi;
		//if particle out of boundary??
		/*
		if(x1 < 1 || x1 >= _X+1 || y1 < 1 || y1 >= _Y+1 || z1 < 1 || z1 >= _Z+1)
		{
			std::cout<<"Particle out of bound"<<std::endl;
			REPORT(x1);
			REPORT(y1);
			REPORT(z1);
			system("pause");
		}
		if(x1 < 1)
			x1 = 1;
		else if(x1 >= _X+1)
			x1 = _X+0.5;//or 0.999 which is better?
		if(y1 < 1)
			y1 = 1;
		else if(y1 >= _Y+1)
			y1 = _Y+0.5;
		if(z1 < 1)
			z1 = 1;
		else if(z1 >= _Z+1)
			z1 = _Z+0.5;
		*/

		if(type[IX(int(x1), int(y1), int(z1))] == AIR)
		{
			particles[i] = Pos(x1, y1, z1);
			velosities[i] = v0;
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
			particles[i] = Pos(x1, y1, z1);
			velosities[i] = v0;
			continue;
		}

		Velo v1 = getVelosity(x1, y1, z1, Vx, Vy, Vz);
		x1 = x0 + dt * 0.5 * (v0.x + v1.x) * hi;
		y1 = y0 + dt * 0.5 * (v0.y + v1.y) * hi;
		z1 = z0 + dt * 0.5 * (v0.z + v1.z) * hi;
		//if particle out of boundary???
		/*
		if(x1 < 1 || x1 >= _X+1 || y1 < 1 || y1 >= _Y+1 || z1 < 1 || z1 >= _Z+1)
		{
			std::cout<<"Particle out of bound"<<std::endl;
			REPORT(x1);
			REPORT(y1);
			REPORT(z1);
			system("pause");
		}
		if(x1 < 1)
			x1 = 1;
		else if(x1 >= _X+1)
			x1 = _X+0.5;
		if(y1 < 1)
			y1 = 1;
		else if(y1 >= _Y+1)
			y1 = _Y+0.5;
		if(z1 < 1)
			z1 = 1;
		else if(z1 >= _Z+1)
			z1 = _Z+0.5;
		*/

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
		
		particles[i] = Pos(x1, y1, z1);
		velosities[i] = Velo((v0.x+v1.x)*0.5, (v0.y+v1.y)*0.5, (v0.z+v1.z)*0.5);
	}
}

void FluidCube3D::updateGrid()
{
	//swap 
	GRIDTYPE *tmp = type;
	type = type0;
	type0 = tmp;

	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
			{
				//if(type0[IX(x, y, z)] != SOLID)
				if(type0[IX(x, y, z)] == AIR || type0[IX(x, y, z)] == FLUID)
				{
					type[IX(x, y, z)] = AIR;
					invertedList[IX(x, y, z)]->clear();
				}
			}

	for(unsigned i = 0; i < particles.size(); i++)
	{
		int x = int(particles[i].x);
		int y = int(particles[i].y);
		int z = int(particles[i].z);
		type[IX(x, y, z)] = FLUID;
		invertedList[IX(x, y, z)]->push_back(i);
	}
	
	fluidNum = 0;
	for(int i = 0; i < size; i++)
		pos2index[i] = -1;
	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
				if(type[IX(x, y, z)] == FLUID)
					pos2index[IX(x, y, z)] = fluidNum++;

	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
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

				//Vy
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


	//init Matrix
	A = Eigen::SparseMatrix<double>(fluidNum, fluidNum);         // default is column major
	A.reserve(Eigen::VectorXi::Constant(fluidNum, 7));
	int index = 0;
	for(int z = 1; z <= _Z; z++)
		for(int y = 1; y <= _Y; y++)
			for(int x = 1; x <= _X; x++)
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
	A.makeCompressed();
	solver.compute(A);

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

	if(x < 0 || x >= _X+1 || y < 0 || y >= _Y+1 || z < 0 || z >= _Z+1)
	{
		std::cout<<"Get velosity out of bound"<<std::endl;
		REPORT(x);
		REPORT(y);
		REPORT(z);
		//system("pause");
	}
	if(x < 0)
		x = 0;
	if(x >= _X+1)
		x = _X+0.999;
	if(y < 0)
		y = 0;
	if(y >= _Y+1)
		y = _Y+0.999;
	if(z < 0)
		z = 0;
	if(z >= _Z+1)
		z = _Z+0.999;


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

Velo FluidCube3D::getVelosity(float x, float y, float z, float *vx, float *vy, float *vz)
{
	return Velo(getVelosity(1,x,y,z,vx), getVelosity(2,x,y,z,vy), getVelosity(3,x,y,z,vz));
}

Pos FluidCube3D::traceParticle(int index, int x, int y, int z, bool backward)
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

	Velo v0 = getVelosity(x0, y0, z0, Vx0, Vy0, Vz0);
	float t = (backward)? -dt : dt;
	//
	//Pos p = Pos(x0 + v0.x*t*hi, y0 + v0.y*t*hi);
	Velo v1 = getVelosity(x0 + v0.x*t*hi, y0 + v0.y*t*hi, z0 + v0.z*t*hi, Vx0, Vy0, Vz0);
	return Pos(x0 + 0.5*t*(v0.x+v1.x)*hi, y0 + 0.5*t*(v0.y+v1.y)*hi, z0 + 0.5*t*(v0.z+v1.z)*hi);
}

void FluidCube3D::errorRemove()
{
	double eps = 1e-12;

	for(int i = 0; i < size; i++)
	{
		if(fabs(Vx[i]) < eps)
			Vx[i] = 0;
		if(fabs(Vy[i]) < eps)
			Vy[i] = 0;
		if(fabs(Vz[i]) < eps)
			Vz[i] = 0;
	}
}

void FluidCube3D::fillParticleInGrid(int x, int y, int z)
{
	srand(time(0));
	int nump = NUMPERGRID;
	int sample = 50;
	float step = 1.0 / nump;
	for(int i = 0; i < nump; i++)
		for(int j = 0; j < nump; j++)
			for(int k = 0; k < nump; k++)
			{
				float x0 = 1.0f * (rand()%sample) / sample;
				float y0 = 1.0f * (rand()%sample) / sample;
				float z0 = 1.0f * (rand()%sample) / sample;
				//float x0 = i * step;
				//float y0 = j * step;
				//float z0 = k * step;
				particles.push_back(Pos(x+x0, y+y0, z+z0) );
				velosities.push_back(Velo(0, 0, 0));
			}
}

void FluidCube3D::report()
{
	REPORT(iteration);
	REPORT(dt);

	REPORT(max_vx);
	REPORT(max_vy);
	REPORT(max_vz);
	REPORT(max_p);
	float fluidShrink = 1.0 * fluidNum / originFluid;
	REPORT(fluidShrink);

	REPORT(totalTime);
	PRINT("======================================");
	PRINT("");
}
#endif