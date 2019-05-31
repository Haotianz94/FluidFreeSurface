#include "stdafx.h"
#include "fluidCube2D.h"

#include <GL/freeglut.h>
#include <memory.h>
#include <math.h>
#include <ctime>
#include <cstdio>


#define IX(x, y) ( (x) + (y) * (NUMGRIDW+2) )
#define IX2(x, y) ( (x) + (y) * (NUMGRIDW+2) * GRIDSIZE )
#define BOUNDED(x, y) ( (type[IX(int(x),int(y))] == SOLID || type[IX(int(x)+1,int(y)+1)] == SOLID)? false : true)
#define DISTANCE(x1, y1, x2, y2) ( sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)) )
#define DISTANCE2(x1, y1, x2, y2) ( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) )


FluidCube2D::FluidCube2D()
{
	float CUBELENGTH;
	assert(Configer::getConfiger()->getInt("Simulation2D", "NumGridW", NUMGRIDW)); 
	assert(Configer::getConfiger()->getInt("Simulation2D", "NumGridH", NUMGRIDH)); 
	assert(Configer::getConfiger()->getFloat("Simulation2D", "CubeLength", CUBELENGTH));

	NUMGRID = (NUMGRIDW + 2) * (NUMGRIDH + 2);
	h = CUBELENGTH / NUMGRIDH;
	h2 = h * h;
	hi = 1 / h;
	frameTime = 1.0 / FRAMERATE;
	ctime = 0;
	totalTime = 0;
	iteration = 0;

	max_vx = 0;
	max_vy = 0;
	max_p = 0;

	Vx = new float [NUMGRID]; 
	Vy = new float [NUMGRID]; 
	Vx0 = new float [NUMGRID]; 
	Vy0 = new float [NUMGRID]; 
	div = new float [NUMGRID];
	memset(Vx, 0, sizeof(float) * NUMGRID);
	memset(Vy, 0, sizeof(float) * NUMGRID);
	memset(Vx0, 0, sizeof(float) * NUMGRID);
	memset(Vy0, 0, sizeof(float) * NUMGRID);
	memset(div, 0, sizeof(float) * NUMGRID);
	type = new GridType [NUMGRID];
	type0 = new GridType [NUMGRID]; 
	invertedList = new std::vector<int>* [NUMGRID];

	dir[0] = Eigen::Vector2i(0, -1);
	dir[1] = Eigen::Vector2i(-1, 0);
	dir[2] = Eigen::Vector2i(1, 0);
	dir[3] = Eigen::Vector2i(0, 1);
	pos2index = new int [NUMGRID]; 
	neighNoneSolid = new int [NUMGRID];
	neighAir = new int [NUMGRID];
	neighbor = new int* [NUMGRID];
	for(int i = 0; i < NUMGRID; i++)
	{
		neighbor[i] = new int[4];
		pos2index[i] = -1;
		type[i] = SOLID;
		type0[i] = SOLID;
		invertedList[i] = new std::vector<int>();
	}
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
			type[IX(x,y)] = AIR;

	//Advection using BFECC
	fai_b = new float [NUMGRID];
	fai_f = new float [NUMGRID];
	memset(fai_b, 0, sizeof(float) * NUMGRID);
	memset(fai_f, 0, sizeof(float) * NUMGRID);

	//Blobby
	pixels = new float [NUMGRID * GRIDSIZE * GRIDSIZE* 3];
	pixelType = new GridType [NUMGRID * GRIDSIZE * GRIDSIZE];
	for(int y = 0; y < (NUMGRIDH+2)*GRIDSIZE; y++)
		for(int x = 0; x < (NUMGRIDW+2)*GRIDSIZE; x++)
		{
			int i = x / GRIDSIZE;
			int j = y / GRIDSIZE;
			int index = IX2(x, y);
			if(type[IX(i, j)] == SOLID)
			{
				pixels[index*3]		= 0; //R
				pixels[index*3 + 1]	= 0.7; //G
				pixels[index*3 + 2]	= 0; //B
				pixelType[index] = SOLID;
			}
			else
			{
				pixels[index*3]		= 0.5; //R
				pixels[index*3 + 1]	= 0.5; //G
				pixels[index*3 + 2]	= 0.5; //B
				pixelType[index] = AIR;
			}
		}
	dir2[0] = Eigen::Vector2i(0, 0);
	dir2[1] = Eigen::Vector2i(0, 1);
	dir2[2] = Eigen::Vector2i(0, -1);
	dir2[3] = Eigen::Vector2i(1, 0);
	dir2[4] = Eigen::Vector2i(-1, 0);
	dir2[5] = Eigen::Vector2i(1, 1);
	dir2[6] = Eigen::Vector2i(1, -1);
	dir2[7] = Eigen::Vector2i(-1, 1);
	dir2[8] = Eigen::Vector2i(-1, -1);


	//init fluid
	srand(time(0));
	originFluid = 0;
	fluidNum = 0;
	if(SCENETYPE.compare("CUBEFALL") == 0)
	{
		for(int y = NUMGRIDH/4.0; y <= NUMGRIDH/2.0; y++)
			for(int x = NUMGRIDW/3.0; x <= NUMGRIDW/3.0*2; x++)
			{
				originFluid ++;
				fillParticleInGrid(x, y);
			}
	}
	else if(SCENETYPE.compare("SPHEREFALL") == 0)
	{
		int cx = NUMGRIDW/2;
		int cy = NUMGRIDH/2;
		float R = NUMGRIDW/6;
		for(int y = 1; y <= NUMGRIDH; y++)
			for(int x = 1; x <= NUMGRIDW; x++)
				if(DISTANCE(x, y, cx, cy) <= R)
				{
					originFluid ++;
					fillParticleInGrid(x, y);
				}
	}
	else if(SCENETYPE.compare("CONTAINER") == 0)
	{
		for(int y = 1; y <= NUMGRIDH/4.0; y++)
			for(int x = 1; x <= NUMGRIDW; x++)
			{
				originFluid ++;
				fillParticleInGrid(x, y);
			}
	}
	else if(SCENETYPE.compare("DAMBREAK") == 0)
	{
		for(int y = 1; y <= NUMGRIDH/3.0*2; y++)
			for(int x = 1; x <= NUMGRIDW/4.0; x++)
			{
				originFluid ++;
				fillParticleInGrid(x, y);
			}
	}
	else if(SCENETYPE.compare("DOUBLEDAM") == 0)
	{
		for(int y = 1; y <= NUMGRIDH/3.0*2; y++)
			for(int x = 1; x <= NUMGRIDW/4.0; x++)
			{
				originFluid ++;
				fillParticleInGrid(x, y);
			}	

		for(int y = 1; y <= NUMGRIDH/3.0*2+8; y++)
			for(int x = NUMGRIDW/4.0*3; x <= NUMGRIDW; x++)
			{
				originFluid ++;
				fillParticleInGrid(x, y);
			}
	}
	else if(SCENETYPE.compare("NONE") == 0)
	{
		
	}
	else
	{
		PRINT("SceneType not known!");
		exit(0);
	}

	if(OBSTACLETYPE.compare("CENTERWALL"))
	{
		//int cx = NUMGRIDW / 2.0;
		//int cy = NUMGRIDH / 4.0;
		//int R = NUMGRIDH * 0.1;
		for(int y = 0; y <= NUMGRIDH/3.0; y++)
			for(int x = NUMGRIDW/8.0*3; x <= NUMGRIDW/8.0*5; x++)
			{
				//if(DISTANCE(x,y,cx,cy) <= R)
				type[IX(x, y)] = SOLID;
				for(int j = y*GRIDSIZE; j < (y+1)*GRIDSIZE; j++)
					for(int i = x*GRIDSIZE; i < (x+1)*GRIDSIZE; i++)
					{
						int index = IX2(i, j);
						pixels[index*3]		= 0; //R
						pixels[index*3 + 1]	= 0.7; //G
						pixels[index*3 + 2]	= 0; //B
						pixelType[index] = SOLID;
					}
			}
		/*for(int y = 1; y <= NUMGRIDH; y++)
			for(int x = 1; x <= NUMGRIDW; x++)
				if(type[IX(x, y)] == SOLID)
				{
					for(int i = 0; i < 4; i++)
						if(type[IX(x+dir[i].x, y+dir[i].y)] == FLUID)
						{
							obstacle.push_back(Pos2D(x, y));
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
	}
}

FluidCube2D::~FluidCube2D()
{
	delete [] Vx;
	delete [] Vy;
	delete [] Vx0;
	delete [] Vy0;
	delete [] pixels;
	delete [] pixelType;
}

void FluidCube2D::vel_step()
{
	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	advectVelosity();
	set_bnd();

	addForce();
	set_bnd();
	/*
	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	diffuseVelosity();
	set_bnd();
	*/
	projectVelosity();
	set_bnd();

	//errorRemove();
}

void FluidCube2D::addForce()
{
#pragma omp parallel for
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
			if(type[IX(x, y)] == FLUID)
				Vy[IX(x, y)] -= dt * GRAVITY;
}

void FluidCube2D::diffuseVelosity()
{
	diffuse(1, Vx0, Vx, VISCOSITY);
	diffuse(2, Vy0, Vy, VISCOSITY);
}

void FluidCube2D::advectVelosity()
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
}

void FluidCube2D::projectVelosity()
{

	max_vx = 0;
	max_vy = 0;
	float *div = fai_b;

#ifdef GAUSS_SEIDEL
	//Gauss_Seidel
	float *p = fai_f;
	//float *div = Vy0;

	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			div[IX(x, y)] = -h * (Vx[IX(x+1,y)]-Vx[IX(x,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y)]);
			p[IX(x, y)] = 0;
		}
	//set_bnd(0, div);
	//set_bnd(0, p);
	
	for(int k = 0; k < ITERATION; k++)
	{
		for(int y = 1; y <= NUMGRIDH; y++)
			for(int x = 1; x <= NUMGRIDW; x++)
				if(type[IX(x, y)] == FLUID)
					p[IX(x, y)] = (div[IX(x,y)] + p[IX(x,y)] + p[IX(x+1,y)] + p[IX(x,y)] + p[IX(x,y+1)]) / 4;
	}

#pragma omp parallel for
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			double p1, p2;
			p2 = p[IX(x,y)];
			//Vx
			if(type[IX(x-1, y)] == AIR)
				p1 = 0;
			else if(type[IX(x-1, y)] == FLUID)
				p1 = p[IX(x-1,y)];
			else
				p1 = p[IX(x,y)];
			Vx[IX(x, y)] -= (p2 - p1) * hi;
			
			//Vy
			if(type[IX(x, y-1)] == AIR)
				p1 = 0;
			else if(type[IX(x, y-1)] == SOLID)
				p1 = p[IX(x,y)];
			else
			{
				if(type[IX(x, y-1)] != FLUID)
					*(int*)0 = 0;
				p1 = p[IX(x,y-1)];
			}
			Vy[IX(x, y)] -= (p2 - p1) * hi;

			if(fabsf(Vx[IX(x, y)]) > max_vx)
				max_vx = fabsf(Vx[IX(x, y)]);
			if(fabsf(Vy[IX(x, y)]) > max_vy)
				max_vy = fabsf(Vy[IX(x, y)]);
		}
	//set_bnd(1, Vx);
	//set_bnd(2, Vy);
#else
	//Conjugate Gradient
	/*
	//check div before project
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;
			div[IX(x, y)] = (Vx[IX(x+1,y)]-Vx[IX(x,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y)]);
		}
	output(div);
	*/
	Eigen::VectorXd b(fluidNum);
	int index = 0;

// #pragma omp parallel for (cannot be parallel)
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			b[index++] =  -h * (Vx[IX(x+1,y)]-Vx[IX(x,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y)]);
		}
	
	//std::cout<<b<<std::endl;

	p.resize(fluidNum);
	p = solver.solve(b);

	float max_p = -9999;
	for(int i = 0; i < fluidNum; i++)
		if(p[i] > max_p)
			max_p = p[i];

#pragma omp parallel for
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			float p1, p2;
			p2 = p[pos2index[IX(x,y)]];
			//Vx
			if(type[IX(x-1, y)] == AIR)
				p1 = 0;
			else if(type[IX(x-1, y)] == FLUID)
				p1 = p[pos2index[IX(x-1,y)]];
			else//solid and flowin
				p1 = p[pos2index[IX(x,y)]];
			Vx[IX(x, y)] -= (p2 - p1) * hi;
			
			//Vy
			if(type[IX(x, y-1)] == AIR)
				p1 = 0;
			else if(type[IX(x, y-1)] == SOLID)
				p1 = p[pos2index[IX(x,y)]];
			else
				p1 = p[pos2index[IX(x,y-1)]];
			Vy[IX(x, y)] -= (p2 - p1) * hi;
		

			if(fabsf(Vx[IX(x, y)]) > max_vx)
				max_vx = fabsf(Vx[IX(x, y)]);
			if(fabsf(Vy[IX(x, y)]) > max_vy)
				max_vy = fabsf(Vy[IX(x, y)]);
		}
	
#endif

}

void FluidCube2D::diffuse(int b, float *u0, float *u, float diffusion)
{
	float a = dt * diffusion / h2;

	//Gauss Seidel relexation
	//in this way, the initcial value for u may be important 
	/*
	for(int k = 0; k < ITERATION; k++)
	{
		for(int y = 1; y <= NUMGRIDH; y++)
			for(int x = 1; x <= NUMGRIDW; x++)
				if(type[IX(x, y)] == FLUID)
				{	
					int fnum = 0;
					float v = 0;
					for(int i = 0; i < 4; i++)
					{
						int xx = x + dir[i].x;
						int yy = y + dir[i].y;
						if(type[IX(xx, yy)] == FLUID)
						{
							v += u[IX(xx, yy)];
							fnum ++;
						}
					}
					u[IX(x, y)] = (u0[IX(x, y)] + a * v) / (1+fnum*a);
				}
	}
	*/

	//can also try unstable way
#pragma omp parallel for
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
			if(type[IX(x, y)] == FLUID)
			{
				u[IX(x, y)] = u0[IX(x, y)];
				for(int i = 0; i < 4; i++)
				{
					int xx = x + dir[i][0];
					int yy = y + dir[i][1];
					if(type[IX(xx, yy)] == FLUID)
					{
						u[IX(x, y)] += a * (u0[IX(xx, yy)] - u0[IX(x, y)]);
					}	
				}
			}
					
}

void FluidCube2D::advect(int b, float *u0, float *u,  bool backward)
{
#pragma omp parallel for
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;
		
			Pos2D pos = traceParticle(b, x, y, backward);
		
			u[IX(x, y)] = getVelosity(b, pos.x, pos.y, u0);
		}
}

void FluidCube2D::set_bnd()
{
	//try to use free-slip condition
	/*
	for(int y = 1; y <= NUMGRIDH; y++)
	{
		Vx[IX(1, y)] = -Vx[IX(2, y)];
		Vx[IX(NUMGRIDW+1, y)] = -Vx[IX(NUMGRIDW, y)];
	}
	for(int x = 1; x <= NUMGRIDW; x++)
	{
		Vy[IX(x, 1)] = -Vy[IX(x, 2)];
		Vy[IX(x, NUMGRIDH+1)] = -Vy[IX(x, NUMGRIDH)];
	}
	*/
#pragma omp parallel for
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			if(type[IX(x-1, y)] == SOLID)
			{
				Vx[IX(x, y)] = 0;
				Vy[IX(x-1, y)] = Vy[IX(x, y)];
			}
			if(type[IX(x+1, y)] == SOLID)
			{
				Vx[IX(x+1, y)] = 0;
				Vy[IX(x+1, y)] = Vy[IX(x, y)];
			}
			if(type[IX(x, y-1)] == SOLID)
			{
				Vy[IX(x, y)] = 0;
				Vx[IX(x, y-1)] = Vx[IX(x, y)];
			}
			if(type[IX(x, y+1)] == SOLID)
			{
				Vy[IX(x, y+1)] = 0;
				Vx[IX(x, y+1)] = Vx[IX(x, y)];
			}

			switch(neighAir[IX(x, y)])
			{
			case 0: 
				break;
			case 1: 
				if(type[IX(x-1, y)] == AIR)
					Vx[IX(x,y)] = Vx[IX(x+1,y)] + Vy[IX(x,y+1)] - Vy[IX(x,y)];
				else if(type[IX(x+1, y)] == AIR)
					Vx[IX(x+1,y)] = Vx[IX(x,y)] - Vy[IX(x,y+1)] + Vy[IX(x,y)];
				else if(type[IX(x, y-1)] == AIR)
					Vy[IX(x,y)] = Vy[IX(x,y+1)] + Vx[IX(x+1,y)] - Vx[IX(x,y)];
				else if(type[IX(x, y+1)] == AIR)
					Vy[IX(x,y+1)] = Vy[IX(x,y)] - Vx[IX(x+1,y)] + Vx[IX(x,y)];
				break;
			case 2:
				if(type[IX(x-1, y)] == AIR)
					if(type[IX(x+1, y)] == AIR)
						break;
					else
						Vx[IX(x,y)] = Vx[IX(x+1,y)];
				else
					if(type[IX(x+1, y)] == AIR)
						Vx[IX(x+1,y)] = Vx[IX(x,y)];
					else
						break;
				if(type[IX(x, y-1)] == AIR)
					Vy[IX(x,y)] = Vy[IX(x,y+1)];
				else
					Vy[IX(x,y+1)] = Vy[IX(x,y)];
				break;
			case 3:
				
				if(type[IX(x+1, y)] != AIR)
					Vx[IX(x,y)] = Vx[IX(x+1,y)] + Vy[IX(x,y+1)] - Vy[IX(x,y)];
				else if(type[IX(x-1, y)] != AIR)
					Vx[IX(x+1,y)] = Vx[IX(x,y)] - Vy[IX(x,y+1)] + Vy[IX(x,y)];
				else if(type[IX(x, y+1)] != AIR)
					Vy[IX(x,y)] = Vy[IX(x,y+1)] + Vx[IX(x+1,y)] - Vx[IX(x,y)];
				else if(type[IX(x, y-1)] != AIR)
					Vy[IX(x,y+1)] = Vy[IX(x,y)] - Vx[IX(x+1,y)] + Vx[IX(x,y)];
				break;
				
				/*
				if(type[IX(x+1, y)] != AIR)
					Vx[IX(x,y)] = Vx[IX(x+1,y)];
				else if(type[IX(x-1, y)] != AIR)
					Vx[IX(x+1,y)] = Vx[IX(x,y)];
				else if(type[IX(x, y+1)] != AIR)
					Vy[IX(x,y)] = Vy[IX(x,y+1)];
				else if(type[IX(x, y-1)] != AIR)
					Vy[IX(x,y+1)] = Vy[IX(x,y)];
				*/
				break;
			case 4:
				break;
			}
		}
}

void FluidCube2D::simulate()
{	
	REPORT("enter simulate");
	clock_t start = clock();

	bool draw = calculateTimeStep();

	addFlowIn();

	updateParticles();
	REPORT("updateParticles");

	updateGrid();
	REPORT("updateGrid");

	set_bnd();
	REPORT("set_bnd");

	vel_step();
	REPORT("vel_step");

	clock_t simTime = clock() - start;
	report(simTime);

	if(draw)
		render();
}

void FluidCube2D::output(float *u)
{
	if(!DEBUGPRINT)
		return; 
	for(int y = 10; y <= 15; y++)
		for(int x = 5; x <= 10; x++)
		{
			std::cout<<u[IX(x, y)]<<' ';
			if(x == 10)
				std::cout<<std::endl;
		}
	LOGSEG;
}

void FluidCube2D::render()
{
	//identify that we are currently modifying the projection matrix
	glMatrixMode(GL_PROJECTION);
	//load the identity matrix (clear old matrix)
	glLoadIdentity();

	//get the current window width/height
	int w = glutGet(GLUT_WINDOW_WIDTH);
	int h = glutGet(GLUT_WINDOW_HEIGHT);

	//initialize the screen coordinates
	//glOrtho(-w/2, w/2-1, -h/2, h/2-1, -1, 1);
	glOrtho(-10, w+10, -10, h+10, -1, 1);

	//now we are editing the modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//clear the screen to a desired color in range [0-1] RGBA
	glClearColor(0.5, 0.5, 0.5, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//enable blending for translucency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

	//calculate divergence
	if(RENDERTYPE == DIVERGENCE)
	{
#pragma omp parallel for
		for(int y = 1; y <= NUMGRIDH; y++)
			for(int x = 1; x <= NUMGRIDW; x++)
			{
				if(type[IX(x, y)] != FLUID)
					continue;
				div[IX(x, y)] = (Vx[IX(x+1,y)]-Vx[IX(x,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y)]);
			}
	}

	if(RENDERTYPE == PARTICLE)
	{
		for(int i = 0; i <= NUMGRIDW+1; i++)
		for(int j = 0; j <= NUMGRIDH+1; j++)
		{
			int x = i;
			int y = j;
			float color;

			if(type[IX(x, y)] == SOLID)
				glColor3f(0, 0.5, 0);
			else if(type[IX(x, y)] == FLUIDIN)
				glColor3f(0, 0, 0.7);
			else
				continue;

			glBegin(GL_QUADS);
			glVertex2f(i*GRIDSIZE, j*GRIDSIZE);
			glVertex2f((i+1)*GRIDSIZE, j*GRIDSIZE);
			glVertex2f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE);
			glVertex2f(i*GRIDSIZE, (j+1)*GRIDSIZE);
			glEnd();
		}

		//draw particles
		glColor3f(0, 0, 0.7);
		glBegin(GL_POINTS);
		for(unsigned i = 0; i < particles.size(); i++)
		{
			glVertex2f(particles[i].x*GRIDSIZE, particles[i].y*GRIDSIZE);
		}
		glEnd();
	}
	else if(RENDERTYPE == BLOBBY)
	{
		double r = 1.0 * GRIDSIZE / PARTICLEPERGRID;
		double h = 3 * r;
		double h2i = 1 / (h*h);
		double thresh = blobbyKernel( (r*r) / (h*h) );
		
		for(int y = 0; y < (NUMGRIDH+2)*GRIDSIZE; y++)
		for(int x = 0; x < (NUMGRIDW+2)*GRIDSIZE; x++)
		{
			int index = IX2(x, y);

			int X = x / GRIDSIZE;
			int Y = y / GRIDSIZE;
			/*
			if(neighAir[IX(i, j)] == 0 || neighAir[IX(i, j)] == 4)
			{
				pixels[index * 3] = 0;
				pixels[index * 3 + 1] = 0;
				pixels[index * 3 + 2] = 0.7;
				continue;
			}
			else
			{
				pixels[index * 3] = 0.5;
				pixels[index * 3 + 1] = 0.5;
				pixels[index * 3 + 2] = 0.5;
			}
			continue;
			*/

			if(pixelType[index] != SOLID && pixelType[index] != FLUIDIN)
			{
				//origin blobby sited by Bridson
				/*
				double F = 0;
				bool inside = false;
				std::vector<int> *list;
				for(int i = 0; i < 9; i++)
				{
					list = invertedList[IX(X+dir2[i][0], Y+dir2[i][1])];
					for(unsigned j = 0; j < list->size(); j++)
						F += blobbyKernel(DISTANCE2(particles[list->at(j)].x*GRIDSIZE, particles[list->at(j)].y*GRIDSIZE, x, y) * h2i);
					if(F >= thresh)
					{
						pixels[index * 3] = 0;
						pixels[index * 3 + 1] = 0;
						pixels[index * 3 + 2] = 0.7;
						inside = true;
						break;
					}
				}
				if(!inside)
				{
					pixels[index * 3] = 0.5;
					pixels[index * 3 + 1] = 0.5;
					pixels[index * 3 + 2] = 0.5;
				}
				*/
				//improved version by Zhu and Bridson

				double Fai = -1;
				Eigen::Vector2d Xu(0, 0);
				double ru = 0, F = 0;
				std::vector<int> *list;
				for(int i = 0; i < 9; i++)
				{
					list = invertedList[IX(X+dir2[i][0], Y+dir2[i][1])];
					for(unsigned j = 0; j < list->size(); j++)
					{
						float x0 = particles[list->at(j)].x*GRIDSIZE;
						float y0 = particles[list->at(j)].y*GRIDSIZE;
						double ker = blobbyKernel(DISTANCE2(x0, y0, x, y) * h2i);
						Xu += ker * Eigen::Vector2d(x0, y0);
						ru += ker * r; //for simple use the average r, but the radii to the closest particle is better
						F += ker;
					}
				}
				if(F > 0)
				{
					Fai = DISTANCE(x, y, Xu[0]/F, Xu[1]/F) - ru/F;
				}
				if(fabs(Fai) < 0.5)
				{
					pixels[index * 3] = 0;
					pixels[index * 3 + 1] = 0;
					pixels[index * 3 + 2] = 0.7;
				}
				/*
				else if(F == 0)
				{
					pixels[index * 3] = 0;
					pixels[index * 3 + 1] = 0.7;
					pixels[index * 3 + 2] = 0.7;
				}
				*/
				else
				{
					pixels[index * 3] = 0.5;
					pixels[index * 3 + 1] = 0.5;
					pixels[index * 3 + 2] = 0.5;
				}
			}
		}
		glDrawPixels((NUMGRIDW+2)*GRIDSIZE, (NUMGRIDH+2)*GRIDSIZE, GL_RGB, GL_FLOAT, pixels);
	}
	else
	{
		for(int i = 0; i <= NUMGRIDW+1; i++)
		for(int j = 0; j <= NUMGRIDH+1; j++)
		{
			int x = i;
			int y = j;
			float color;

			if(type[IX(x, y)] == SOLID)
				glColor3f(0, 0.5, 0);
			else if(type[IX(x, y)] == FLUIDIN)
				glColor3f(0, 0, 0.7);
			else if(type[IX(x, y)] == FLUID)
			{
				switch(RENDERTYPE)
				{
				case FLUIDGRID:
					glColor3f(0, 0, 0.7);
					break;
				case PRESSURE:
					glColor3f(p[pos2index[IX(x,y)]]/max_p, 0, 0);
					break;
				case VELOSITYY:
					if(Vy[IX(x, y)] >= 0)
						glColor3f(Vy[IX(x, y)]/max_vy, 0, 0);
					else
						glColor3f(0, -Vy[IX(x, y)]/max_vy, 0);
					break;
				case VELOSITYX:
					if(Vx[IX(x, y)] >= 0)
						glColor3f(Vx[IX(x, y)]/max_vx, 0, 0);
					else
						glColor3f(0, -Vx[IX(x, y)]/max_vx, 0);
					break;
				case DIVERGENCE:
					if(fabs(div[IX(x, y)]) > 1e-5)
						glColor3f(1, 0, 0);
					else
						glColor3f(0, 0, 0.7);
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

			glBegin(GL_QUADS);
			glVertex2f(i*GRIDSIZE, j*GRIDSIZE);
			glVertex2f((i+1)*GRIDSIZE, j*GRIDSIZE);
			glVertex2f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE);
			glVertex2f(i*GRIDSIZE, (j+1)*GRIDSIZE);
			glEnd();
			//if(GRIDSIZE >= 10 && type[IX(x, y)] == FLUID)
			//	draw_velo(i, j, Vx[IX(x, y)], Vy[IX(x, y)]);
		}
	}

	glutSwapBuffers();
}

void FluidCube2D::draw_velo(int i, int j, float vx, float vy)
{
	float x1 = i * GRIDSIZE;
	float y1 = (j+0.5) * GRIDSIZE;
	float x2 = (i+0.5) * GRIDSIZE;
	float y2 = j * GRIDSIZE;
	float vl = sqrtf(vx*vx + vy*vy);
	//float seita = atanf(fabs(vy / vx));

	//float max_v = max_vx > max_vy? max_vx : max_vy;
	//if(max_v == 0)
	//	return;

	/*if(fabs(vx) < eps && fabs(vy) < eps)
		return;
	int LOG = 10;
	max_v = logf(max_v) + LOG;
	float dx = 0.6 * GRIDSIZE * (logf(fabsf(vx))+LOG) / max_v * ((vx>0)? 1:-1);
	float dy = 0.6 * GRIDSIZE * (logf(fabsf(vy))+LOG) / max_v * ((vy>0)? 1:-1);
	*/

	float dx = 0.5 * GRIDSIZE * vx;
	float dy = 0.5 * GRIDSIZE * vy;

	//float dx = 20 * GRIDSIZE * vx ;
	//float dy = 20 * GRIDSIZE * vy ;

	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex2f(x1, y1);
	glVertex2f(x1 + dx, y1);
	
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex2f(x2, y2);
	glVertex2f(x2, y2 + dy);
	glEnd();

	glBegin(GL_POINTS);
	glColor3f(0, 0, 0);
	glVertex2f(x1, y1);
	glVertex2f(x2, y2);
	glEnd();
}

bool FluidCube2D::calculateTimeStep()
{
	iteration ++;

	float max_v;
	if(max_vx > max_vy)
		max_v = max_vx;
	else
		max_v = max_vy;
	if(max_v == 0)
		dt = frameTime;
	else
		dt = h / max_v;
	if(dt > frameTime)
		dt = frameTime;
	
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

void FluidCube2D::updateParticles()
{
	for(unsigned i = 0; i < particles.size(); i++)
	{
		float x0 = particles[i].x;
		float y0 = particles[i].y;
		float vx0 = getVelosity(1, x0, y0, Vx);
		float vy0 = getVelosity(2, x0, y0, Vy);

		float x1 = x0 + dt * vx0 * hi;
		float y1 = y0 + dt * vy0 * hi;
		//if particle out of boundary??
		/*
		if(x1 < 1 || x1 >= NUMGRIDW+1 || y1 < 1 || y1 >= NUMGRIDH+1)
		{
			std::cout<<"Particle out of bound"<<std::endl;
			REPORT(x1);
			REPORT(y1);
			//system("pause");
		}
		if(x1 < 1)
			x1 = 1;
		else if(x1 >= NUMGRIDW+1)
			x1 = NUMGRIDW+0.999;
		if(y1 < 1)
			y1 = 1;
		else if(y1 >= NUMGRIDH+1)
			y1 = NUMGRIDH+0.999;
		*/

		if(type[IX(int(x1), int(y1))] == AIR)
		{
			particles[i] = Pos2D(x1, y1);
			velosities[i] = Velo2D(vx0, vy0);
			continue;
		}
		else if(type[IX(int(x1), int(y1))] == SOLID)
		{
			if(int(x1) > int(x0))
				x1 = int(x0)+0.99;
			else if(int(x1) < int(x0))
				x1 = int(x0)+0.01;
			if(int(y1) > int(y0))
				y1 = int(y0)+0.99;
			else if(int(y1) < int(y0))
				y1 = int(y0)+0.01;

			particles[i] = Pos2D(x1, y1);
			velosities[i] = Velo2D(vx0, vy0);
			continue;
		}

		float vx1 = getVelosity(1, x1, y1, Vx);
		float vy1 = getVelosity(2, x1, y1, Vy);

		x1 = x0 + dt * 0.5 * (vx0 + vx1) * hi;
		y1 = y0 + dt * 0.5 * (vy0 + vy1) * hi;
		//if particle out of boundary???
		/*
		if(x1 < 1 || x1 >= NUMGRIDW+1 || y1 < 1 || y1 >= NUMGRIDH+1)
		{
			std::cout<<"Particle out of bound"<<std::endl;
			REPORT(x1);
			REPORT(y1);
			//system("pause");
		}
		if(x1 < 1)
			x1 = 1;
		else if(x1 >= NUMGRIDW+1)
			x1 = NUMGRIDW+0.999;
		if(y1 < 1)
			y1 = 1;
		else if(y1 >= NUMGRIDH+1)
			y1 = NUMGRIDH+0.999;
		*/
		if(type[IX(int(x1), int(y1))] == SOLID)
		{
			if(int(x1) > int(x0))
				x1 = int(x0)+0.99;
			else if(int(x1) < int(x0))
				x1 = int(x0)+0.01;
			if(int(y1) > int(y0))
				y1 = int(y0)+0.99;
			else if(int(y1) < int(y0))
				y1 = int(y0)+0.01;
		}

		particles[i] = Pos2D(x1, y1);
		velosities[i] = Velo2D((vx0+vx1)*0.5, (vy0+vy1)*0.5);
	}
}

void FluidCube2D::updateGrid()
{
	//swap 
	GridType *tmp = type;
	type = type0;
	type0 = tmp;

#pragma omp parallel for
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
		{
			if(type0[IX(x, y)] == AIR || type0[IX(x, y)] == FLUID)
			{
				type[IX(x, y)] = AIR;
				invertedList[IX(x, y)]->clear();
			}
		}
	for(unsigned i = 0; i < particles.size(); i++)
	{
		int x = particles[i].x;
		int y = particles[i].y;
		type[IX(x, y)] = FLUID;
		invertedList[IX(x, y)]->push_back(i);
	}
	
	fluidNum = 0;
	for(int i = 0; i < NUMGRID; i++)
		pos2index[i] = -1;

// #pragma omp parallel for (cannot be parallel)
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
			if(type[IX(x,y)] == FLUID)
				pos2index[IX(x, y)] = fluidNum++;

#pragma omp parallel for
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
		{
			if(type[IX(x,y)] == FLUID)
			{
				neighNoneSolid[IX(x, y)] = 0;
				neighAir[IX(x, y)] = 0;
				for(int i = 0; i < 4; i++)
				{
					int xx = x+dir[i][0];
					int yy = y+dir[i][1];
					neighbor[IX(x,y)][i] = pos2index[IX(xx, yy)];
					if(type[IX(xx,yy)] == AIR || type[IX(xx,yy)] == FLUID)
						neighNoneSolid[IX(x, y)] ++;
					if(type[IX(xx,yy)] == AIR)
						neighAir[IX(x, y)] ++;
				}
			}

			// set velosity for new fluid cell
			float downRate = 0.7;
			if(type0[IX(x, y)] == AIR && type[IX(x, y)] == FLUID)
			{
				float vsum = 0;
				int nsum = 0;
				float dist = 100;
				int pid = -1;
				//Vx
				std::vector<int> *list = invertedList[IX(x, y)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					float x0 = particles[list->at(i)].x;
					float y0 = particles[list->at(i)].y;
					if( x0 <= x + 0.5)
					{
						nsum ++;
						vsum += velosities[list->at(i)].x;
					}
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
				list = invertedList[IX(x-1, y)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					float x0 = particles[list->at(i)].x;
					float y0 = particles[list->at(i)].y;
					if(x0 >= x - 0.5)
					{
						nsum ++;
						vsum += velosities[list->at(i)].x;
					}
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
				//if(nsum > 0)
				//	Vx[IX(x, y)] = vsum / nsum * downRate;
				//else
				{
					//for simple only take the nearest particle in the cell
					Vx[IX(x, y)] = velosities[pid].x * downRate;
				}

				dist = 100;
				vsum = 0;
				nsum = 0;
				//Vy
				list = invertedList[IX(x, y)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					float x0 = particles[list->at(i)].x;
					float y0 = particles[list->at(i)].y;
					if( y0 <= y + 0.5)
					{
						nsum ++;
						vsum += velosities[list->at(i)].y;
					}
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
				list = invertedList[IX(x, y-1)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					float x0 = particles[list->at(i)].x;
					float y0 = particles[list->at(i)].y;
					if(y0 >= y - 0.5)
					{
						nsum ++;
						vsum += velosities[list->at(i)].y;
					}
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
				//if(nsum > 0)
				//	Vy[IX(x, y)] = vsum / nsum;
				//else
				{
					//for simple only take the nearest particle in the cell
					Vy[IX(x, y)] = velosities[pid].y * downRate;
				}
			}
		}


	//init Matrix
// #pragma omp parallel for (cannot be parallel)
	A = Eigen::SparseMatrix<double>(fluidNum, fluidNum);         // default is column major
	A.reserve(Eigen::VectorXi::Constant(fluidNum, 5));
	int index = 0;
	for(int y = 1; y <= NUMGRIDH; y++)
		for(int x = 1; x <= NUMGRIDW; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			//SparseVecf A0(fluidNum);
			//A0.Begin();
			for(int i = 0; i < 2; i++)
			{
				int neighid = neighbor[IX(x,y)][i];
				if(neighid != -1)
				{
					//A0.AddNZElt(neighid, -1);
					A.insert(index, neighid) = -1;
				}
			}
			//A0.AddNZElt(index, neighNum[IX(x, y)]);
			A.insert(index, index) = neighNoneSolid[IX(x, y)];
			for(int i = 2; i < 4; i++)
			{
				int neighid = neighbor[IX(x,y)][i];
				if(neighid != -1)
				{
					//A0.AddNZElt(neighid, -1);
					A.insert(index, neighid) = -1;
				}
			}
			//A0.End();
			//cout<<A0<<endl;
			//A0.SetElts(index, neighs, VL_SV_END);
			//(*A)[index++] = A0;
			index ++;
		}
	//std::cout<<A<<std::endl;
	//check symmetry
	//for(int i = 0; i < fluidNum; i++)
	//	for(int j = i+1; j < fluidNum; j++)
	//		A.
	A.makeCompressed();
	solver.compute(A);

}

float FluidCube2D::getVelosity(int index, float x, float y, float *u)
{
	if(index == 1)
	{
		y -= 0.5;
	}
	else
	{
		x -= 0.5;
	}
	/*
	if(x < 0 || x >= NUMGRIDW+1 || y < 0 || y >= NUMGRIDH+1)
	{
		std::cout<<"Get velosity out of bound"<<std::endl;
		REPORT(x);
		REPORT(y);
		//system("pause");
	}
	*/
	if(x < 0)
		x = 0;
	if(x >= NUMGRIDW+1)
		x = NUMGRIDW+0.999;
	if(y < 0)
		y = 0;
	if(y >= NUMGRIDH+1)
		y = NUMGRIDH+0.999;


	int i0 = int(x), i1 = i0 + 1;
	int j0 = int(y), j1 = j0 + 1;
	float s1 = x - i0, s0 = 1 - s1;
	float t1 = y - j0, t0 = 1 - t1;

	return s0 * (t0*u[IX(i0,j0)] + t1*u[IX(i0,j1)]) +
		   s1 * (t0*u[IX(i1,j0)] + t1*u[IX(i1,j1)]);
}

Velo2D FluidCube2D::getVelosity(float x, float y, float *vx, float *vy)
{
	return Velo2D(getVelosity(1,x,y,vx), getVelosity(2,x,y,vy));
}

Pos2D FluidCube2D::traceParticle(int index, int x, int y, bool backward)
{
	float x0, y0;
	if(index == 1)
	{
		x0 = x;
		y0 = y + 0.5;
	}
	else
	{
		x0 = x + 0.5;
		y0 = y;
	}
	Velo2D v0 = getVelosity(x0, y0, Vx0, Vy0);
	float t = (backward)? -dt : dt;
	//
	//Pos2D p = Pos2D(x0 + v0.x*t*hi, y0 + v0.y*t*hi);
	Velo2D v1 = getVelosity(x0 + v0.x*t*hi, y0 + v0.y*t*hi, Vx0, Vy0);
	return Pos2D(x0 + 0.5*t*(v0.x+v1.x)*hi, y0 + 0.5*t*(v0.y+v1.y)*hi);
}

void FluidCube2D::errorRemove()
{
	double eps = 1e-12;

	for(int i = 0; i < NUMGRID; i++)
	{
		if(fabs(Vx[i]) < eps)
			Vx[i] = 0;
		if(fabs(Vy[i]) < eps)
			Vy[i] = 0;
	}
}

void FluidCube2D::fillParticleInGrid(int x, int y)
{
	int nump = PARTICLEPERGRID;
	int sample = 10000;
	float subSize = 1.0 / nump;
	for(int i = 0; i < nump; i++)
		for(int j = 0; j < nump; j++)
		{
			float x0 = subSize * (rand()%sample) / sample + i * subSize;
			float y0 = subSize * (rand()%sample) / sample + j * subSize;
			particles.push_back(Pos2D(x+x0, y+y0) );
			velosities.push_back(Velo2D(0, 0));
		}
}

void FluidCube2D::report(clock_t simTime)
{
	REPORT(iteration);
	REPORT(dt);

	REPORT(max_vx);
	REPORT(max_vy);
	REPORT(max_p);
	float fluidShrink = 1.0 * fluidNum / originFluid;
	REPORT(fluidShrink);

	REPORT(totalTime);
	LOGSEG;
}

void FluidCube2D::addFlowIn()
{
	if(FLOWINTYPE.compare("TOP") == 0)
	{
		for(int y = NUMGRIDH/8.0*6; y <= NUMGRIDH/8.0*7; y++)
		{
			type[IX(0, y)] = type0[IX(0, y)] = FLUIDIN;
			fillParticleInGrid(1, y);
			Vx[IX(0, y)] = Vx[IX(1, y)] = 1;
			Vy[IX(0, y)] = 0;

			//Blobby
			for(int j = y*GRIDSIZE; j < (y+1)*GRIDSIZE; j++)
				for(int i = 0; i < GRIDSIZE; i++)
				{
					int index = IX2(i, j);
					pixels[index*3]		= 0; //R
					pixels[index*3 + 1]	= 0; //G
					pixels[index*3 + 2]	= 0.7; //B
					pixelType[index] = FLUIDIN;
				}
		}
	}
}

double FluidCube2D::blobbyKernel(double s2)
{
	if(s2 < 1)
		return (1-s2)*(1-s2)*(1-s2);
	else
		return 0;
}