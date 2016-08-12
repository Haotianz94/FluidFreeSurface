#include "stdafx.h"
#include "fluidCube.h"
#include <memory.h>
#include <math.h>
#include <Eigen\Eigen>
//#include "vl\VLf.h"
//#include "vl\Solve.h"

#ifdef SIMULATION_2D

FluidCube2D::FluidCube2D(float viscosity, float fr, SCENETYPE sc)
{
	size = (_W+2) * (_H+2);
	h = _L / _H;
	h2 = h * h;
	hi = 1 / h;
	visc = viscosity;
	frameTime = 1.0 / fr;
	scene = sc;
	ctime = 0;

	max_vx = 0;
	max_vy = 0;

	Vx = new float [size]; 
	Vy = new float [size]; 
	Vx0 = new float [size]; 
	Vy0 = new float [size]; 
	type = new GRIDTYPE [size];

	//Advection using BFECC
	fai_b = new float [size];
	fai_f = new float [size];

	//Projection using Conjugate Gradient
	dir[0] = Eigen::Vector2i(0, -1);
	dir[1] = Eigen::Vector2i(-1, 0);
	dir[2] = Eigen::Vector2i(1, 0);
	dir[3] = Eigen::Vector2i(0, 1);
	pos2index = new int [size]; 
	neighNoneSolid = new int [size];
	neighAir = new int [size];
	neighbor = new int* [size];
	for(int i = 0; i < size; i ++)
		neighbor[i] = new int[4];

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
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
			type[IX(x,y)] = AIR;

	//init fluid
	//cube fall
	if(scene == CUBEFALL)
	{
		for(int y = _H/4.0; y <= _H/2.0; y++)
			for(int x = _W/3.0; x <= _W/3.0*2; x++)
				fillParticleInGrid(x, y);
	}
	//sphere fall
	else if(scene == SPHEREFALL)
	{
		int cx = _W/2;
		int cy = _H/2;
		float R = _W/6;
		for(int y = 1; y <= _H; y++)
			for(int x = 1; x <= _W; x++)
				fillParticleInGrid(x, y);
	}
	//contain bottom
	else if(scene == CONTAINER)
	{
		
		for(int y = _H/3.0; y <= _H/2.0; y++)
			for(int x = _W/3.0; x <= _W/3.0*2; x++)
				fillParticleInGrid(x, y);
		

		for(int y = 1; y <= _H/4.0; y++)
			for(int x = 1; x <= _W; x++)
				fillParticleInGrid(x, y);
	}
	//dam break
	else if(scene == DAMBREAK)
	{
		for(int y = 1; y <= _H/3.0*2; y++)
			for(int x = 1; x <= _W/4.0; x++)
				fillParticleInGrid(x, y);
	}

#ifdef OBSTACLE
	int cx = _W / 2.0;
	int cy = _H / 3.0;
	int R = _H * 0.1;
	fluidNum = 0;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(DISTANCE(x,y,cx,cy) <= R)
				type[IX(x, y)] = SOLID;
		}
	/*for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
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
	memset(Vx0, 0, sizeof(float) * size);
	memset(Vy0, 0, sizeof(float) * size);

	memset(fai_b, 0, sizeof(float) * size);
	memset(fai_f, 0, sizeof(float) * size);
}

FluidCube2D::~FluidCube2D()
{
	delete [] Vx;
	delete [] Vy;
	delete [] Vx0;
	delete [] Vy0;
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

}

void FluidCube2D::vel_step()
{
	/*
	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	diffuseVelosity();
	set_bnd();
	*/

	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	advectVelosity();
	set_bnd();

	addForce();
	set_bnd();

	projectVelosity();
	set_bnd();

	errorRemove();
}

void FluidCube2D::addForce()
{
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
			if(type[IX(x, y)] == FLUID)
				Vy[IX(x, y)] -= dt * GRAVITY;
}

void FluidCube2D::diffuseVelosity()
{
	diffuse(1, Vx0, Vx, visc);
	diffuse(2, Vy0, Vy, visc);
}

void FluidCube2D::advectVelosity()
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

	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
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
		for(int y = 1; y <= _H; y++)
			for(int x = 1; x <= _W; x++)
				if(type[IX(x, y)] == FLUID)
					p[IX(x, y)] = (div[IX(x,y)] + p[IX(x,y)] + p[IX(x+1,y)] + p[IX(x,y)] + p[IX(x,y+1)]) / 4;
	}
 
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
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
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;
			div[IX(x, y)] = (Vx[IX(x+1,y)]-Vx[IX(x,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y)]);
		}
	output(div);
	*/
	Eigen::VectorXd b(fluidNum);
	int index = 0;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			b[index++] =  -h * (Vx[IX(x+1,y)]-Vx[IX(x,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y)]);
		}
	
	//std::cout<<b<<std::endl;

	p.resize(fluidNum);
	p = solver.solve(b);
	
	//for(int i = 100; i < 110; i++)
	//	cout<<p[i]<<' ';
	//cout<<endl;
	/*
	std::cout<<"bbbbbbb"<<std::endl;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			std::cout<<b[pos2index[IX(x,y)]]<<std::endl;
			if(x == _W)
				std::cout<<std::endl;
		}
	
	std::cout<<"ppppppp"<<std::endl;
	for(int y = 1; y <= _H; y++)
	{
		float lastp = -1;
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			//std::cout<<p[pos2index[IX(x,10)]]<<std::endl;
			//if(x == _W)
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

	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			double p1, p2;
			p2 = p[pos2index[IX(x,y)]];
			//Vx
			if(type[IX(x-1, y)] == AIR)
				p1 = 0;
			else if(type[IX(x-1, y)] == FLUID)
				p1 = p[pos2index[IX(x-1,y)]];
			else
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
		for(int y = 1; y <= _H; y++)
			for(int x = 1; x <= _W; x++)
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
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
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
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;
		
			Pos pos = traceParticle(b, x, y, backward);
		
			u[IX(x, y)] = getVelosity(b, pos.x, pos.y, u0);

			/*if(x0 < 0.5)
				x0 = 0.5;
			else if(x0 > _W + 0.5)
				x0 = _W + 0.5;

			if(y0 < 0.5)
				y0 = 0.5;
			else if(y0 > _H + 0.5)
				y0 = _H + 0.5;

			int i0 = int(x0), i1 = i0 + 1;
			int j0 = int(y0), j1 = j0 + 1;
			float s1 = x0 - i0, s0 = 1 - s1;
			float t1 = y0 - j0, t0 = 1 - t1;

			//if trace into solid, up it unchanged
			int nonFluidNum = 0;
			if(type[IX(i0, j0)] != FLUID)
				nonFluidNum ++;
			if(type[IX(i0, j1)] != FLUID)
				nonFluidNum ++;
			if(type[IX(i1, j0)] != FLUID)
				nonFluidNum ++;
			if(type[IX(i1, j1)] != FLUID)
				nonFluidNum ++;
			if(nonFluidNum >= 3)
			{
				u[IX(x,y)] = u0[IX(x,y)];
				continue;
			}

			u[IX(x,y)] = s0 * (t0*u0[IX(i0,j0)] + t1*u0[IX(i0,j1)]) +
						 s1 * (t0*u0[IX(i1,j0)] + t1*u0[IX(i1,j1)]);*/

		}
}

void FluidCube2D::set_bnd()
{
	//try to use free-slip condition
	/*
	for(int y = 1; y <= _H; y++)
	{
		Vx[IX(1, y)] = -Vx[IX(2, y)];
		Vx[IX(_W+1, y)] = -Vx[IX(_W, y)];
	}
	for(int x = 1; x <= _W; x++)
	{
		Vy[IX(x, 1)] = -Vy[IX(x, 2)];
		Vy[IX(x, _H+1)] = -Vy[IX(x, _H)];
	}
	*/
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			if(type[IX(x-1, y)] == SOLID)
				Vx[IX(x, y)] = 0;
			if(type[IX(x, y-1)] == SOLID)
				Vy[IX(x, y)] = 0;

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
	while(true)
	{
		bool draw = calculateTimeStep();
	
		updateParticles();
	
		updateGrid();
	
		set_bnd();
	
		vel_step();

		if(draw)
			draw_dens();
	}
}

void FluidCube2D::output(float *u)
{
#ifdef OUTPUT 
	for(int y = 10; y <= 15; y++)
		for(int x = 5; x <= 10; x++)
		{
			std::cout<<u[IX(x, y)]<<' ';
			if(x == 10)
				std::cout<<std::endl;
		}
	PRINT("=================================\n");
#endif
}

void FluidCube2D::draw_dens()
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

	REPORT(max_vx);
	REPORT(max_vy);

	float max_p = -9999;
	for(int i = 0; i < fluidNum; i++)
		if(p[i] > max_p)
			max_p = p[i];
	REPORT(max_p);
	REPORT(fluidNum);
	
	/*
	float *p = fai_f;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
			if(p[IX(x, y)] > max_p)
				max_p = p[IX(x, y)];
	REPORT(max_p);
	*/

	/*
	//calculate divergence
	float *div = fai_b;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;
			div[IX(x, y)] = (Vx[IX(x+1,y)]-Vx[IX(x,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y)]);
		}
	*/

	for(int i = 0; i <= _W+1; i++)
		for(int j = 0; j <= _H+1; j++)
		{
			int x = i;
			int y = j;
			float color;

			if(type[IX(x, y)] == SOLID)
				glColor3f(0, 0.5, 0);
			
			//else if(type[IX(x, y)] == FLUID)
				//glColor3f(0, 0, 0.7);
				//glColor3f(p[pos2index[IX(x,y)]]/max_p, 0, 0);
				
				//show velosity
				/*
				if(Vy[IX(x, y)] >= 0)
					glColor3f(Vy[IX(x, y)]/max_vy, 0, 0);
				else
					glColor3f(0, -Vy[IX(x, y)]/max_vy, 0);
				*/
				//show diverence
				/*
				if(fabs(div[IX(x, y)]) > 1e-5)
					glColor3f(1, 0, 0);
				else
					glColor3f(0, 0, 0.7);
				*/
			else
				glColor3f(0.5, 0.5, 0.5);
				//vorticity
				//float w = 0.5 * (Vy[IX(x+1, y)] - Vy[IX(x-1, y)]);
				//		  - 0.5 * (Vx[IX(x, y+1)] - Vx[IX(x, y-1)]);
			
			glBegin(GL_QUADS);
			glVertex2f(i*GRIDSIZE, j*GRIDSIZE);
			glVertex2f((i+1)*GRIDSIZE, j*GRIDSIZE);
			glVertex2f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE);
			glVertex2f(i*GRIDSIZE, (j+1)*GRIDSIZE);
			glEnd();

			//if(GRIDSIZE >= 10 && type[IX(x, y)] == FLUID)
			//	draw_velo(i, j, Vx[IX(x, y)], Vy[IX(x, y)]);
		}

	//draw particles

	glColor3f(0, 0, 0.7);
	glBegin(GL_POINTS);
	for(unsigned i = 0; i < particles.size(); i++)
	{
		glVertex2f(particles[i].x*GRIDSIZE, particles[i].y*GRIDSIZE);
	}
	glEnd();


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
	//
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
	REPORT(dt);
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
		if(x1 < 1 || x1 >= _W+1 || y1 < 1 || y1 >= _H+1)
		{
			std::cout<<"Particle out of bound"<<std::endl;
			REPORT(x1);
			REPORT(y1);
			//system("pause");
		}
		if(x1 < 1)
			x1 = 1;
		if(x1 >= _W+1)
			x1 = _W+0.999;
		if(y1 < 1)
			y1 = 1;
		if(y1 >= _H+1)
			y1 = _H+0.999;

		if(type[IX(int(x1), int(y1))] != FLUID)
		{
			particles[i] = Pos(x1, y1);
			continue;
		}

		float vx1 = getVelosity(1, x1, y1, Vx);
		float vy1 = getVelosity(2, x1, y1, Vy);
		//if particle out of boundary???
		if(x1 < 1 || x1 >= _W+1 || y1 < 1 || y1 >= _H+1)
		{
			std::cout<<"Particle out of bound"<<std::endl;
			REPORT(x1);
			REPORT(y1);
			//system("pause");
		}

		x1 = x0 + dt * 0.5 * (vx0 + vx1) * hi;
		y1 = y0 + dt * 0.5 * (vy0 + vy1) * hi;
		if(x1 < 1)
			x1 = 1;
		if(x1 >= _W+1)
			x1 = _W+0.999;
		if(y1 < 1)
			y1 = 1;
		if(y1 >= _H+1)
			y1 = _H+0.999;
		particles[i].x = x1;
		particles[i].y = y1;
	}
}

void FluidCube2D::updateGrid()
{
	//swap 
	GRIDTYPE *tmp = type;
	type = type0;
	type0 = tmp;

	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type0[IX(x, y)] != SOLID)
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
	for(int i = 0; i < size; i++)
		pos2index[i] = -1;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
			if(type[IX(x,y)] == FLUID)
				pos2index[IX(x, y)] = fluidNum++;

	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
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
					if(type[IX(xx,yy)] != SOLID)
						neighNoneSolid[IX(x, y)] ++;
					if(type[IX(xx,yy)] == AIR)
						neighAir[IX(x, y)] ++;
				}
			}

			//if(type[IX(x, y)] == AIR)
			//	Vx[IX(x, y)] = Vy[IX(x, y)] = 0;

			if(type0[IX(x, y)] == AIR && type[IX(x, y)] == FLUID)
			{
				float vx = 0, vy = 0;
				int nx = 0, ny = 0;
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
						nx ++;
						vx += getVelosity(1, x0, y0, Vx);
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
						nx ++;
						vx += getVelosity(1, x0, y0, Vx);
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
				//if(nx > 0)
				//	Vx[IX(x, y)] = vx / nx;
				//else
				{
					//for simple only take the nearest particle in the cell
					Vx[IX(x, y)] = getVelosity(1, particles[pid].x, particles[pid].y, Vx);
				}

				dist = 100;
				//Vy
				list = invertedList[IX(x, y)];
				for(unsigned i = 0; i < list->size(); i++)
				{
					float x0 = particles[list->at(i)].x;
					float y0 = particles[list->at(i)].y;
					if( y0 <= y + 0.5)
					{
						ny ++;
						vy += getVelosity(2, x0, y0, Vy);
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
						ny ++;
						vy += getVelosity(2, x0, y0, Vy);
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
				//if(ny > 0)
				//	Vy[IX(x, y)] = vy / ny;
				//else
				{
					//for simple only take the nearest particle in the cell
					Vy[IX(x, y)] = getVelosity(2, particles[pid].x, particles[pid].y, Vy);
				}
			}
		}


	//init Matrix
	A = Eigen::SparseMatrix<double>(fluidNum, fluidNum);         // default is column major
	A.reserve(Eigen::VectorXi::Constant(fluidNum, 5));
	int index = 0;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
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

	if(x < 0 || x >= _W+1 || y < 0 || y >= _H+1)
	{
		std::cout<<"Get velosity out of bound"<<std::endl;
		REPORT(x);
		REPORT(y);
		//system("pause");
	}

	int i0 = int(x), i1 = i0 + 1;
	int j0 = int(y), j1 = j0 + 1;
	float s1 = x - i0, s0 = 1 - s1;
	float t1 = y - j0, t0 = 1 - t1;

	return s0 * (t0*u[IX(i0,j0)] + t1*u[IX(i0,j1)]) +
		   s1 * (t0*u[IX(i1,j0)] + t1*u[IX(i1,j1)]);
}

Velo FluidCube2D::getVelosity(float x, float y, float *vx, float *vy)
{
	return Velo(getVelosity(1,x,y,vx), getVelosity(2,x,y,vy));
}

Pos FluidCube2D::traceParticle(int index, int x, int y, bool backward)
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
	Velo v0 = getVelosity(x0, y0, Vx0, Vy0);
	float t = (backward)? -dt : dt;
	//
	//Pos p = Pos(x0 + v0.x*t*hi, y0 + v0.y*t*hi);
	Velo v1 = getVelosity(x0 + v0.x*t*hi, y0 + v0.y*t*hi, Vx0, Vy0);
	return Pos(x0 + 0.5*t*(v0.x+v1.x)*hi, y0 + 0.5*t*(v0.y+v1.y)*hi);
}

void FluidCube2D::errorRemove()
{
	double eps = 1e-12;

	for(int i = 0; i < size; i++)
	{
		if(fabs(Vx[i]) < eps)
			Vx[i] = 0;
		if(fabs(Vy[i]) < eps)
			Vy[i] = 0;
	}
}

void FluidCube2D::fillParticleInGrid(int x, int y)
{
	int nump = NUMPERGRID;
	float step = 1.0 / nump;
	for(int i = 0; i < nump; i++)
		for(int j = 0; j < nump; j++)
			particles.push_back(Pos((x+step*i), (y+step*j)));
}

#endif