#ifndef _STDAFX_H_
#define _STDAFX_H_

#include <iostream>

#define SIMULATION_2D 1
//#define CONNECTED 1
//#define OBSTACLE 1
#define OUTPUT 1
//#define GAUSS_SEIDEL 1

#ifdef SIMULATION_2D
	#define _W 80 
	#define _H 80
	#define _L 1.0
 	#define GRIDSIZE 10
	#define VISCOSITY 0.00001
	//#define TIMESTEP 0.05
	#define ITERATION 30
	#define FRAMERATE 32
	#define NUMPERGRID 6
	#define GRAVITY 9.8
	#define MYSCENE CONTAINER
	#define MYRENDER PARTICLE

	#define IX(x, y) ( (x) + (y) * (_W+2) )
	#define BOUNDED(x, y) ( (type[IX(int(x),int(y))] == SOLID || type[IX(int(x)+1,int(y)+1)] == SOLID)? false : true)
	#define DISTANCE(x1, y1, x2, y2) ( sqrtf((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)) )
#else
	#define _N 40
	#define GRIDSIZE 20
	#define DIFFUSION 0.0001
	#define VISCOSITY 0.0001
	#define TIMESTEP 0.01
	#define ITERATION 10
	#define FRAMERATE 32

	#define IX(x, y, z) ((x) + (y)*(_N+2) + (z)*(_N+2)*(_N+2) )
#endif

#define PI 3.14159265
#define LENGTH _N*GRIDSIZE
#define SWAP(x0, x) {float *tmp = x0; x0 = x; x = tmp;}

//system output defines
#define REPORT(X) std::cout << #X << ": " << X << std::endl
#define PRINT(X) std::cout << X << std::endl
#define PRINT_ONELINE(X) std::cout << X

#endif

//Best for vortex street connected
/*
	#define _W 1000 
	#define VISBLEW 1600
	#define _H 50 
	#define GRIDSIZE 2
	#define DIFFUSION 0.01
	#define VISCOSITY 0.01
	#define TIMESTEP 0.01
	#define ITERATION 30
	#define FRAMERATE 32
	#define DRAGSCALE 100
	#define FLOWTIME 10
	#define DENSITY 100
	#define SPEED 10000
	#define OBSTACLEX 30
*/

//Best for vortex street not connected
/*
	#define _W 150
	#define VISBLEW 1500
	#define _H 30
	#define GRIDSIZE 10
	#define DIFFUSION 0.01
	#define VISCOSITY 0.01
	#define TIMESTEP 0.01
	#define ITERATION 10
	#define FRAMERATE 32
	#define DRAGSCALE 100
	#define FLOWTIME 20
	#define DENSITY 100
	#define SPEED 200
*/