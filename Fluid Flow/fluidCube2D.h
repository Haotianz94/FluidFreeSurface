#ifndef _FLUIDCUBE2D_H_
#define _FLUIDCUBE2D_H_

#include "freeglut.h"
#include <vector>
#include <Eigen/Eigen>

enum GRIDTYPE
{
	FLUID,
	AIR,
	SOLID
};

struct Pos
{
	float x;
	float y;
	Pos() {}
	Pos(float _x, float _y): x(_x), y(_y) {}
};
typedef Pos Velo;

enum SCENETYPE
{
	CUBEFALL,
	SPHEREFALL,
	CONTAINER,
	DAMBREAK,
	DOUBLEDAM
};

enum RENDERTYPE
{
	VELOSITYX,
	VELOSITYY,
	PRESSURE,
	DIVERGENCE,
	PARTICLE,
	FLUIDGRID
};

class FluidCube2D
{
private:
	float h;
	float hi;
	float h2;
	float dt;
	float visc;

	int size;
	float *Vx;
	float *Vy;
	float *Vx0;
	float *Vy0;
	float *div;
	GRIDTYPE *type;
	std::vector<Pos> obstacle;
	float max_vx;
	float max_vy;
	float max_p;
	
	//Projection using Conjugate Gradient
	Eigen::Vector2i dir[4];
	int fluidNum;
	int originFluid;
	int **neighbor;
	int *neighNoneSolid;
	int *neighAir;
	int *pos2index;
	Eigen::SparseMatrix<double> A;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	Eigen::VectorXd p;

	//Advection using BFECC
	float *fai_b;
	float *fai_f;

	//MAC
	std::vector<Pos> particles;
	GRIDTYPE *type0;
	std::vector<int> **invertedList;
	
	SCENETYPE scene;
	RENDERTYPE renderType;
	float totalTime;
	int iteration;
	float ctime;
	float frameTime;

private:
	bool calculateTimeStep();
	void updateParticles();
	void updateGrid();
	float getVelosity(int index, float x, float y, float *u);
	Velo getVelosity(float x, float y, float *vx, float *vy);
	Pos traceParticle(int index, int x, int y, bool backward);

	void errorRemove();
	void fillParticleInGrid(int x, int y);

	void vel_step();
	void addForce();
	void diffuseVelosity();
	void advectVelosity();
	void projectVelosity();
	
	void diffuse(int b, float *u0, float *u, float diffusion);
	void advect(int b, float *u0, float *u, bool backward);
	void swap(float *x0, float *x);
	void set_bnd();

	void draw_velo(int i, int j, float vx, float vy);

	void output(float *u);
	void report();

public:
	FluidCube2D(float viscosity, float fr, SCENETYPE sc = CONTAINER, RENDERTYPE rt = PARTICLE);
	~FluidCube2D();
	void simulate();
	void render();
};

#endif