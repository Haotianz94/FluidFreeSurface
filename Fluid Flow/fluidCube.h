#ifndef _FLUID_H_
#define _FLUID_H_

#include "freeglut.h"
#include <vector>
#include <Eigen/Eigen>

#ifdef SIMULATION_2D

struct Pos
{
	float x;
	float y;
	Pos() {}
	Pos(float _x, float _y): x(_x), y(_y) {}
};

typedef Pos Velo;

class FluidCube2D
{
private:
	float h;
	float hi;
	float h2;
	float dt;
	float visc;

	float *Vx;
	float *Vy;
	GRIDTYPE *type;
	std::vector<Pos> obstacle;
	float max_vx;
	float max_vy;
	
	//Projection using Conjugate Gradient
	Pos dir[4];
	int fluidNum;
	int **neighbor;
	int *neighNoneSolid;
	int *neighAir;
	int *pos2index;
	Eigen::SparseMatrix<double> A;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

	//Advection using BFECC
	float *fai_b;
	float *fai_f;

	//MAC
	std::vector<Pos> particles;
	GRIDTYPE *type0;
	std::vector<int> **invertedList;

public:
	int size;
	float *Vx0;
	float *Vy0;

private:
	void calculateTimeStep();
	void updateParticles();
	void updateGrid();
	float getVelosity(int index, float x, float y, float *u);
	Velo getVelosity(float x, float y, float *vx, float *vy);
	Pos traceParticle(int index, int x, int y, bool backward);

	void vel_step();
	void addForce();
	void diffuseVelosity();
	void advectVelosity();
	void projectVelosity();
	
	void add_source(float *x, float *amount);
	void diffuse(int b, float *u0, float *u, float diffusion);
	void advect(int b, float *u0, float *u, bool backward);
	void swap(float *x0, float *x);
	void set_bnd();

	void draw_dens();
	void draw_velo(int i, int j, float vx, float vy);

	void output(float *u);

public:
	FluidCube2D(float viscosity, float dt);
	~FluidCube2D();
	void simulate();
};

#else

class FluidCube3D
{
private:
	float h;
	float dt;
	float diff;
	float visc;

	float *d;
	float *Vx;
	float *Vy;
	float *Vz;

	float max_d;
	float max_vx;
	float max_vy;
	float max_vz;

public:
	int size;
	float *d0;
	float *Vx0;
	float *Vy0;
	float *Vz0;

private:
	void dens_step(float *amount);
	void addDensity(float *amount);
	void diffuseDensity();
	void advectDensity();
	
	void vel_step(float *amountX, float *amountY, float *amountZ);
	void addVelocity(float *amountX, float *amountY, float *amountZ);
	void diffuseVelosity();
	void advectVelosity();
	void projectVelosity();
	
	void add_source(float *x, float *amount);
	void diffuse(int b, float *u0, float *u, float diffusion);
	void advect(int b, float *u0, float *u);
	void swap(float *x0, float *x);
	void set_bnd(int b, float *x);

	void draw_velo(int i, int j, int k, float vx, float vy, float vz);
	void draw_grid();

	void output(float *u);

public:
	FluidCube3D(float diffusion, float viscosity, float dt);
	~FluidCube3D();
	void simulate(bool idle);
	void draw_dens();
};

#endif

#endif