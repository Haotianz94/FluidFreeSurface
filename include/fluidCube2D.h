#ifndef _FLUIDCUBE2D_H_
#define _FLUIDCUBE2D_H_

#include "fluidCube.h"


struct Pos2D
{
	float x;
	float y;
	Pos2D() {}
	Pos2D(float _x, float _y): x(_x), y(_y) {}
};
typedef Pos2D Velo2D;


class FluidCube2D: public FluidCube
{
private:
	// cube
	int NUMGRIDH;
	int NUMGRIDW;
	Eigen::Vector2i dir[4];

	// fluid field
	float *Vx;
	float *Vy;
	float *Vx0;
	float *Vy0;
	std::vector<Pos2D> particles;
	std::vector<Velo2D> velosities;

	// status
	float max_vx;
	float max_vy;
	
	// blobby
	float* pixels;
	GridType* pixelType;
	Eigen::Vector2i dir2[9];

private:
	// simulate
	bool calculateTimeStep();
	void updateParticles();
	void updateGrid();
	void set_bnd();
	void vel_step();

	// vel_step
	void addForce();
	void diffuseVelosity();
	void advectVelosity();
	void projectVelosity();
	void errorRemove();

	// help function
	void diffuse(int b, float *u0, float *u, float diffusion);
	void advect(int b, float *u0, float *u, bool backward);
	float getVelosity(int index, float x, float y, float *u);
	Velo2D getVelosity(float x, float y, float *vx, float *vy);
	Pos2D traceParticle(int index, int x, int y, bool backward);
	void fillParticleInGrid(int x, int y);
	
	// scene
	void addFlowIn();

	// output
	void draw_velo(int i, int j, float vx, float vy);
	void output(float *u);
	void report(clock_t);

	// blobby
	double blobbyKernel(double s2);

public:
	FluidCube2D();
	~FluidCube2D();
	void simulate();
	void render();
};

#endif