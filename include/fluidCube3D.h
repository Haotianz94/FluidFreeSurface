#ifndef _FLUIDCUBE3D_H_
#define _FLUIDCUBE3D_H_

#include "fluidCube.h"
#include "quadmesh.h"

struct Pos3D
{
	float x;
	float y;
	float z;
	Pos3D() {}
	Pos3D(float _x, float _y, float _z): x(_x), y(_y), z(_z){}
};
typedef Pos3D Velo3D;


class FluidCube3D: public FluidCube
{
private:
	// cube
	int NUMGRIDX;
	int NUMGRIDY;
	int NUMGRIDZ;
	Eigen::Vector3i dir[6];

	// fluid field
	float *Vx;
	float *Vy;
	float *Vz;
	float *Vx0;
	float *Vy0;
	float *Vz0;
	std::vector<Pos3D> particles;
	std::vector<Velo3D> velosities;

	// status
	float max_vx;
	float max_vy;
	float max_vz;
	
	//Blobby
	Eigen::Vector3i dir2[27];
	char obj_dir[100];
	int blobbyFrameStride;

	//Extrapolate
	int *layer;

	//volcano
	QuadMesh* volcano;
	QuadMesh* ground;

private:
	// simulate
	bool calculateTimeStep();
	void updateParticles();
	void updateGrid();
	void set_bnd();
	void vel_step();
	void extrapolate();

	// vel_step
	void addForce();
	void diffuseVelosity();
	void advectVelosity();
	void projectVelosity();
	void errorRemove();

	// help function
	void diffuse(int b, float *u0, float *u, float diffusion);
	void advect(int b, float *u0, float *u, bool backward);
	float getVelosity(int index, float x, float y, float z, float *u);
	Velo3D getVelosity(float x, float y, float z, float *vx, float *vy, float *vz);
	Pos3D traceParticle(int index, int x, int y, int z, bool backward);
	void fillParticleInGrid(int x, int y, int z);

	// scene
	void addFlowIn();
	void initFluid();
	void initSolid();

	// output
	void output(float *u);
	void report(clock_t);

	// blobby
	void createBlobbySurface();
	double blobbyKernel(double s2);

public:
	FluidCube3D();
	~FluidCube3D();
	void simulate();
	void render();
	void createBlobby();
};

#endif