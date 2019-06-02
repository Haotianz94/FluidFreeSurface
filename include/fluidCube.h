#ifndef _FLUIDCUBE_H_
#define _FLUIDCUBE_H_

#include <Eigen/Eigen>
#include <vector>
#include <ctime>


enum GridType
{
	FLUID,
	AIR,
	SOLID,
	FLUIDIN
};


enum RenderType
{
	VELOSITYX,
	VELOSITYY,
	PRESSURE,
	DIVERGENCE,
	PARTICLE,
	FLUIDGRID,
	BLOBBY
};


class FluidCube
{
protected:
	// cube
	int GRIDSIZE;
	int PARTICLEPERGRID;
	int NUMGRID;
	float GRAVITY;
	float VISCOSITY;
	float h;
	float hi;
	float h2;

	// fluid field
	float *div;
	GridType *type;
	GridType *type0;
	int **neighbor;
	int *neighNoneSolid;
	int *neighAir;
	int *pos2index;
	std::vector<int> **invertedList;

	// Conjugate Gradient
	Eigen::SparseMatrix<double> A;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
	Eigen::VectorXd p;

	// advection using BFECC
	float *fai_b;
	float *fai_f;

	// scene
	std::string SCENETYPE;
	std::string FLOWINTYPE;
	std::string OBSTACLETYPE;
	int originFluid;

	// simulate/render
	RenderType RENDERTYPE;
	int MAXITERATION;
	float FRAMERATE;
	float frameTime;

	// dynamic status
	int fluidNum;
	float max_v;
	float max_p;
	float totalTime;
	int iteration;
	float dt;
	float ctime;

	// output
	bool CREATEBLOBBY;
	bool DEBUGPRINT;

public:
	FluidCube();
	~FluidCube();
	virtual void simulate() = 0;
	virtual void render() = 0;
};

#endif