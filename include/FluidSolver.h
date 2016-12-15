#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include <FluidDomain.h>

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <glm/glm.hpp>

class FluidSolver
{
public:
	FluidSolver();
	~FluidSolver();

	void step(FluidDomain& fluid_domain, double dt);
private:
	void enforceDirichlet(MacGrid& mac_grid);
	void pressureSolve(MacGrid& mac_grid, double dt);
	void extendVelocity(MacGrid& mac_grid);

	void advectVelocity(MacGrid& mac_grid, double dt);
	void advectParticles(MarkerParticleSet& particle_set, MacGrid& mac_grid, double dt);
	void advectLevelSet(LevelSet& levelSet, MacGrid& mac_grid, double dt);

	void getAdvectedPositionRK3(
		MacGrid& mac_grid,
		double x_pos,
		double y_pos,
		double dt,
		double* x,
		double* y);
	void getAdvectedPositionForwardEuler(
		MacGrid& mac_grid,
		double x_pos,
		double y_pos,
		double dt,
		double* x,
		double* y);

	// Solver of linear system
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > _cg_solver;
};

#endif
