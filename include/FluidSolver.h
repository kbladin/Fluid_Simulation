#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include <FluidDomain.h>

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <glm/glm.hpp>

#include "MathDefinitions.h"

class FluidSolver
{
public:
	FluidSolver();
	~FluidSolver();

	void step(FluidDomain& fluid_domain, MyFloat dt);
private:
	void enforceDirichlet(MacGrid& mac_grid);
	void pressureSolve(
		MacGrid& mac_grid,
		MarkerParticleSet& particle_set,
		MyFloat density);
	void extendVelocity(MacGrid& mac_grid);

	void advectVelocitySemiLagrangian(MacGrid& mac_grid, MyFloat dt);
	void advectParticles(MarkerParticleSet& particle_set, MacGrid& mac_grid, MyFloat dt);
	void advectLevelSet(LevelSet& levelSet, MacGrid& mac_grid, MyFloat dt);

	void getAdvectedPositionRK3(
		MacGrid& mac_grid,
		MyFloat x_pos,
		MyFloat y_pos,
		MyFloat dt,
		MyFloat* x,
		MyFloat* y);
	void getAdvectedPositionForwardEuler(
		MacGrid& mac_grid,
		MyFloat x_pos,
		MyFloat y_pos,
		MyFloat dt,
		MyFloat* x,
		MyFloat* y);

	// Solver of linear system
	Eigen::ConjugateGradient<Eigen::SparseMatrix<MyFloat> > _cg_solver;
};

#endif
