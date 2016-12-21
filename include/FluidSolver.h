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
	FluidSolver(int x_size, int y_size);
	~FluidSolver();

	void step(FluidDomain& fluid_domain, MyFloat dt);
private:
	void addExternalForce(
		FluidDomain& fluid_domain,
		MyFloat F_x,
		MyFloat F_y,
		MyFloat dt);
	void addExternalAcceleration(
		MacGrid& mac_grid,
		MyFloat a_x,
		MyFloat a_y,
		MyFloat dt);
	void enforceDirichlet(MacGrid& mac_grid);
	void pressureSolve(
		MacGrid& mac_grid,
		MarkerParticleSet& particle_set,
		MyFloat density);
	void extendVelocity(MacGrid& mac_grid);

	void advectVelocitySemiLagrangian(MacGrid& mac_grid, MyFloat dt);
	void advectParticles(
		MarkerParticleSet& particle_set,
		MacGrid& mac_grid,
		MyFloat dt);
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

	// Cached
	Grid<int> _fluid_indices;
	//Grid<int> _n_particles;
	Grid<int> _valid_mask;
	Grid<int> _valid_mask_back_buffer;
	const int _SIZE_X, _SIZE_Y;
};

#endif
