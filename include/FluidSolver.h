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
	FluidSolver(int size_x, int size_y, MyFloat delta_x, MyFloat delta_y);
	~FluidSolver();

	void stepSemiLagrangian(FluidDomain& fluid_domain, MyFloat dt);
	void stepPIC(FluidDomain& fluid_domain, MyFloat dt);
	void stepFLIP(FluidDomain& fluid_domain, MyFloat dt);
	void stepPICFLIP(FluidDomain& fluid_domain, MyFloat dt, MyFloat pic_ratio);
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
	void extendVelocity(MacGrid& mac_grid, int n_iterations);

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

	void transferVelocityToGridGather(
		MarkerParticleSet& particle_set,
		MacGrid& mac_grid);
	void transferVelocityToGridSpread(
		MarkerParticleSet& particle_set,
		MacGrid& mac_grid);
	void transferVelocityToParticlesPIC(
		MacGrid& mac_grid,
		MarkerParticleSet& particle_set);
	void transferVelocityToParticlesFLIP(
		MacGrid& mac_grid,
		MarkerParticleSet& particle_set);
	void transferVelocityToParticlesPICFLIP(
		MacGrid& mac_grid,
		MarkerParticleSet& particle_set,
        MyFloat pic_ratio);

	// Solver of linear system
	Eigen::ConjugateGradient<Eigen::SparseMatrix<MyFloat> > _cg_solver;

	// Cached
	Grid<int> _fluid_indices;
	//Grid<int> _n_particles;
	Grid<char> _valid_mask_x;
	Grid<char> _valid_mask_x_back_buffer;
    Grid<char> _valid_mask_y;
	Grid<char> _valid_mask_y_back_buffer;
    Eigen::SparseMatrix<MyFloat> A;

    // Used for PIC solve
    SizedGrid<MyFloat> _vel_x_sum;
	SizedGrid<MyFloat> _vel_y_sum;
	SizedGrid<MyFloat> _weight_vel_x_sum;
	SizedGrid<MyFloat> _weight_vel_y_sum;

	const int _SIZE_X, _SIZE_Y;
	const MyFloat _DELTA_X, _DELTA_Y;
};

#endif
