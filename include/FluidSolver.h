#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include <FluidDomain.h>
#include <OdeSolver.h>

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include "MathDefinitions.h"

#include <memory>

class FluidSolverMemoryPool : public GridInterface
{
public:
	FluidSolverMemoryPool(
		int size_x,
		int size_y,
		MyFloat delta_x,
		MyFloat delta_y);
	FluidSolverMemoryPool(const FluidDomain& fluid_domain);
	FluidSolverMemoryPool(const FluidSolverMemoryPool& other);
	~FluidSolverMemoryPool();

	void swapValidMaskBuffer();

	// Used for pressure solve
    Grid<int> fluid_indices;
    Grid<unsigned char> n_particles;

    // Used for velocity extension
    std::unique_ptr< Grid<unsigned char> > valid_mask_x_front_buffer;
    std::unique_ptr< Grid<unsigned char> > valid_mask_x_back_buffer;
    std::unique_ptr< Grid<unsigned char> > valid_mask_y_front_buffer;
    std::unique_ptr< Grid<unsigned char> > valid_mask_y_back_buffer;
    
    // Used for PIC solve
    Grid<MyFloat> vel_x_sum;
    Grid<MyFloat> vel_y_sum;
    Grid<MyFloat> weight_vel_x_sum;
    Grid<MyFloat> weight_vel_y_sum;
};

class FluidSolver
{
public:
	FluidSolver(FluidSolverMemoryPool mem_pool);
	~FluidSolver();

	void stepSemiLagrangian(FluidDomain& fluid_domain, MyFloat dt);
	void stepPIC(FluidDomain& fluid_domain, MyFloat dt);
	void stepFLIP(FluidDomain& fluid_domain, MyFloat dt);
	void stepPICFLIP(FluidDomain& fluid_domain, MyFloat dt);
private:
	bool validate(const FluidDomain& fluid_domain);
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
		MyFloat density,
		MyFloat dt);
	void extendVelocityIndividual(MacGrid& mac_grid, int n_iterations);
	void extendVelocityAvarageing(MacGrid& mac_grid, int n_iterations);

	void advectVelocitySemiLagrangian(MacGrid& mac_grid, MyFloat dt);
	void advectParticlesWithGrid(
		MarkerParticleSet& particle_set,
		MacGrid& mac_grid,
		MyFloat dt);
	void advectParticles(
		MarkerParticleSet& particle_set,
		MacGrid& mac_grid,
		MyFloat dt);
	void advectLevelSet(LevelSet& levelSet, MacGrid& mac_grid, MyFloat dt);

	void getAdvectedPosition(
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
    // Laplacian matrix
    Eigen::SparseMatrix<MyFloat> _A;
    FluidSolverMemoryPool _mem_pool;

    // Struct needed by the solver
   	struct Vec2
    {
    	MyFloat x;
    	MyFloat y;

    	Vec2 operator+(const MyFloat& h)
    	{
    		return {x + h, y + h};
    	}
    	Vec2 operator+(const Vec2& v)
    	{
    		return {x + v.x, y + v.y};
    	}
    	Vec2 operator*(const MyFloat& s)
    	{
    		return {x * s, y * s};
    	}
    	Vec2 operator/(const MyFloat& s)
    	{
    		return {x / s, y / s};
    	}
    };

    //EulerExplicit<Vec2, Vec2> _ode_solver;
	RK3<Vec2, Vec2, MyFloat> _ode_solver;
	//RK4<Vec2, Vec2> _ode_solver;
};

#endif
