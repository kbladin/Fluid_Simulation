#ifndef FLUID_DOMAIN_H
#define FLUID_DOMAIN_H

#include <MacGrid.h>
#include <LevelSet.h>
#include <MarkerParticleSet.h>

class FluidDomain
{
public:
	FluidDomain(int size_x, int size_y, double length_x, double length_y);
	~FluidDomain();

	// To be implemented
	void addFluidSource();
	void addExternalForce(double F_x, double F_y, double dt);
	void advectParticles(double dt);
	void advectLevelSet(double dt);

	void classifyCells(LevelSet& levelSet);
	void classifyCells(MarkerParticleSet& particle_set);
	
	MacGrid& macGrid();
	LevelSet& levelSet();
	MarkerParticleSet& markerParticleSet();

private:
	MacGrid _mac_grid;
	LevelSet _level_set;
	MarkerParticleSet _particle_set;
};

#endif
