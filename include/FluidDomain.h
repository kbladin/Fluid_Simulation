#ifndef FLUID_DOMAIN_H
#define FLUID_DOMAIN_H

#include <MacGrid.h>
#include <LevelSet.h>
#include <MarkerParticleSet.h>

#include "MathDefinitions.h"

class FluidDomain
{
public:
	FluidDomain(
		int size_x,
		int size_y,
		MyFloat length_x,
		MyFloat length_y,
		MyFloat density = 5);
	~FluidDomain();

	// To be implemented
	void addFluidSource();
	void addExternalForce(MyFloat F_x, MyFloat F_y, MyFloat dt);
	void addExternalAcceleration(MyFloat a_x, MyFloat a_y, MyFloat dt);
	void advectParticles(MyFloat dt);
	void advectLevelSet(MyFloat dt);

	void classifyCells(LevelSet& levelSet);
	void classifyCells(MarkerParticleSet& particle_set);
	
	MacGrid& macGrid();
	LevelSet& levelSet();
	MarkerParticleSet& markerParticleSet();
	const MyFloat density();

private:
	MyFloat _density;
	MacGrid _mac_grid;
	LevelSet _level_set;
	MarkerParticleSet _particle_set;
};

#endif
