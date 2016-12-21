#ifndef FLUID_DOMAIN_H
#define FLUID_DOMAIN_H

#include <MacGrid.h>
#include <LevelSet.h>
#include <MarkerParticleSet.h>

#include "MathDefinitions.h"

class FluidSource
{
public:
	FluidSource(
		MyFloat x_min,
		MyFloat x_max,
		MyFloat y_min,
		MyFloat y_max,
		MyFloat time_step,
		int max_creations);
	~FluidSource();
	void update(MarkerParticleSet& particle_set, MyFloat dt);
	bool isFinished();
private:
	MyFloat _x_min;
	MyFloat _x_max;
	MyFloat _y_min;
	MyFloat _y_max;
	MyFloat _time_step;
	MyFloat _time_since_last;
	int _max_creations;
	int _n_creations;
};

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
	void addFluidSource(FluidSource fluid_source);
	void update(MyFloat dt);
	
	MacGrid& macGrid();
	LevelSet& levelSet();
	MarkerParticleSet& markerParticleSet();
	const MyFloat density();

private:
	void advectParticles(MyFloat dt);
	void advectLevelSet(MyFloat dt);

	void classifyCells(LevelSet& levelSet);
	void classifyCells(MarkerParticleSet& particle_set);
	
	MyFloat _density;
	MacGrid _mac_grid;
	LevelSet _level_set;
	MarkerParticleSet _particle_set;
	std::vector<FluidSource> _fluid_sources;
};

#endif
