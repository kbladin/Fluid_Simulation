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
		BBox<MyFloat> area,
		MyFloat delta_x,
		MyFloat delta_y,
		MyFloat x_velocity,
		MyFloat y_velocity,
		MyFloat time_step,
		int max_spawns);
	~FluidSource();
	void update(MarkerParticleSet& particle_set, MyFloat dt);
	inline bool isFinished() { return _max_spawns != -1 && _n_spawns >= _max_spawns; };
	inline void resetSpawns() { _n_spawns = 0; };
private:
	BBox<MyFloat> _area;
	MyFloat _x_velocity;
	MyFloat _y_velocity;
	MyFloat _delta_x;
	MyFloat _delta_y;
	MyFloat _time_step;
	MyFloat _time_since_last;
	int _max_spawns;
	int _n_spawns;
};

class FluidDomain : public GridInterface
{
public:
	FluidDomain(
		int size_x,
		int size_y,
		MyFloat length_x,
		MyFloat length_y,
		MyFloat density,
		MyFloat pic_ratio);
	~FluidDomain();

	void addFluidSource(FluidSource fluid_source);
	void update(MyFloat dt);
	void clearFluidSources();
	void resetParticleSet();
	void setPicRatio(MyFloat pic_ratio);
	
	MacGrid& macGrid();
	LevelSet& levelSet();
	MarkerParticleSet& markerParticleSet();
	const MyFloat density() const;
	const MyFloat picRatio() const;

	void advectParticles(MyFloat dt);
	void advectLevelSet(MyFloat dt);

	void classifyCells(LevelSet& levelSet);
	void classifyCells(MarkerParticleSet& particle_set);
private:
	MyFloat _density;
	MyFloat _pic_ratio;
	MacGrid _mac_grid;
	LevelSet _level_set;
	MarkerParticleSet _particle_set;
	std::vector<FluidSource> _fluid_sources;
};

#endif
