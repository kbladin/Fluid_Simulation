#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <MacGrid.h>
#include <Renderer.h>
#include <MarkerParticleSet.h>

#include <string>
#include <sstream>

class Simulator
{
public:
	Simulator();
	~Simulator();
	
	void advectParticles(double dt);
	void updateCellTypesWithParticles();
private:
	MacGrid* _grid;
	Renderer* _renderer;
	MarkerParticleSet* _particle_set;
};

#endif