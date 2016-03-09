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
	std::unique_ptr<MacGrid> _grid;
	std::unique_ptr<Renderer> _renderer;
	std::unique_ptr<MarkerParticleSet> _particle_set;
};

#endif