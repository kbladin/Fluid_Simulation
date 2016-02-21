#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <MacGrid.h>
#include <Renderer.h>

#include <string>
#include <sstream>

class Simulator
{
public:
	Simulator();
	~Simulator();
	
private:
	MacGrid* _grid;
	Renderer* _renderer;
};

#endif