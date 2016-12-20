#ifndef RENDERER_H
#define RENDERER_H

#include <fstream>

#include <MacGrid.h>
#include <LevelSet.h>
#include <Canvas.h>
#include <MarkerParticleSet.h>

#include "MathDefinitions.h"

class Renderer
{
public:
	Renderer(MyFloat x_min, MyFloat y_min, MyFloat x_max, MyFloat y_max);
	~Renderer();

	void clearCanvas();
	void renderGridVelocitiesToCanvas(const MacGrid& grid);
	void renderGridCellsToCanvas(const MacGrid& grid);
	void renderLevelSetFunctionValuesToCanvas(const LevelSet& level_set);
	void renderParticlesToCanvas(const MarkerParticleSet& particle_set);

	void writeCanvasToPpm( const char* file_path);
private:
	Canvas _canvas;

	// World coordinates, defines a quad to render
	MyFloat _x_min;
	MyFloat _y_min;
	MyFloat _x_max;
	MyFloat _y_max;
};

#endif
