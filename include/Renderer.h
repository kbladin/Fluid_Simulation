#ifndef RENDERER_H
#define RENDERER_H

#include <fstream>

#include <MacGrid.h>
#include <LevelSet.h>
#include <Canvas.h>
#include <MarkerParticleSet.h>

class Renderer
{
public:
	Renderer(double x_min, double y_min, double x_max, double y_max);
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
	double _x_min;
	double _y_min;
	double _x_max;
	double _y_max;
};

#endif
