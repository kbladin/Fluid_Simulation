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
	Renderer(BBox<MyFloat> area);
	~Renderer();

	void clearCanvas(Canvas& canvas);
	void renderGridVelocitiesToCanvas(const MacGrid& grid, Canvas& canvas);
	void renderGridCellsToCanvas(const MacGrid& grid, Canvas& canvas);
	void renderLevelSetFunctionValuesToCanvas(const LevelSet& level_set, Canvas& canvas);
	void renderParticlesToCanvas(const MarkerParticleSet& particle_set, Canvas& canvas);
	void renderMetaBallsToCanvas(const MarkerParticleSet& particle_set, Canvas& canvas);

	void writeCanvasToPpm( const char* file_path, Canvas& canvas);
private:
	// World coordinates, defines a quad to render
	BBox<MyFloat> _area;
};

#endif
