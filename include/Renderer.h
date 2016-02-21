#ifndef RENDERER_H
#define RENDERER_H

#include <fstream>

#include <MacGrid.h>
#include <Canvas.h>

class Renderer
{
public:
	Renderer(MacGrid* grid);
	~Renderer();

	void renderCanvas();

	void writeToPpm( const char* file_path, float min_val, float max_val);
	void writeCanvasToPpm( const char* file_path);
private:
	MacGrid* _grid;
	Canvas* _canvas;
};

#endif