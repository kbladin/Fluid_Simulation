#ifndef FLUID_INTERACTION_HANDLER
#define FLUID_INTERACTION_HANDLER

#include <gl/glew.h>
#include <gl/glfw3.h>

#include <FluidDomain.h>
#include "FluidRendererGL.h"

class FluidInteractionHandler
{
public:
	FluidInteractionHandler(FluidDomain* fluid_domain, FluidRendererGL* renderer);
	~FluidInteractionHandler();

	// NDC coordinates
	void mousePosCallback(double x, double y);
	void mouseButtonCallback(int button, int action, int mods);
	void mouseScrollCallback(double dx, double dy);
	void keyCallback(int key, int scancode, int action, int mods);
  	void windowSizeCallback(int width, int height);
	
private:
	// Data
	FluidDomain* _fluid_domain;
	FluidRendererGL* _renderer;
	double _mouse_x, _mouse_y;
	float _emitter_radius;
};

#endif
