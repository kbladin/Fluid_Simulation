#ifndef FLUID_INTERACTION_HANDLER
#define FLUID_INTERACTION_HANDLER

#include <gl/glew.h>
#include <gl/glfw3.h>

#include <FluidDomain.h>
#include "FluidRendererGL.h"

class FluidInteractionHandler : public Controller
{
public:
	FluidInteractionHandler(FluidDomain& fluid_domain, FluidRendererGL& renderer);
	~FluidInteractionHandler();
	
  virtual void windowSizeCallback(int width, int height) override;
	virtual void keyCallback(Key key, KeyAction action) override;
	virtual void step(float dt) override;	
private:
	// Data
	FluidDomain& _fluid_domain;
	FluidRendererGL& _renderer;
	float _emitter_radius;
};

#endif
