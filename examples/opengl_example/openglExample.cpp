#include <FluidSolver.h>
#include <FluidDomain.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>

#include "FluidRendererGL.h"
#include "FluidInteractionHandler.h"

#include "sge/window/application_window_glfw.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

#include <gl/glew.h>
#include <gl/glfw3.h>

static const int GRID_X_SIZE = 60;
static const int GRID_Y_SIZE = 40;
static const double WORLD_X_SIZE = 1.5;
static const double WORLD_Y_SIZE = 1;

int main(int argc, char const *argv[])
{
  int window_width = 512 * WORLD_X_SIZE;
  int window_height = 512 * WORLD_Y_SIZE;
  sge::window::ApplicationWindowGLFW window("PIC / FLIP Fluid Simulation",window_width, window_height);
	FluidRendererGL renderer(GRID_X_SIZE, GRID_Y_SIZE, WORLD_X_SIZE, WORLD_Y_SIZE, window_width, window_height);
  
  FluidInteractionHandler interaction_handler(renderer);
  SphericalController controller(renderer.camera());
  WindowSizeController window_controller(renderer.renderer());
  
	window.addController(interaction_handler);
	window.addController(controller);
	window.addController(window_controller);
	
	std::function<void(double)> loop = [&](double dt)
  {
    renderer.update(dt);
  };

	try
  {
    window.run(loop);
  }
	catch (const std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
