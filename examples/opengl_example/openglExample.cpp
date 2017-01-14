#include <FluidSolver.h>
#include <FluidDomain.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>

#include "FluidRendererGL.h"
#include "ApplicationWindowGLFW.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

#include <gl/glew.h>
#include <gl/glfw3.h>

static const int GRID_X_SIZE = 90;
static const int GRID_Y_SIZE = 60;
static const double WORLD_X_SIZE = 1.5;
static const double WORLD_Y_SIZE = 1;
static const double DELTA_X = WORLD_X_SIZE / double(GRID_X_SIZE);
static const double DELTA_Y = WORLD_Y_SIZE / double(GRID_Y_SIZE);

int main(int argc, char const *argv[])
{
	// Initialize
	FluidDomain fluid_domain(GRID_X_SIZE, GRID_Y_SIZE, WORLD_X_SIZE, WORLD_Y_SIZE, 0.005, 0.02);
    FluidSolverMemoryPool mem_pool(fluid_domain);
	FluidSolver fluid_solver(mem_pool);

	// Setup
    fluid_domain.addFluidSource(FluidSource( { 4.0 / GRID_X_SIZE, 0.1, 3.0 / GRID_Y_SIZE, 1 - 4.0 / GRID_Y_SIZE }, DELTA_X, DELTA_Y, 0.0, 0.0, 0.0, 1));
    
	ApplicationWindowGLFW window(512 * WORLD_X_SIZE, 512 * WORLD_Y_SIZE);
	FluidRendererGL renderer(512 * WORLD_X_SIZE, 512 * WORLD_Y_SIZE, WORLD_X_SIZE, WORLD_Y_SIZE);
    FluidInteractionHandler interaction_handler(&fluid_domain, &renderer);
	window.setInteractionHandler(&interaction_handler);
	

    // Prepare simulation
	MyFloat seconds_per_frame = 1 / 60.0;
	
	std::function<void(void)> loop = [&](){
		MyFloat dt;
		// Simulate
		for (MyFloat frame_time = 0; frame_time < seconds_per_frame; frame_time += dt)
		{
			// Calculate dt (for now just set it)
			dt = 0.005;
            dt = CLAMP(dt, 0, seconds_per_frame - frame_time);

			// Update fluid domain (creates new fluid from sources)
			fluid_domain.update(dt);

			// Solve
			//fluid_solver.stepSemiLagrangian(fluid_domain, dt);
			fluid_solver.stepPICFLIP(fluid_domain, dt);
		}
		// Render
		renderer.renderParticles(fluid_domain.markerParticleSet());
        return;
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
