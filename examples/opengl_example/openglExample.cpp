#include <FluidSolver.h>
#include <FluidDomain.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>

#include "FluidRendererGL.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

#include <gl/glew.h>
#include <gl/glfw3.h>

static const int GRID_X_SIZE = 60;
static const int GRID_Y_SIZE = 60;
static const int WORLD_X_SIZE = 1;
static const int WORLD_Y_SIZE = 1;
static const double DELTA_X = WORLD_X_SIZE / double(GRID_X_SIZE);
static const double DELTA_Y = WORLD_Y_SIZE / double(GRID_Y_SIZE);

int main(int argc, char const *argv[])
{
	// Initialize
	FluidDomain fluid_domain(GRID_X_SIZE, GRID_Y_SIZE, WORLD_X_SIZE, WORLD_Y_SIZE, 0.001);
    FluidSolverMemoryPool mem_pool(fluid_domain);
	FluidSolver fluid_solver(mem_pool);

	// Setup
    fluid_domain.addFluidSource(FluidSource( { 4.0 / GRID_X_SIZE, 0.5, 3.0 / GRID_Y_SIZE, 1 - 4.0 / GRID_Y_SIZE }, DELTA_X, DELTA_Y, 0.0, 0.0, 0.0, 1));
    
	// Initialize glfw
	if (!glfwInit())
		return EXIT_FAILURE;
	// Modern OpenGL
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	// Create a windowed mode window and its OpenGL context
	GLFWwindow* window = glfwCreateWindow(512 * WORLD_X_SIZE, 512 * WORLD_Y_SIZE, "Fluid Simulation", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return EXIT_FAILURE;
	}
	// Make the window's context current
	glfwMakeContextCurrent(window);
	printf("%s\n", glGetString(GL_VERSION));

	FluidRendererGL renderer(512 * WORLD_X_SIZE, 512 * WORLD_Y_SIZE, WORLD_X_SIZE, WORLD_Y_SIZE);

    // Prepare simulation
	MyFloat seconds_per_frame = 0.02;
	
	// Start simulation
	while (!glfwWindowShouldClose(window))
	{
		MyFloat dt;
		// Simulate
		for (MyFloat frame_time = 0; frame_time < seconds_per_frame; frame_time += dt)
		{
			// Calculate dt (for now just set it)
			dt = 0.001;
            dt = CLAMP(dt, 0, seconds_per_frame - frame_time);

			// Update fluid domain (creates new fluid from sources)
			fluid_domain.update(dt);

			// Solve
			try
			{
				//fluid_solver.stepSemiLagrangian(fluid_domain, dt);
				fluid_solver.stepPICFLIP(fluid_domain, dt, 0.05);
			}
			catch (const std::runtime_error& e)
			{
				std::cout << e.what() << std::endl;
				return EXIT_FAILURE;
			}
		}
		// Render
		renderer.renderParticles(fluid_domain.markerParticleSet());
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	return EXIT_SUCCESS;
}
