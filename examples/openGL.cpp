#include <FluidSolver.h>
#include <FluidDomain.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>

#include <gl/glew.h>
#include <gl/glfw3.h>

#include <SGE/SimpleGraphicsEngine.h>

// This class extends SimpleGraphicsEngine,
// Before initializing this object, an OpenGL context must be created
class MyEngine : public SimpleGraphicsEngine
{
public:
  MyEngine();
  ~MyEngine();

  void update();
};

MyEngine::MyEngine() : SimpleGraphicsEngine()
{
  
}

MyEngine::~MyEngine()
{
  
}

void MyEngine::update()
{
  SimpleGraphicsEngine::render();
}

static const int GRID_X_SIZE = 64;
static const int GRID_Y_SIZE = 64;
static const int WORLD_X_SIZE = 1;
static const int WORLD_Y_SIZE = 1;
static const double DELTA_X = WORLD_X_SIZE / double(GRID_X_SIZE);
static const double DELTA_Y = WORLD_Y_SIZE / double(GRID_Y_SIZE);

int main(int argc, char const *argv[])
{
	// Initialize
	FluidDomain fluid_domain(GRID_X_SIZE, GRID_Y_SIZE, WORLD_X_SIZE, WORLD_Y_SIZE, 0.01);
    FluidSolverMemoryPool mem_pool(fluid_domain);
	FluidSolver fluid_solver(mem_pool);

	// Setup
    fluid_domain.addFluidSource(FluidSource( { 2.0 / GRID_X_SIZE, 0.35, 2.0 / GRID_Y_SIZE, 1 - 2.0 / GRID_Y_SIZE }, DELTA_X, DELTA_Y, 0.0, 0.0, 0.0, 1));
    //fluid_domain.addFluidSource(FluidSource( { 1.0 / GRID_X_SIZE, 0.2 }, DELTA_X, DELTA_Y, 0.1, 0.4, 1.5, 5));
    //fluid_domain.addFluidSource(FluidSource( { 0.8, 1.0 - 1.0 / GRID_Y_SIZE, 0.1, 0.4 }, DELTA_X, DELTA_Y, 0.0, 0.0, 1.5, 5));
    //fluid_domain.addFluidSource(FluidSource( { 0.8, 1.2, 0.1, 0.5 }, DELTA_X, DELTA_Y, 0.0, 0.0, 2.0, 2));
    
    //fluid_domain.addFluidSource(FluidSource( { 1.0 / GRID_X_SIZE, 2.0 / GRID_X_SIZE, 0.3, 0.35 }, DELTA_X, DELTA_Y, 3.0, 0.0, 0.003, 300));
    //fluid_domain.addFluidSource(FluidSource( { 1.0 - 2.0 / GRID_X_SIZE, 1.0 - 1.0 / GRID_X_SIZE, 0.3, 0.35 }, DELTA_X, DELTA_Y, -3.0, 0.0, 0.003, 300));
    //fluid_domain.addFluidSource(FluidSource( { 0.17, 0.23, 1.0 / GRID_Y_SIZE, 2.0 / GRID_Y_SIZE }, DELTA_X, DELTA_Y, 0.0, 3.0, 0.003, 300));
    //fluid_domain.addFluidSource(FluidSource( { 1-0.23, 1-0.17, 1.0 / GRID_Y_SIZE, 2.0 / GRID_Y_SIZE }, DELTA_X, DELTA_Y, 0.0, 3.0, 0.003, 300));
    
	// Initialize glfw
	if (!glfwInit())
		return EXIT_FAILURE;
	// Modern OpenGL
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	// Create a windowed mode window and its OpenGL context
	GLFWwindow* window = glfwCreateWindow(512, 512, "Fluid Simulation", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return EXIT_FAILURE;
	}
	// Make the window's context current
	glfwMakeContextCurrent(window);
	printf("%s\n", glGetString(GL_VERSION));

	MyEngine e;

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
			dt = 0.01;
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

		e.update();
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	return EXIT_SUCCESS;
}
