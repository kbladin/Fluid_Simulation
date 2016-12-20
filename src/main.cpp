#include <FluidSolver.h>
#include <FluidDomain.h>
#include <Renderer.h>

#include <iostream>
#include <sstream>
#include <string>

static const int GRID_SIZE = 40;
static const int WORLD_SIZE = 1;

int main(int argc, char const *argv[])
{
	FluidDomain fluid_domain(GRID_SIZE, GRID_SIZE, WORLD_SIZE, WORLD_SIZE);
	FluidSolver fluid_solver;
	Renderer renderer(0,0, WORLD_SIZE, WORLD_SIZE);

	// Setup
	int n_frames = 300;
	MyFloat seconds_per_frame = 0.02;

	// Classify the cells of the domain (AIR, LIQUID or SOLID)
	fluid_domain.classifyCells(fluid_domain.markerParticleSet());
	//fluid_domain.classifyCells(fluid_domain.levelSet());

	// Start simulation
	for (int i = 0; i < n_frames; ++i)
	{
		std::cout << "iteration " << i << std::endl;

		MyFloat dt;
		// Simulate
		for (MyFloat frame_time = 0; frame_time < seconds_per_frame; frame_time += dt)
		{
			// Calculate dt (for now just set it)
			dt = 0.0005;

			if (dt > seconds_per_frame - frame_time)
			{
				dt = seconds_per_frame - frame_time;
			}

			// Solve
			fluid_solver.step(fluid_domain, dt);

			// Advect the fluid through the newly updated field
			fluid_domain.advectParticles(dt);
			//fluid_domain.advectLevelSet(dt);

			// Classify the cells of the domain (AIR, LIQUID or SOLID)
			fluid_domain.classifyCells(fluid_domain.markerParticleSet());
			//fluid_domain.classifyCells(fluid_domain.levelSet());
		}
		// Render
		renderer.clearCanvas();
		renderer.renderGridCellsToCanvas(fluid_domain.macGrid());
		//renderer.renderGridVelocitiesToCanvas(fluid_domain.macGrid());
		renderer.renderParticlesToCanvas(fluid_domain.markerParticleSet());
		
		std::stringstream file_name;
		file_name << "test" << i << ".ppm";
		renderer.writeCanvasToPpm(file_name.str().c_str());
	}

	return 0;
}
