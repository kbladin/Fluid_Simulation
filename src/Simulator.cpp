#include <Simulator.h>

static const int GRID_SIZE = 40;
static const int WORLD_SIZE = 10;
static const int N_PARTICLES = 500;

Simulator::Simulator()
{
	_grid = std::unique_ptr<MacGrid>(new MacGrid(
		GRID_SIZE, 		// Size_x
		GRID_SIZE, 		// Size_y
		WORLD_SIZE, 	// Length_x
		WORLD_SIZE)); 	// Length_y
	_particle_set = std::unique_ptr<MarkerParticleSet>(
		new MarkerParticleSet(N_PARTICLES));
	_renderer = std::unique_ptr<Renderer>(
		new Renderer(0,0, WORLD_SIZE, WORLD_SIZE));
	
	for (int j = 0; j < _grid->sizeY(); ++j)
	{
		for (int i = 0; i < _grid->sizeX(); ++i)
		{
			_grid->setColor(i, j, j%2);
		}
	}

	// Set positions for all N_PARTICLES particles
	for (auto it = _particle_set->begin(); it != _particle_set->end(); it++)
	{
		it->setPosition(
			rand() / double(INT_MAX) * WORLD_SIZE,
			rand() / double(INT_MAX) * WORLD_SIZE);
	}

	// Setup
	int n_frames = 100;
	double seconds_per_frame = 1;
	
	// Start simulation
	for (int i = 0; i < n_frames; ++i)
	{
		std::cout << "iteration " << i << std::endl;

		double dt;
		// Simulate
		for (double frame_time = 0; frame_time < seconds_per_frame; frame_time += dt)
		{
			// Calculate dt (for now just set it)
			dt = 0.01;
			// Update the fluid grid
			//updateCellTypesWithParticles();
			if (i < 10)
			_grid->addExternalForce(dt, 0, 10);
			_grid->advect(dt);
			_grid->enforceDirichlet();
			_grid->pressureSolve(dt);
			_grid->enforceDirichlet();

			// Advect particles through fluid
			advectParticles(dt);
		}

		// Render
		std::stringstream str;
		str << "test" << i << ".ppm";


		_renderer->clearCanvas();
		//_renderer->renderColorToCanvas(_grid);
		//_renderer->renderGridCellsToCanvas(_grid);
		_renderer->renderParticlesToCanvas(_particle_set.get());
		_renderer->renderGridVelocitiesToCanvas(_grid.get());
		
		_renderer->writeCanvasToPpm(str.str().c_str());
	}
}

void Simulator::advectParticles(double dt)
{
	for (auto it = _particle_set->begin(); it != _particle_set->end(); it++)
	{
		// Position in world
		double pos_x = it->posX();
		double pos_y = it->posY();

		// Calculate new position
		pos_x += _grid->velXInterpolated(pos_x, pos_y) * dt;
		pos_y += _grid->velYInterpolated(pos_x, pos_y) * dt;
		
		// Write data
		it->setPosition(pos_x, pos_y);
	}
}

void Simulator::updateCellTypesWithParticles()
{
	// First reset types (set all cells to AIR)
	_grid->clearCellTypeBuffer();

	// Loop through all particles and set cells to liquid if they
	// contain a particle
	for (auto it = _particle_set->begin(); it != _particle_set->end(); it++)
	{
		// Find the particles position in the grid
		int x = (it->posX() / _grid->lengthX()) * _grid->sizeX();
		int y = (it->posY() / _grid->lengthY()) * _grid->sizeY();

		_grid->setCellType(x, y, LIQUID);
	}
}

Simulator::~Simulator()
{

}