#include <Simulator.h>

static const int SIZE = 40;

Simulator::Simulator()
{
	_grid = new MacGrid(
		SIZE, 	// Size_x
		SIZE, 	// Size_y
		SIZE, 	// Length_x
		SIZE); 	// Length_y

	for (int j = 0; j < SIZE; ++j)
	{
		for (int i = 0; i < SIZE ; ++i)
		{
			_grid->setColor(i, j, (1 - (i/2)%3) && (1 - (j/2)%3));
		}
	}

	_renderer = new Renderer(_grid);
		
	for (int i = 0; i < 300; ++i)
	{
		std::cout << "iteration " << i << std::endl;

		// Simulate
		double dt = 0.1;

		_grid->addExternalForce(dt);
		_grid->advect(dt);
		_grid->pressureSolve(dt);

		// Render
		std::stringstream str;
		str << "test" << i << ".ppm";

		_renderer->renderCanvas();
		_renderer->writeCanvasToPpm( str.str().c_str());

		//_renderer->writeToPpm(
		//str.str().c_str(),	// File path
		//0,					// Min val
		//1);					// Max val
	}
}

Simulator::~Simulator()
{
	delete _grid;
	delete _renderer;
}