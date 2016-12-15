#include <Simulator.h>

static const int GRID_SIZE = 40;
static const int WORLD_SIZE = 40;
static const int N_PARTICLES = 10000;

Simulator::Simulator()
{
	_grid = std::unique_ptr<MacGrid>(new MacGrid(
		GRID_SIZE, 		// Size_x
		GRID_SIZE, 		// Size_y
		WORLD_SIZE, 	// Length_x
		WORLD_SIZE)); 	// Length_y
	_level_set = std::unique_ptr<LevelSet>(
		new LevelSet( GRID_SIZE, GRID_SIZE, WORLD_SIZE, WORLD_SIZE));
	_particle_set = std::unique_ptr<MarkerParticleSet>(
		new MarkerParticleSet(N_PARTICLES));
	_renderer = std::unique_ptr<Renderer>(
		new Renderer(0,0, WORLD_SIZE, WORLD_SIZE));
	
	/*
	for (int j = 0; j < _grid->sizeY(); ++j)
	{
		for (int i = 0; i < _grid->sizeX(); ++i)
		{
			_grid->setColor(i, j, j%2);
		}
	}
	*/

	for (int j = 0; j < _level_set->sizeY(); ++j)
	{
		for (int i = 0; i < _level_set->sizeX(); ++i)
		{
			(*_level_set)(i, j) = 10;
		}
	}
	for (int j = _level_set->sizeY() / 4; j < _level_set->sizeY() * 3 / 4.0; ++j)
	{
		for (int i = _level_set->sizeX() / 4; i < _level_set->sizeX() * 3 / 4.0; ++i)
		{
			(*_level_set)(i, j) = -10;
		}
	}

	// Set positions for all N_PARTICLES particles
	for (auto it = _particle_set->begin(); it != _particle_set->end(); it++)
	{
		it->setPosition(
			(rand() / double(INT_MAX) / 2 + 0.25) * WORLD_SIZE,
			(rand() / double(INT_MAX) / 2 + 0.25) * WORLD_SIZE);
	}

	// Setup
	int n_frames = 100;
	double seconds_per_frame = 0.2;
	
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
			updateCellTypesWithParticles();
			//updateCellTypesWithLevelSet();
			if(i < 15)
			_grid->addExternalForce(dt, 0, 10);
			_grid->advect(dt);
			_grid->enforceDirichlet();
			_grid->pressureSolve(dt);
			_grid->enforceDirichlet();
			ExtendVelocity();

			// Advect particles through fluid
			advectParticles(dt);
			//advectLevelSet(dt);
		}

		// Render
		_renderer->clearCanvas();
		//_renderer->renderColorToCanvas(_grid.get());
		//_renderer->renderGridCellsToCanvas(_grid.get());
		//_renderer->renderLevelSetFunctionValuesToCanvas(_level_set.get());
		//_renderer->renderGridVelocitiesToCanvas(_grid.get());
		_renderer->renderParticlesToCanvas(_particle_set.get());
		
		std::stringstream str;
		str << "test" << i << ".ppm";		
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

void Simulator::advectLevelSet(double dt)
{
	LevelSet new_level_set(
		_level_set->sizeX(),
		_level_set->sizeY(),
		_level_set->lengthX(),
		_level_set->lengthY());
	for (int j = 0; j < _level_set->sizeX(); ++j)
	{
		for (int i = 0; i < _level_set->sizeY(); ++i)
		{
			double vel_x = _grid->velX(i, j);
			double vel_y = _grid->velY(i, j);
			double grad_x = _level_set->computeUpwindGradientX(i, j, vel_x);
			double grad_y = _level_set->computeUpwindGradientY(i, j, vel_y);

			double change_rate = - (vel_x * grad_x + vel_y * grad_y);

			// Forward Euler
			new_level_set(i,j) = (*_level_set)(i,j) + change_rate * dt * 5;
		}
	}
	*_level_set = std::move(new_level_set);
	
}

/**
	Extends the velocity of the mac grid from the interface of the level set.
*/
/*
void Simulator::ExtendVelocity()
{
	double dt = 0.7;
	int iterations = 5;
	for (int iter = 0; iter < iterations; ++iter)
	{
		for (int j = 0; j < _level_set->sizeX(); ++j)
		{
			for (int i = 0; i < _level_set->sizeY(); ++i)
			{
				// Currently doing extension everywhere. Not just narrow band.
				// This is bad for efficiency
                // Extend only in air or solid
                if (_grid->cellType(i, j) == LIQUID)
                {
                    _grid->setVelXBackBuffer(i, j, _grid->velX(i, j));
                    _grid->setVelYBackBuffer(i, j, _grid->velY(i, j));
                }
                else {
					double vel_x = _grid->velX(i, j);
					double vel_y = _grid->velY(i, j);
					double grad_x = _level_set->computeUpwindGradientX(i, j, vel_x);
					double grad_y = _level_set->computeUpwindGradientY(i, j, vel_y);

					double level_set_val = (*_level_set)(i, j);

					glm::dvec2 normal(grad_x, grad_y);
					normal = glm::normalize(normal);

					glm::dmat2 vel_grad = _grid->computeVelocityGradientMatrix(i, j);

					glm::dvec2 change_rate =
						- (vel_grad * normal) * glm::sign(level_set_val);

					// Forward Euler
					glm::dvec2 new_vel = glm::dvec2(vel_x, vel_y) + change_rate * dt;
					_grid->setVelX(i, j, new_vel.x);
					_grid->setVelY(i, j, new_vel.y);
				}
			}
		}
		_grid->swapBuffers();
	}
}
*/


/**
	Extends the velocity of the mac grid from the interface of the level set.
*/
/*
void Simulator::ExtendVelocity()
{
    for (int j = 0; j < _level_set->sizeX(); ++j)
    {
        for (int i = 0; i < _level_set->sizeY(); ++i)
        {     
            _grid->setVelXBackBuffer(i, j, 0);
            _grid->setVelYBackBuffer(i, j, 0);
        }
    }

    double dt = 0.7;
    int iterations = 10;
    for (int iter = 0; iter < iterations; ++iter)
    {
        for (int j = 0; j < _level_set->sizeX(); ++j)
        {
            for (int i = 0; i < _level_set->sizeY(); ++i)
            {
                // Extend only in air or solid
                if (_grid->cellType(i, j) == LIQUID)
                {
                    _grid->setVelXBackBuffer(i, j, _grid->velX(i, j));
                    _grid->setVelYBackBuffer(i, j, _grid->velY(i, j));
                }
                else
                {
                    glm::dvec2 normal((*_level_set)(i+1, j) - (*_level_set)(i-1, j),
                                      (*_level_set)(i, j+1) - (*_level_set)(i, j-1));
                    
                    if (glm::length(normal) > 0)
                        normal = glm::normalize(normal);
                    
                    //double vel_x = _grid->velX(i, j);
                    //double vel_y = _grid->velY(i, j);
                    
                    glm::dvec2 v  = glm::dvec2(_grid->velX(i, j), _grid->velY(i, j));

                    glm::dvec2 ip = glm::vec2(_grid->velXHalfIndexed(i+1, j), _grid->velY(i+1, j));
                    glm::dvec2 im = glm::vec2(_grid->velXHalfIndexed(i, j), _grid->velY(i-1, j));
                    glm::dvec2 jp = glm::vec2(_grid->velX(i, j+1), _grid->velYHalfIndexed(i, j+1));
                    glm::dvec2 jm = glm::vec2(_grid->velX(i, j-1), _grid->velYHalfIndexed(i, j));
                    
                    glm::dvec2 v0 =  v - dt * (glm::max(0.0, normal.x)*(v - im) + glm::min(0.0, normal.x)*(ip - v) +
                                               glm::max(0.0, normal.y)*(v - jm) + glm::min(0.0, normal.y)*(jp - v));
                  
                    _grid->setVelXBackBuffer(i, j, v0.x);
                    _grid->setVelYBackBuffer(i, j, v0.y);
                }
            }
        }
        _grid->swapBuffers();
    }   
}
*/

void Simulator::ExtendVelocity()
{
	Grid<int> valid_mask(_grid->sizeX(), _grid->sizeY());
    for (int j = 0; j < _grid->sizeX(); ++j)
    {
        for (int i = 0; i < _grid->sizeX(); ++i)
        {
        	if(_grid->cellType(i, j) == LIQUID)
        	{
        		valid_mask(i,j) = 1;
				_grid->setVelXBackBufferHalfIndexed(i, j, _grid->velXHalfIndexed(i, j));
	            _grid->setVelYBackBufferHalfIndexed(i, j, _grid->velYHalfIndexed(i, j));
        	}
        }
    }

    int iterations = 2;
    for (int iter = 0; iter < iterations; ++iter)
    {
        for (int j = 0; j < _grid->sizeX(); ++j)
        {
            for (int i = 0; i < _grid->sizeY(); ++i)
            {
            	if (valid_mask.value(i, j) == 0)
            	{
	            	double new_vel_x = 0;
	            	double new_vel_y = 0;
	            	int n_valid_neighbors = 0;

	            	// Get values of all neighbors
	            	if(valid_mask.value(i-1, j) == 1)
	            	{
	            		new_vel_x += _grid->velX(i-1, j);
	            		new_vel_y += _grid->velY(i-1, j);
	            		n_valid_neighbors++;
	            	}
	            	if(valid_mask.value(i+1, j) == 1)
	            	{
	            		new_vel_x += _grid->velX(i+1, j);
	            		new_vel_y += _grid->velY(i+1, j);
	            		n_valid_neighbors++;
	            	}
	            	if(valid_mask.value(i, j-1) == 1)
	            	{
	            		new_vel_x += _grid->velX(i, j-1);
	            		new_vel_y += _grid->velY(i, j-1);
	            		n_valid_neighbors++;
	            	}
	            	if(valid_mask.value(i, j+1) == 1)
	            	{
	            		new_vel_x += _grid->velX(i, j+1);
	            		new_vel_y += _grid->velY(i, j+1);
	            		n_valid_neighbors++;
	            	}

	            	// Average the value for the current cell
	            	if (n_valid_neighbors > 0)
	            	{
	            		new_vel_x /= n_valid_neighbors;
	            		new_vel_y /= n_valid_neighbors;

	            		_grid->setVelXBackBufferHalfIndexed(i, j, new_vel_x);
						_grid->setVelYBackBufferHalfIndexed(i, j, new_vel_y);
	            		valid_mask(i, j) = 1;
	            	}
            	}
            }
        }
    }
    _grid->swapBuffers();
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
    for (int j = 0; j < _grid->sizeY(); ++j)
    {
        for (int i = 0; i < _grid->sizeX(); ++i)
        {
            if (i == 0 || j == 0 || i == _grid->sizeX() - 1 || j == _grid->sizeY() - 1)
            {
                _grid->setCellType(i, j, SOLID);
            }
        }
    }
}

void Simulator::updateCellTypesWithLevelSet()
{
	// First reset types (set all cells to AIR)
	_grid->clearCellTypeBuffer();

	for (int j = 0; j < _grid->sizeY(); ++j)
	{
		for (int i = 0; i < _grid->sizeX(); ++i)
		{
			if ((*_level_set)(i,j) < 0)
				_grid->setCellType(i, j, LIQUID);
            
            if (i == 0 || j == 0 || i == _grid->sizeX() - 1 || j == _grid->sizeY() - 1)
            {
                _grid->setCellType(i, j, SOLID);
            }
		}
	}
}

Simulator::~Simulator()
{

}
