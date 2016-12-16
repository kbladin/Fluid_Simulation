#include <FluidDomain.h>

#include <iostream>
#include <random>

static const int N_PARTICLES = 20000;

FluidDomain::FluidDomain(int size_x, int size_y, double length_x, double length_y)
	: _mac_grid(size_x, size_y, length_x, length_y)
	, _level_set(size_x, size_y, length_x, length_y)
	, _particle_set(N_PARTICLES)
{
	// Set positions for all N_PARTICLES particles
	for (auto it = _particle_set.begin(); it != _particle_set.end(); it++)
	{
		//it->setPosition(0.5 * size_x, 0.5 * size_y);
		it->setPosition(
			(rand() / double(INT_MAX) / 2 ) * length_x,
			(rand() / double(INT_MAX) ) * length_y);
	}
}

FluidDomain::~FluidDomain()
{

}

void FluidDomain::addFluidSource()
{
	std::cout << "FluidDomain::addFluidSource() not implemented!" << std::endl;
}

void FluidDomain::addExternalForce(double F_x, double F_y, double dt)
{
	/*
	for (int j = _SIZE_Y * 3 / 7; j < _SIZE_Y * 4 / 7; ++j)
	{
		for (int i = _SIZE_X * 3 / 7; i < _SIZE_X * 4 / 7; ++i)
		{
			if (cellTypeXHalfIndexed(i, j) == LIQUID)
			{ // Only add force to the liquid cells
				// Euler integration (here write directly to front buffer for now)
				_vel_x_front_buffer(i,j) = 
					_vel_x_front_buffer(i,j) + F_x * dt;
				_vel_y_front_buffer(i,j) = 
					_vel_y_front_buffer(i,j) + F_y * dt;
			}
		}
	}
*/
	for (int j = 0; j < _mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < _mac_grid.sizeY(); ++i)
		{
			if (_mac_grid.cellType(i, j) == LIQUID)
			{ // Only add force to the liquid cells
                _mac_grid.setVelXHalfIndexed(i,j, _mac_grid.velXHalfIndexed(i,j) + F_x * dt);
                _mac_grid.setVelYHalfIndexed(i,j, _mac_grid.velYHalfIndexed(i,j) + F_y * dt);
			}
		}
	}
	//_swapBuffers();
}

void FluidDomain::advectParticles(double dt)
{
	for (auto it = _particle_set.begin(); it != _particle_set.end(); it++)
	{
		// Position in world
		double pos_x = it->posX();
		double pos_y = it->posY();

		// Calculate new position (forward Euler)
		pos_x += _mac_grid.velXInterpolated(pos_x, pos_y) * dt;
		pos_y += _mac_grid.velYInterpolated(pos_x, pos_y) * dt;
		
		// Write data
		it->setPosition(pos_x, pos_y);
	}
}

void FluidDomain::advectLevelSet(double dt)
{
	LevelSet new_level_set(
		_level_set.sizeX(),
		_level_set.sizeY(),
		_level_set.lengthX(),
		_level_set.lengthY());
	for (int j = 0; j < _level_set.sizeX(); ++j)
	{
		for (int i = 0; i < _level_set.sizeY(); ++i)
		{
            double vel_x = _mac_grid.velX(i, j);
            double vel_y = _mac_grid.velY(i, j);
			double grad_x = _level_set.computeUpwindGradientX(i, j, vel_x);
			double grad_y = _level_set.computeUpwindGradientY(i, j, vel_y);

			double change_rate = - (vel_x * grad_x + vel_y * grad_y);

			// Forward Euler
			new_level_set(i,j) = _level_set(i,j) + change_rate * dt * 5;
		}
	}
	_level_set = std::move(new_level_set);
}

void FluidDomain::classifyCells(LevelSet& levelSet)
{
	// First reset types (set all cells to AIR)
	_mac_grid.clearCellTypeBuffer();

	for (int j = 0; j < _mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < _mac_grid.sizeX(); ++i)
		{
			if (_level_set(i,j) < 0)
				_mac_grid.setCellType(i, j, LIQUID);
            
            if (i == 0 || j == 0 || i == _mac_grid.sizeX() - 1 ||
            	j == _mac_grid.sizeY() - 1)
            {
                _mac_grid.setCellType(i, j, SOLID);
            }
		}
	}
}

void FluidDomain::classifyCells(MarkerParticleSet& particle_set)
{
	// First reset types (set all cells to AIR)
	_mac_grid.clearCellTypeBuffer();

	// Loop through all particles and set cells to liquid if they
	// contain a particle
	for (auto it = _particle_set.begin(); it != _particle_set.end(); it++)
	{
		// Find the particles position in the grid
		int x = (it->posX() / _mac_grid.lengthX()) * _mac_grid.sizeX();
		int y = (it->posY() / _mac_grid.lengthY()) * _mac_grid.sizeY();

		_mac_grid.setCellType(x, y, LIQUID);
	}
	// Reset border (ugly hack)
    for (int j = 0; j < _mac_grid.sizeY(); ++j)
    {
        for (int i = 0; i < _mac_grid.sizeX(); ++i)
        {
            if (i == 0 || j == 0 || i == _mac_grid.sizeX() - 1 ||
            	j == _mac_grid.sizeY() - 1)
            {
                _mac_grid.setCellType(i, j, SOLID);
            }
        }
    }
}

MacGrid& FluidDomain::macGrid()
{
	return _mac_grid;
}

LevelSet& FluidDomain::levelSet()
{
	return _level_set;
}

MarkerParticleSet& FluidDomain::markerParticleSet()
{
	return _particle_set;
}
