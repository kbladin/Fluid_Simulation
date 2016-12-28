#include <FluidDomain.h>

#include <iostream>
#include <random>

FluidSource::FluidSource(
	MyFloat x_min,
	MyFloat x_max,
	MyFloat y_min,
	MyFloat y_max,
	MyFloat time_step,
	int max_creations)
	: _x_min(x_min)
	, _x_max(x_max)
	, _y_min(y_min)
	, _y_max(y_max)
	, _time_step(time_step)
	, _time_since_last(0.0)
	, _max_creations(max_creations)
	, _n_creations(0)
{

}

FluidSource::~FluidSource()
{

}

void FluidSource::update(MarkerParticleSet& particle_set, MyFloat dt)
{
	if (isFinished())
		return;

	if (_time_since_last >= _time_step)
	{
		MyFloat x_incr = 1.0 / 40.0 / 4.0;
		MyFloat y_incr = 1.0 / 40.0 / 4.0;
		for (MyFloat y = _y_min; y < _y_max; y += y_incr)
		{
			for (MyFloat x = _x_min; x < _x_max; x += x_incr)
			{
				particle_set.addParticle(x, y);
			}
		}

		_time_since_last = 0;
		_n_creations++;
	}
	_time_since_last += dt;
}

bool FluidSource::isFinished()
{
	return _max_creations != -1 && _n_creations >= _max_creations;
}

static const int N_PARTICLES = 0;

FluidDomain::FluidDomain(
	int size_x,
	int size_y,
	MyFloat length_x,
	MyFloat length_y,
	MyFloat density)
	: _density(density)
    , _mac_grid(size_x, size_y, length_x, length_y)
    , _level_set(size_x, size_y, length_x, length_y)
    , _particle_set(N_PARTICLES)
{
	// Set positions for all N_PARTICLES particles
	for (auto it = _particle_set.begin(); it != _particle_set.end(); it++)
	{
		it->setPosition(
			(rand() / MyFloat(INT_MAX) / 2 ) * (length_x - length_x / size_x * 2) + length_x / size_x,
			(rand() / MyFloat(INT_MAX) / 2 ) * (length_y - length_y / size_y * 2) + length_y / 2);
	}
}

FluidDomain::~FluidDomain()
{

}

void FluidDomain::addFluidSource(FluidSource fluid_source)
{
	_fluid_sources.push_back(fluid_source);
}

void FluidDomain::update(MyFloat dt)
{
	// Update all fluid sources
	for (int i = 0; i < _fluid_sources.size(); ++i)
	{
		_fluid_sources[i].update(_particle_set, dt);
	}

	// Advect the fluid through the newly updated field
	advectParticles(dt);
	//fluid_domain.advectLevelSet(dt);

	// Classify the cells of the domain (AIR, LIQUID or SOLID)
	classifyCells(_particle_set);
	//fluid_domain.classifyCells(_level_set);
}

void FluidDomain::advectParticles(MyFloat dt)
{
	for (auto it = _particle_set.begin(); it != _particle_set.end(); it++)
	{
		// Position in world
		MyFloat pos_x = it->posX();
		MyFloat pos_y = it->posY();

		// Calculate new position (forward Euler)
		pos_x += _mac_grid.velXInterpolated(pos_x, pos_y) * dt;
		pos_y += _mac_grid.velYInterpolated(pos_x, pos_y) * dt;
		
		// Write data
		it->setPosition(pos_x, pos_y);
	}
}

void FluidDomain::advectLevelSet(MyFloat dt)
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
            MyFloat vel_x = _mac_grid.velX(i, j);
            MyFloat vel_y = _mac_grid.velY(i, j);
			MyFloat grad_x = _level_set.computeUpwindGradientX(i, j, vel_x);
			MyFloat grad_y = _level_set.computeUpwindGradientY(i, j, vel_y);

			MyFloat change_rate = - (vel_x * grad_x + vel_y * grad_y);

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

const MyFloat FluidDomain::density()
{
	return _density;
}
