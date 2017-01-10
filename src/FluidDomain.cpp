#include <FluidDomain.h>

FluidSource::FluidSource(
	BBox<MyFloat> area,
	MyFloat delta_x,
	MyFloat delta_y,
	MyFloat x_velocity,
	MyFloat y_velocity,
	MyFloat time_step,
	int max_spawns) :
	_area(area),
	_x_velocity(x_velocity),
	_y_velocity(y_velocity),
	_delta_x(delta_x),
	_delta_y(delta_y),
	_time_step(time_step),
	_time_since_last(0.0),
	_max_spawns(max_spawns),
	_n_spawns(0)
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
		// Roughly 6.25 particles per cell
		MyFloat x_incr = _delta_x / 2.5;
		MyFloat y_incr = _delta_y / 2.5;
		for (MyFloat y = _area.y_min; y < _area.y_max; y += y_incr)
		{
			for (MyFloat x = _area.x_min; x < _area.x_max; x += x_incr)
			{
				particle_set.addParticle(
					MarkerParticle(x, y, _x_velocity, _y_velocity));
			}
		}

		_time_since_last = 0;
		_n_spawns++;
	}
	_time_since_last += dt;
}

bool FluidSource::isFinished()
{
	return _max_spawns != -1 && _n_spawns >= _max_spawns;
}

FluidDomain::FluidDomain(
	int size_x,
	int size_y,
	MyFloat length_x,
	MyFloat length_y,
	MyFloat density) :
	GridInterface(size_x, size_y, length_x / size_x, length_y / size_y),
	_density(density),
    _mac_grid(size_x, size_y, length_x, length_y),
    _level_set(size_x, size_y, length_x, length_y)
{

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
        
        x = CLAMP(x, 0, _mac_grid.sizeX()-1);
        y = CLAMP(y, 0, _mac_grid.sizeY()-1);
        
		_mac_grid.setCellType(x, y, LIQUID);
	}
	// Reset border
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
