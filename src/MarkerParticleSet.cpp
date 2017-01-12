#include <MarkerParticleSet.h>

MarkerParticle::MarkerParticle() :
	_pos_x(0),
	_pos_y(0),
	_vel_x(0),
	_vel_y(0)
{

}

MarkerParticle::MarkerParticle(
	MyFloat pos_x,
	MyFloat pos_y,
	MyFloat vel_x,
	MyFloat vel_y) :
	_pos_x(pos_x),
	_pos_y(pos_y),
	_vel_x(vel_x),
	_vel_y(vel_y)
{

}

MarkerParticle::~MarkerParticle()
{
	
}

MarkerParticleSet::MarkerParticleSet(int size)
{
	_particles.resize(size);
}

MarkerParticleSet::~MarkerParticleSet()
{

}

void MarkerParticleSet::addParticle(MarkerParticle p)
{
	_particles.push_back(p);
}

void MarkerParticleSet::reserve(int particle_count)
{
	_particles.reserve(particle_count);
}

void MarkerParticleSet::advect(MyFloat dt)
{
	for (auto it = _particles.begin(); it != _particles.end(); it++)
	{
		it->advect(dt);
	}
}

void MarkerParticleSet::advectAndEnsureOutsideObstacles(
	MyFloat dt,
	const MacGrid& mac_grid)
{
	for (auto it = _particles.begin(); it != _particles.end(); it++)
	{
		it->advect(dt);
		int x = it->posX() / mac_grid.deltaX();
		int y = it->posY() / mac_grid.deltaY();
		if (mac_grid.cellType(x, y) == SOLID)
		{
			it->advect(-dt);
		}
	}
}
