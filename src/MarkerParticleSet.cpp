#include <MarkerParticleSet.h>

MarkerParticle::MarkerParticle()
{
	_pos_x = 0;
	_pos_y = 0;
}

MarkerParticle::MarkerParticle(double pos_x, double pos_y)
{
	_pos_x = pos_x;
	_pos_y = pos_y;
}

MarkerParticle::~MarkerParticle()
{
	
}

// Getters
double MarkerParticle::posX() const
{
	return _pos_x;
}

double MarkerParticle::posY() const
{
	return _pos_y;
}

// Setters
void MarkerParticle::setPosition(double pos_x, double pos_y)
{
	_pos_x = pos_x;
	_pos_y = pos_y;
}

MarkerParticleSet::MarkerParticleSet(int size)
{
	_particles.resize(size);
}

MarkerParticleSet::~MarkerParticleSet()
{

}

void MarkerParticleSet::addParticle(double pos_x, double pos_y)
{
	_particles.push_back(MarkerParticle(pos_x, pos_y));
}

MarkerParticleSet::iterator MarkerParticleSet::begin()
{
	return _particles.begin();
}

MarkerParticleSet::iterator MarkerParticleSet::end()
{
	return _particles.end();
}

MarkerParticleSet::const_iterator MarkerParticleSet::begin() const
{
	return _particles.begin();
}

MarkerParticleSet::const_iterator MarkerParticleSet::end() const
{
	return _particles.end();
}
