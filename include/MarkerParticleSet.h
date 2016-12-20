#ifndef PARTICLE_SET_H
#define PARTICLE_SET_H

#include <vector>

#include "MathDefinitions.h"

class MarkerParticle
{
public:
	MarkerParticle();
	MarkerParticle(MyFloat pos_x, MyFloat pos_y);
	~MarkerParticle();

	// Getters
	MyFloat posX() const;
	MyFloat posY() const;

	// Setters
	void setPosition(MyFloat pos_x, MyFloat pos_y);
private:
	MyFloat _pos_x;
	MyFloat _pos_y;
};

class MarkerParticleSet
{
public:
	MarkerParticleSet(int size);
	~MarkerParticleSet();

	void addParticle(MyFloat pos_x, MyFloat pos_y);

	typedef std::vector<MarkerParticle>::iterator iterator;
	typedef std::vector<MarkerParticle>::const_iterator const_iterator;

	MarkerParticleSet::iterator begin();
	MarkerParticleSet::iterator end();

	MarkerParticleSet::const_iterator begin() const;
	MarkerParticleSet::const_iterator end() const;

private:
	std::vector<MarkerParticle> _particles;
};

#endif