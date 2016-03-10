#ifndef PARTICLE_SET_H
#define PARTICLE_SET_H

#include <vector>

class MarkerParticle
{
public:
	MarkerParticle();
	MarkerParticle(double pos_x, double pos_y);
	~MarkerParticle();

	// Getters
	double posX() const;
	double posY() const;

	// Setters
	void setPosition(double pos_x, double pos_y);
private:
	double _pos_x;
	double _pos_y;
};

class MarkerParticleSet
{
public:
	MarkerParticleSet(int size);
	~MarkerParticleSet();

	void addParticle(double pos_x, double pos_y);

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