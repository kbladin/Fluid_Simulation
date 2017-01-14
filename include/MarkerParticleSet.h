#ifndef PARTICLE_SET_H
#define PARTICLE_SET_H

#include <vector>

#include "MacGrid.h"
#include "MathDefinitions.h"

class MarkerParticle
{
public:
	MarkerParticle();
	MarkerParticle(
		MyFloat pos_x,
		MyFloat pos_y,
		MyFloat vel_x = 0,
		MyFloat vel_y = 0);
	~MarkerParticle();

	// Getters
	inline MyFloat posX() const { return _pos_x; };
	inline MyFloat posY() const { return _pos_y; };
	inline MyFloat velX() const { return _vel_x; };
	inline MyFloat velY() const { return _vel_y; };

	// Setters
	inline void setPosition(MyFloat pos_x, MyFloat pos_y)
	{
		_pos_x = pos_x;
		_pos_y = pos_y;
	};
	inline void setVelocity(MyFloat vel_x, MyFloat vel_y)
	{
		_vel_x = vel_x;
		_vel_y = vel_y;
	};
	inline void advect(MyFloat dt)
	{
		_pos_x += _vel_x * dt;
		_pos_y += _vel_y * dt;
	};
	
private:
	// Position
	MyFloat _pos_x;
	MyFloat _pos_y;

	// Velocity
	MyFloat _vel_x;
	MyFloat _vel_y;
};

class MarkerParticleSet
{
public:
	MarkerParticleSet(int size = 0);
	~MarkerParticleSet();

	void addParticle(MarkerParticle p) 	{ _particles.push_back(p); };
	void reserve(int particle_count)	{ _particles.reserve(particle_count); };
	void clear() 						{ _particles.clear(); };
	void advect(MyFloat dt);
	void advectAndEnsureOutsideObstacles(MyFloat dt, const MacGrid& mac_grid);

	typedef std::vector<MarkerParticle>::iterator iterator;
	typedef std::vector<MarkerParticle>::const_iterator const_iterator;

	inline int size() const { return _particles.size(); };
	
	inline MarkerParticleSet::iterator begin() { return _particles.begin(); };
	inline MarkerParticleSet::iterator end() { return _particles.end(); };

	inline MarkerParticleSet::const_iterator begin() const { return _particles.begin(); };
	inline MarkerParticleSet::const_iterator end() const { return _particles.end(); };

private:
	std::vector<MarkerParticle> _particles;
};

#endif
