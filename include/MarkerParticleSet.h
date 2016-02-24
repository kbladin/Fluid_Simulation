#ifndef PARTICLE_SET_H
#define PARTICLE_SET_H

// Particles stored as linked list
class MarkerParticle
{
public:
	MarkerParticle();
	MarkerParticle(double pos_x, double pos_y);
	~MarkerParticle();
	// Connects p after current particle
	void connect(MarkerParticle* p);

	// Getters
	double posX();
	double posY();

	// Setters
	void setPosition(double pos_x, double pos_y);

	MarkerParticle* next;
private:
	double _pos_x;
	double _pos_y;
};

class MarkerParticleSet
{
public:
	MarkerParticleSet();
	~MarkerParticleSet();

	void addParticle(double pos_x, double pos_y);

	// Returns the first real particle
	MarkerParticle* getFirst() const;
	
private:
	MarkerParticle _dummy_head;
};

#endif