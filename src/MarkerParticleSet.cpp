#include <MarkerParticleSet.h>

MarkerParticle::MarkerParticle()
{
	_pos_x = 0;
	_pos_y = 0;
	next = nullptr;
}

MarkerParticle::MarkerParticle(double pos_x, double pos_y)
{
	_pos_x = pos_x;
	_pos_y = pos_y;
	next = nullptr;
}

MarkerParticle::~MarkerParticle()
{
	
}

void MarkerParticle::connect(MarkerParticle* p)
{
	p->next = next;
	next = p;
}

// Getters
double MarkerParticle::posX()
{
	return _pos_x;
}

double MarkerParticle::posY()
{
	return _pos_y;
}

// Setters
void MarkerParticle::setPosition(double pos_x, double pos_y)
{
	_pos_x = pos_x;
	_pos_y = pos_y;
}

MarkerParticleSet::MarkerParticleSet()
{

}

MarkerParticleSet::~MarkerParticleSet()
{
	while (_dummy_head.next)
	{
		MarkerParticle* to_delete = _dummy_head.next;
		_dummy_head.next = to_delete->next;
		delete to_delete;
	}
}

void MarkerParticleSet::addParticle(double pos_x, double pos_y)
{
	MarkerParticle* p = new MarkerParticle(pos_x, pos_y);
	_dummy_head.connect(p);
}

MarkerParticle* MarkerParticleSet::getFirst() const
{
	return _dummy_head.next;
}


