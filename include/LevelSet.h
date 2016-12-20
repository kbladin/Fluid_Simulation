#ifndef LEVELSET_H
#define LEVELSET_H

#include <math.h>

#include <Grid.h>

#include "MathDefinitions.h"

class LevelSet : public SizedGrid<MyFloat>
{
public:
	LevelSet(int size_x, int size_y, MyFloat length_x, MyFloat length_y);
	//LevelSet& operator=(LevelSet to_copy);
	~LevelSet();

	MyFloat distance(int from_i, int from_j, int to_i, int to_j);

	MyFloat computeUpwindGradientX(int i, int j, MyFloat vel_x);
	MyFloat computeUpwindGradientY(int i, int j, MyFloat vel_y);

	MyFloat lengthX() const;
	MyFloat lengthY() const;
private:
	MyFloat _LENGTH_X;
	MyFloat _LENGTH_Y;
};

#endif