#ifndef LEVELSET_H
#define LEVELSET_H

#include <math.h>

#include <Grid.h>

#include "MathDefinitions.h"

class LevelSet : public Grid<MyFloat>
{
public:
	LevelSet(int size_x, int size_y, MyFloat length_x, MyFloat length_y);
	~LevelSet();

	MyFloat distance(int from_i, int from_j, int to_i, int to_j);

	MyFloat computeUpwindGradientX(int i, int j, MyFloat vel_x);
	MyFloat computeUpwindGradientY(int i, int j, MyFloat vel_y);
};

#endif
