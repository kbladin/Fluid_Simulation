#ifndef LEVELSET_H
#define LEVELSET_H

#include <math.h>

#include <Grid.h>

class LevelSet : public SizedGrid<double>
{
public:
	LevelSet(int size_x, int size_y, double length_x, double length_y);
	//LevelSet& operator=(LevelSet to_copy);
	~LevelSet();

	double distance(int from_i, int from_j, int to_i, int to_j);

	double computeUpwindGradientX(int i, int j, double vel_x);
	double computeUpwindGradientY(int i, int j, double vel_y);

	double lengthX() const;
	double lengthY() const;
private:
	double _LENGTH_X;
	double _LENGTH_Y;
};

#endif