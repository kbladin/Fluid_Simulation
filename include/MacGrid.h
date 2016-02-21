#ifndef MAC_GRID_H
#define MAC_GRID_H

#include <vector>
#include <iostream>
#include <random>

#include <assert.h>

#include <armadillo>

//! Let this have dense arrays for now
class MacGrid
{
public:
	MacGrid(int size_x, int size_y, double length_x, double length_y);
	~MacGrid();

	// Simple semi Lagrangian advection with Euler integration
	void advect(double dt);
	void addExternalForce(double dt);
	void pressureSolve(double dt);

	// Getters
	double velX(int x, int y);
	double velY(int x, int y);
	double velXHalfIndexed(int x, int y);
	double velYHalfIndexed(int x, int y);
	double color(int x, int y);
	double velXInterpolated(double x, double y);
	double velYInterpolated(double x, double y);
	double colorInterpolated(double x, double y);
	double divVelX(int x, int y);
	double divVelY(int x, int y);
	int sizeX();
	int sizeY();
	double lengthX();
	double lengthY();
	double deltaX();
	double deltaY();

	bool isSolid(int i, int j);

	// Setters
	void setVelX(int x, int y, double vel_x);
	void setVelY(int x, int y, double vel_y);
	void setColor(int x, int y, double color);
	void addToVelXInterpolated(double x, double y, double vel_x);
	void addToVelYInterpolated(double x, double y, double vel_y);
	void addToColorInterpolated(double x, double y, double color);

	void _swapBuffers();

private:
	void _advectVelX(double dt);
	void _advectVelY(double dt);
	void _advectColor(double dt);

	// Constants
	const int _SIZE_X;
	const int _SIZE_Y;
	const double _LENGTH_X;
	const double _LENGTH_Y;
	const double _DELTA_X;
	const double _DELTA_Y;

	// Always render to back bufer from front buffer, then swap them
	// Since advection can not be done in place, another set of data is needed
	// (except for when adding forces)
	double* _vel_x_front_buffer;
	double* _vel_y_front_buffer;
	double* _color_front_buffer;

	double* _vel_x_back_buffer;
	double* _vel_y_back_buffer;
	double* _color_back_buffer;
};

#endif