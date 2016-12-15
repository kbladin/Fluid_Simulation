#ifndef MAC_GRID_H
#define MAC_GRID_H

#include <Grid.h>

#include <glm/glm.hpp>

#include <iostream>
#include <assert.h>

enum CellType
{
	LIQUID, AIR, SOLID
};

template class Grid<double>;
template class Grid<CellType>;

class MacGrid
{
public:
	MacGrid(int size_x, int size_y, double length_x, double length_y);
	~MacGrid();

	void clearCellTypeBuffer();

	// Index transforms
	void worldToCell(double x, double y, int* i, int* j) const;
	void cellToWorld(int i, int j, double* x, double* y) const;
	void linearTo2DCellCenter(int idx, int* i, int* j) const;
	int twoDToLinearCellCenter(int i, int j) const;

	// Getters
	double velX(int i, int j) const;
	double velY(int i, int j) const;
	double velXHalfIndexed(int i, int j) const;
	double velYHalfIndexed(int i, int j) const;
	double velXBackBufferHalfIndexed(int i, int j) const;
	double velYBackBufferHalfIndexed(int i, int j) const;
	double velXInterpolated(double x, double y) const;
	double velYInterpolated(double x, double y) const;
	double divVelX(int i, int j) const;
	double divVelY(int i, int j) const;
	glm::dmat2 computeVelocityGradientMatrix(int i, int j);
	int sizeX() const;
	int sizeY() const;
	double lengthX() const;
	double lengthY() const;
	double deltaX() const;
	double deltaY() const;
	CellType cellType(int i, int j) const;
	
	// Setters
	void setVelXHalfIndexed(int i, int j, double vel_x);
    void setVelYHalfIndexed(int i, int j, double vel_y);
    void setVelXBackBuffer(int i, int j, double vel_x);
	void setVelYBackBuffer(int i, int j, double vel_y);
	void setVelXBackBufferHalfIndexed(int i, int j, double vel_x);
	void setVelYBackBufferHalfIndexed(int i, int j, double vel_y);
	void setCellType(int i, int j, CellType cell_type);
	void addToVelXInterpolated(double x, double y, double vel_x);
	void addToVelYInterpolated(double x, double y, double vel_y);
	
	void swapBuffers();

private:
	// Constants
	const int _SIZE_X;
	const int _SIZE_Y;
	const double _LENGTH_X;
	const double _LENGTH_Y;
	const double _DELTA_X;
	const double _DELTA_Y;

	// Always render to back buffer from front buffer, then swap them
	// since advection can not be done in place, another set of data is needed
	// (except for when adding forces)
	SizedGrid<double> _vel_x_front_buffer;
	SizedGrid<double> _vel_y_front_buffer;

	SizedGrid<double> _vel_x_back_buffer;
	SizedGrid<double> _vel_y_back_buffer;

	Grid<CellType> _cell_type_buffer;
};

#endif
