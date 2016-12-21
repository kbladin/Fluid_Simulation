#ifndef MAC_GRID_H
#define MAC_GRID_H

#include <Grid.h>

#include <glm/glm.hpp>

#include <iostream>
#include <assert.h>

#include "MathDefinitions.h"

enum CellType
{
	LIQUID, AIR, SOLID
};

class MacGrid
{
public:
	MacGrid(int size_x, int size_y, MyFloat length_x, MyFloat length_y);
	~MacGrid();

	void clearCellTypeBuffer();

	// Getters
	MyFloat velX(int i, int j) const;
	MyFloat velY(int i, int j) const;
	MyFloat velXHalfIndexed(int i, int j) const;
	MyFloat velYHalfIndexed(int i, int j) const;
	MyFloat velXBackBufferHalfIndexed(int i, int j) const;
	MyFloat velYBackBufferHalfIndexed(int i, int j) const;
	MyFloat velXInterpolated(MyFloat x, MyFloat y) const;
	MyFloat velYInterpolated(MyFloat x, MyFloat y) const;
	MyFloat divVelX(int i, int j) const;
	MyFloat divVelY(int i, int j) const;
	int sizeX() const;
	int sizeY() const;
	MyFloat lengthX() const;
	MyFloat lengthY() const;
	MyFloat deltaX() const;
	MyFloat deltaY() const;
	CellType cellType(int i, int j) const;
	
	glm::dmat2 computeVelocityGradientMatrix(int i, int j);
	
	// Setters
	void setVelXHalfIndexed(int i, int j, MyFloat vel_x);
    void setVelYHalfIndexed(int i, int j, MyFloat vel_y);
    void setVelXBackBuffer(int i, int j, MyFloat vel_x);
	void setVelYBackBuffer(int i, int j, MyFloat vel_y);
	void setVelXBackBufferHalfIndexed(int i, int j, MyFloat vel_x);
	void setVelYBackBufferHalfIndexed(int i, int j, MyFloat vel_y);
	void setCellType(int i, int j, CellType cell_type);
	void addToVelXInterpolated(MyFloat x, MyFloat y, MyFloat vel_x);
	void addToVelYInterpolated(MyFloat x, MyFloat y, MyFloat vel_y);
	
	void swapBuffers();

private:
	// Constants
	const int _SIZE_X;
	const int _SIZE_Y;
	const MyFloat _LENGTH_X;
	const MyFloat _LENGTH_Y;
	const MyFloat _DELTA_X;
	const MyFloat _DELTA_Y;

	// Velocity grids need front buffers and back buffers
	// Normally read from front buffers and write to back buffers, then swap.
	SizedGrid<MyFloat> _vel_x_front_buffer;
	SizedGrid<MyFloat> _vel_y_front_buffer;

	SizedGrid<MyFloat> _vel_x_back_buffer;
	SizedGrid<MyFloat> _vel_y_back_buffer;

	Grid<CellType> _cell_type_buffer;
};

#endif
