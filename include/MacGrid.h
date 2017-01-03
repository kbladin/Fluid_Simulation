#ifndef MAC_GRID_H
#define MAC_GRID_H

#include <Grid.h>

//#include <glm/glm.hpp>

#include <iostream>
#include <assert.h>
#include <memory>

#include "MathDefinitions.h"

enum CellType
{
	LIQUID, AIR, SOLID
};

class MacGrid : public GridInterface
{
public:
	MacGrid(int size_x, int size_y, MyFloat length_x, MyFloat length_y);
	~MacGrid();

    void clearCellTypeBuffer();
    void updatePreviousVelocityBuffer();
    void updateVelocityDiffBuffer();

	// Getters
    inline MyFloat velX(int i, int j) const
    {
        return (_vel_x_front_buffer->value(i,j) +
                _vel_x_front_buffer->value(i+1,j)) / 2;
    };
	inline MyFloat velY(int i, int j) const
    {
        return (_vel_y_front_buffer->value(i,j) +
                _vel_y_front_buffer->value(i,j+1)) / 2;
    };
	inline MyFloat velXHalfIndexed(int i, int j) const
    {
        return _vel_x_front_buffer->value(i,j);
    };
	inline MyFloat velYHalfIndexed(int i, int j) const
    {
        return _vel_y_front_buffer->value(i,j);
    };
	inline MyFloat velXBackBufferHalfIndexed(int i, int j) const
    {
        return _vel_x_back_buffer->value(i,j);
    };
	inline MyFloat velYBackBufferHalfIndexed(int i, int j) const
    {
        return _vel_y_back_buffer->value(i,j);
    };
    inline MyFloat velXBackBuffer(int i, int j) const
    {
        return (_vel_x_back_buffer->value(i,j) +
            _vel_x_back_buffer->value(i + 1, j)) / 2;
    };
    inline MyFloat velYBackBuffer(int i, int j) const
    {
        return (_vel_y_back_buffer->value(i,j) +
            _vel_y_back_buffer->value(i, j + 1)) / 2;
    };
    inline MyFloat velXInterpolated(MyFloat x, MyFloat y) const
    {
        MyFloat v_x =
            _vel_x_front_buffer->valueInterpolated(x, y - _DELTA_Y * 0.5);
        // -0.5 Due to the MAC grid structure
        return v_x;
    };
	inline MyFloat velYInterpolated(MyFloat x, MyFloat y) const
    {
        MyFloat v_y =
            _vel_y_front_buffer->valueInterpolated(x - _DELTA_X * 0.5, y);
        // -0.5 Due to the MAC grid structure
        return v_y;
    };
    inline MyFloat velXDiffInterpolated(MyFloat x, MyFloat y) const
    {
        MyFloat v_x = _vel_x_diff.valueInterpolated(x, y - _DELTA_Y * 0.5);
        // -0.5 Due to the MAC grid structure
        return v_x;
    };
    inline MyFloat velYDiffInterpolated(MyFloat x, MyFloat y) const
    {
        MyFloat v_y = _vel_y_diff.valueInterpolated(x - _DELTA_X * 0.5, y);
        // -0.5 Due to the MAC grid structure
        return v_y;
    };
    inline CellType cellType(int i, int j) const
    {
        i = CLAMP(i, 0, _SIZE_X - 1);
        j = CLAMP(j, 0, _SIZE_Y - 1);
        return _cell_type_buffer.value(i, j);
    };
	inline MyFloat divVelX(int i, int j) const
    {
        return
        (_vel_x_front_buffer->value(i + 1, j)
         - _vel_x_front_buffer->value(i, j))
        / _DELTA_X;
    };
	inline MyFloat divVelY(int i, int j) const
    {
        return
        (_vel_y_front_buffer->value(i, j + 1)
         - _vel_y_front_buffer->value(i, j))
        / _DELTA_Y;
    };
	
	//glm::dmat2 computeVelocityGradientMatrix(int i, int j);
	
	// Setters
	inline void setVelXHalfIndexed(int i, int j, MyFloat vel_x)
    {
        (*_vel_x_front_buffer)(i, j) = vel_x;
    };
    inline void setVelYHalfIndexed(int i, int j, MyFloat vel_y)
    {
        (*_vel_y_front_buffer)(i, j) = vel_y;
    };
    inline void setVelXBackBuffer(int i, int j, MyFloat vel_x)
    {
        (*_vel_x_back_buffer)(i, j) = vel_x;
        (*_vel_x_back_buffer)(i + 1, j) = vel_x;
    };
	inline void setVelYBackBuffer(int i, int j, MyFloat vel_y)
    {
        (*_vel_y_back_buffer)(i, j) = vel_y;
        (*_vel_y_back_buffer)(i, j + 1) = vel_y;
    };
	inline void setVelXBackBufferHalfIndexed(int i, int j, MyFloat vel_x)
    {
        (*_vel_x_back_buffer)(i, j) = vel_x;
    };
	inline void setVelYBackBufferHalfIndexed(int i, int j, MyFloat vel_y)
    {
        (*_vel_y_back_buffer)(i, j) = vel_y;
    };
	inline void setCellType(int i, int j, CellType cell_type)
    {
        _cell_type_buffer(i, j) = cell_type;
    };
	inline void addToVelXInterpolated(MyFloat x, MyFloat y, MyFloat vel_x)
    {
        _vel_x_back_buffer->addToValueInterpolated(
            x, y - 0.5 * _DELTA_Y, vel_x);
        // -0.5 Due to the MAC grid structure
    };
	inline void addToVelYInterpolated(MyFloat x, MyFloat y, MyFloat vel_y)
    {
        _vel_y_back_buffer->addToValueInterpolated(
            x - 0.5 * _DELTA_X, y, vel_y);
        // -0.5 Due to the MAC grid structure
    };
	
	void swapVelocityBuffers();

private:
	// Velocity grids need front buffers and back buffers
	// Normally read from front buffers and write to back buffers, then swap
    std::unique_ptr< Grid<MyFloat> > _vel_x_front_buffer;
	std::unique_ptr< Grid<MyFloat> > _vel_y_front_buffer;

	std::unique_ptr< Grid<MyFloat> > _vel_x_back_buffer;
	std::unique_ptr< Grid<MyFloat> > _vel_y_back_buffer;

    Grid<MyFloat> _vel_x_previous;
    Grid<MyFloat> _vel_y_previous;

    Grid<MyFloat> _vel_x_diff;
    Grid<MyFloat> _vel_y_diff;

	Grid<CellType> _cell_type_buffer;
};

#endif
