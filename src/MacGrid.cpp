#include <MacGrid.h>

MacGrid::MacGrid(
	int size_x,
	int size_y,
	MyFloat length_x,
	MyFloat length_y) :

	_SIZE_X(size_x),
	_SIZE_Y(size_y),
	_LENGTH_X(length_x),
	_LENGTH_Y(length_y),
	_DELTA_X(length_x / size_x),
	_DELTA_Y(length_y / size_y),

    _vel_x_previous(size_x, size_y, _DELTA_Y, _DELTA_Y),
    _vel_y_previous(size_x, size_y, _DELTA_Y, _DELTA_Y),
    _vel_x_diff(size_x, size_y, _DELTA_Y, _DELTA_Y),
    _vel_y_diff(size_x, size_y, _DELTA_Y, _DELTA_Y),
    _cell_type_buffer(size_x, size_y)
{
    _vel_x_front_buffer = std::make_unique< SizedGrid<MyFloat> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    _vel_y_front_buffer = std::make_unique< SizedGrid<MyFloat> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    _vel_x_back_buffer = std::make_unique< SizedGrid<MyFloat> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    _vel_y_back_buffer = std::make_unique< SizedGrid<MyFloat> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    
	clearCellTypeBuffer();
}

MacGrid::~MacGrid()
{
    
}

void MacGrid::clearCellTypeBuffer()
{
	// Set all to SOLID
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			_cell_type_buffer(i, j) = SOLID;
		}
	}
	// Set center to AIR
	for (int j = 1; j < _SIZE_Y - 1; ++j)
	{
		for (int i = 1; i < _SIZE_X - 1; ++i)
		{
			_cell_type_buffer(i, j) = AIR;
		}
	}
}

void MacGrid::updatePreviousVelocityBuffer()
{
	_vel_x_previous = *_vel_x_front_buffer;
    _vel_y_previous = *_vel_y_front_buffer;
}

void MacGrid::updateVelocityDiffBuffer()
{
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			_vel_x_diff(i, j) = _vel_x_front_buffer->value(i, j) - _vel_x_previous.value(i, j);
			_vel_y_diff(i, j) = _vel_y_front_buffer->value(i, j) - _vel_y_previous.value(i, j);
		}
	}
}

// Getters

glm::dmat2 MacGrid::computeVelocityGradientMatrix(int i, int j)
{
	glm::dmat2 vel_grad;
	vel_grad[0][0] = divVelX(i, j);
	vel_grad[0][1] = (_vel_x_front_buffer->value(i, j + 1)
		- _vel_x_front_buffer->value(i, j))
		/ _DELTA_Y;
	vel_grad[1][0] = (_vel_y_front_buffer->value(i + 1, j)
		- _vel_y_front_buffer->value(i, j))
		/ _DELTA_X;
	vel_grad[1][1] = divVelY(i, j);
	return vel_grad;
}

void MacGrid::swapBuffers()
{
    _vel_x_front_buffer.swap(_vel_x_back_buffer);
    _vel_y_front_buffer.swap(_vel_y_back_buffer);
}
