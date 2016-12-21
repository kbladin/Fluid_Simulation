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

	_vel_x_front_buffer(size_x, size_y, _DELTA_X, _DELTA_Y),
	_vel_y_front_buffer(size_x, size_y, _DELTA_X, _DELTA_Y),
	_vel_x_back_buffer(	size_x, size_y, _DELTA_X, _DELTA_Y),
	_vel_y_back_buffer(	size_x, size_y, _DELTA_X, _DELTA_Y),
    _cell_type_buffer(	size_x, size_y)
{
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

// Getters

MyFloat MacGrid::velX(int i, int j) const
{
	return (_vel_x_front_buffer.value(i,j) + _vel_x_front_buffer.value(i+1,j)) / 2;
}

MyFloat MacGrid::velY(int i, int j) const
{
	return (_vel_y_front_buffer.value(i,j) + _vel_y_front_buffer.value(i,j+1)) / 2;
}

MyFloat MacGrid::velXHalfIndexed(int i, int j) const
{
	return _vel_x_front_buffer.value(i,j);
}

MyFloat MacGrid::velYBackBufferHalfIndexed(int i, int j) const
{
	return _vel_y_back_buffer.value(i,j);
}

MyFloat MacGrid::velXBackBufferHalfIndexed(int i, int j) const
{
	return _vel_x_back_buffer.value(i,j);
}

MyFloat MacGrid::velYHalfIndexed(int i, int j) const
{
	return _vel_y_front_buffer.value(i,j);
}

MyFloat MacGrid::velXInterpolated(MyFloat x, MyFloat y) const
{
	MyFloat v_x = _vel_x_front_buffer.valueInterpolated(x, y - _DELTA_Y * 0.5);
	// -0.5 Due to the MAC grid structure
	return v_x;
}

MyFloat MacGrid::velYInterpolated(MyFloat x, MyFloat y) const
{
	MyFloat v_y = _vel_y_front_buffer.valueInterpolated(x - _DELTA_X * 0.5, y);
	// -0.5 Due to the MAC grid structure
	return v_y;
}

CellType MacGrid::cellType(int i, int j) const
{
	// Outside it is always SOLID
	//bool i_outside = i < 0 ? true : (i > _SIZE_X - 1 ? true : false);
	//bool j_outside = j < 0 ? true : (j > _SIZE_Y - 1 ? true : false);
	//bool boundary = i_outside || j_outside;
	return _cell_type_buffer.value(i, j);
}

MyFloat MacGrid::divVelX(int i, int j) const
{
	return 
		(_vel_x_front_buffer.value(i + 1, j)
		- _vel_x_front_buffer.value(i, j))
		/ _DELTA_X;
}

MyFloat MacGrid::divVelY(int i, int j) const
{
	return
		(_vel_y_front_buffer.value(i, j + 1)
		- _vel_y_front_buffer.value(i, j))
		/ _DELTA_Y;
}

MyFloat MacGrid::lengthX() const
{
	return _LENGTH_X;
}

MyFloat MacGrid::lengthY() const
{
	return _LENGTH_Y;
}

int MacGrid::sizeX() const
{
	return _SIZE_X;
}

int MacGrid::sizeY() const
{
	return _SIZE_Y;
}

MyFloat MacGrid::deltaX() const
{
	return _DELTA_X;
}

MyFloat MacGrid::deltaY() const
{
	return _DELTA_Y;
}

glm::dmat2 MacGrid::computeVelocityGradientMatrix(int i, int j)
{
	glm::dmat2 vel_grad;
	vel_grad[0][0] = divVelX(i, j);
	vel_grad[0][1] = (_vel_x_front_buffer.value(i, j + 1)
		- _vel_x_front_buffer.value(i, j))
		/ _DELTA_Y;
	vel_grad[1][0] = (_vel_y_front_buffer.value(i + 1, j)
		- _vel_y_front_buffer.value(i, j))
		/ _DELTA_X;
	vel_grad[1][1] = divVelY(i, j);
	return vel_grad;
}

// Setters

void MacGrid::setVelXHalfIndexed(int i, int j, MyFloat vel_x)
{
    _vel_x_front_buffer(i, j) = vel_x;
}

void MacGrid::setVelYHalfIndexed(int i, int j, MyFloat vel_y)
{
    _vel_y_front_buffer(i, j) = vel_y;
}

void MacGrid::setVelXBackBufferHalfIndexed(int i, int j, MyFloat vel_x)
{
	_vel_x_back_buffer(i, j) = vel_x;
}

void MacGrid::setVelYBackBufferHalfIndexed(int i, int j, MyFloat vel_y)
{
	_vel_y_back_buffer(i, j) = vel_y;
}

void MacGrid::setVelXBackBuffer(int i, int j, MyFloat vel_x)
{
	_vel_x_back_buffer(i, j) = vel_x;
	_vel_x_back_buffer(i + 1, j) = vel_x;
}

void MacGrid::setVelYBackBuffer(int i, int j, MyFloat vel_y)
{
	_vel_y_back_buffer(i, j) = vel_y;
	_vel_y_back_buffer(i, j + 1) = vel_y;
}

void MacGrid::setCellType(int i, int j, CellType cell_type)
{
	_cell_type_buffer(i, j) = cell_type;
}

void MacGrid::addToVelXInterpolated(MyFloat x, MyFloat y, MyFloat vel_x)
{
	_vel_x_back_buffer.addToValueInterpolated(x, y - 0.5 * _DELTA_Y, vel_x);
	// -0.5 Due to the MAC grid structure
}

void MacGrid::addToVelYInterpolated(MyFloat x, MyFloat y, MyFloat vel_y)
{
	_vel_y_back_buffer.addToValueInterpolated(x - 0.5 * _DELTA_X, y, vel_y);
	// -0.5 Due to the MAC grid structure
}

void MacGrid::swapBuffers()
{
	SizedGrid<MyFloat>* vel_x_tmp = &_vel_x_front_buffer;
	SizedGrid<MyFloat>* vel_y_tmp = &_vel_y_front_buffer;

	_vel_x_front_buffer = std::move(_vel_x_back_buffer);
	_vel_y_front_buffer = std::move(_vel_y_back_buffer);

	_vel_x_back_buffer = std::move(*vel_x_tmp);
	_vel_y_back_buffer = std::move(*vel_y_tmp);
}
