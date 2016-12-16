#include <MacGrid.h>

MacGrid::MacGrid(
	int size_x,
	int size_y,
	double length_x,
	double length_y) :

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

// Index transforms
void MacGrid::worldToCell(double x, double y, int* i, int* j) const
{
	*i = x / _DELTA_X;
	*j = y / _DELTA_Y;
}

void MacGrid::cellToWorld(int i, int j, double* x, double* y) const
{
	*x = i * _DELTA_X;
	*y = j * _DELTA_Y;
}

void MacGrid::linearTo2DCellCenter(int idx, int* i, int* j) const
{
	*i = idx % _SIZE_X;
	*j = idx / _SIZE_X;
}

int MacGrid::twoDToLinearCellCenter(int i, int j) const
{
	return i + j * _SIZE_X;
}

// Getters
/**
	Returns an interpolated value due to the mac grid structure.
*/
double MacGrid::velX(int i, int j) const
{
	return (_vel_x_front_buffer.value(i,j) + _vel_x_front_buffer.value(i+1,j)) / 2;
}
/**
	Returns an interpolated value due to the mac grid structure.
*/
double MacGrid::velY(int i, int j) const
{
	return (_vel_y_front_buffer.value(i,j) + _vel_y_front_buffer.value(i,j+1)) / 2;
}

/**
	Returns u[x - 1/2][y], velocities are stored in the borders of the cells.
*/
double MacGrid::velXHalfIndexed(int i, int j) const
{
	return _vel_x_front_buffer.value(i,j);
}

/**
	Returns u[x][y - 1/2], velocities are stored in the borders of the cells.
*/
double MacGrid::velYHalfIndexed(int i, int j) const
{
	return _vel_y_front_buffer.value(i,j);
}

/**
	Input is in world coordinates, not necessary on cell borders.
	Currently simple linear interpolation.
*/
double MacGrid::velXInterpolated(double x, double y) const
{
	double v_x = _vel_x_front_buffer.valueInterpolated(x, y - _DELTA_Y * 0.5);
	// -0.5 Due to the MAC grid structure
	return v_x;
}

/**
	Input is in world coordinates, not necessary on cell borders.
	Currently simple linear interpolation.
*/
double MacGrid::velYInterpolated(double x, double y) const
{
	double v_y = _vel_y_front_buffer.valueInterpolated(x - _DELTA_X * 0.5, y);
	// -0.5 Due to the MAC grid structure
	return v_y;
}

CellType MacGrid::cellType(int i, int j) const
{
	// Outside it is always SOLID
	bool i_outside = i < 0 ? true : (i > _SIZE_X - 1 ? true : false);
	bool j_outside = j < 0 ? true : (j > _SIZE_Y - 1 ? true : false);
	bool boundary = i_outside || j_outside;
	if (boundary)
	{
		return SOLID;
	}
	else
	{
		return _cell_type_buffer.value(i, j);
	}
}

double MacGrid::divVelX(int i, int j) const
{
	return 
		(_vel_x_front_buffer.value(i + 1, j)
		- _vel_x_front_buffer.value(i, j))
		/ _DELTA_X;
}

double MacGrid::divVelY(int i, int j) const
{
	return
		(_vel_y_front_buffer.value(i, j + 1)
		- _vel_y_front_buffer.value(i, j))
		/ _DELTA_Y;
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

double MacGrid::lengthX() const
{
	return _LENGTH_X;
}

double MacGrid::lengthY() const
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

double MacGrid::deltaX() const
{
	return _DELTA_X;
}

double MacGrid::deltaY() const
{
	return _DELTA_Y;
}

// Setters
/**
*/

void MacGrid::setVelXHalfIndexed(int i, int j, double vel_x)
{
    _vel_x_front_buffer(i, j) = vel_x;
}

void MacGrid::setVelYHalfIndexed(int i, int j, double vel_y)
{
    _vel_y_front_buffer(i, j) = vel_y;
}

void MacGrid::setVelXBackBufferHalfIndexed(int i, int j, double vel_x)
{
	_vel_x_back_buffer(i, j) = vel_x;
}

/**
*/
void MacGrid::setVelYBackBufferHalfIndexed(int i, int j, double vel_y)
{
	_vel_y_back_buffer(i, j) = vel_y;
}

void MacGrid::setVelXBackBuffer(int i, int j, double vel_x)
{
	_vel_x_back_buffer(i, j) = vel_x;
	_vel_x_back_buffer(i + 1, j) = vel_x;
}

/**
*/
void MacGrid::setVelYBackBuffer(int i, int j, double vel_y)
{
	_vel_y_back_buffer(i, j) = vel_y;
	_vel_y_back_buffer(i, j + 1) = vel_y;
}

void MacGrid::setCellType(int i, int j, CellType cell_type)
{
	_cell_type_buffer(i, j) = cell_type;
}

/**
	Writes to the four closest grid points, writes to back buffer.
*/
void MacGrid::addToVelXInterpolated(double x, double y, double vel_x)
{
	_vel_x_back_buffer.addToValueInterpolated(x, y - 0.5 * _DELTA_Y, vel_x);
	// -0.5 Due to the MAC grid structure
}

/**
	Writes to the four closest grid points, writes to back buffer.
*/
void MacGrid::addToVelYInterpolated(double x, double y, double vel_y)
{
	_vel_y_back_buffer.addToValueInterpolated(x - 0.5 * _DELTA_X, y, vel_y);
	// -0.5 Due to the MAC grid structure
}

void MacGrid::swapBuffers()
{
	SizedGrid<double> vel_x_tmp = std::move(_vel_x_front_buffer);
	SizedGrid<double> vel_y_tmp = std::move(_vel_y_front_buffer);

	_vel_x_front_buffer = std::move(_vel_x_back_buffer);
	_vel_y_front_buffer = std::move(_vel_y_back_buffer);

	_vel_x_back_buffer = std::move(vel_x_tmp);
	_vel_y_back_buffer = std::move(vel_y_tmp);
}
