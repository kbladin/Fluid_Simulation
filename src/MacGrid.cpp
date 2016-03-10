#include <MacGrid.h>

#define CLAMP(x, low, high) x < low ? low : (x > high ? high : x)

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

	_vel_x_front_buffer(size_x + 1, size_y, _DELTA_X, _DELTA_Y),
	_vel_y_front_buffer(size_x, size_y + 1, _DELTA_X, _DELTA_Y),
	_color_front_buffer(size_x, size_y, 	_DELTA_X, _DELTA_Y),
	_vel_x_back_buffer(size_x + 1, size_y, 	_DELTA_X, _DELTA_Y),
	_vel_y_back_buffer(size_x, size_y + 1, 	_DELTA_X, _DELTA_Y),
	_color_back_buffer(size_x, size_y, 		_DELTA_X, _DELTA_Y),
	_cell_type_buffer(size_x, size_y)
{
	// Allocate matrix
	A = Eigen::SparseMatrix<double>(_SIZE_X * _SIZE_Y, _SIZE_X * _SIZE_Y);

	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			int idx = 			twoDToLinearCellCenter(i,j);
			int idx_i_minus1 = 	twoDToLinearCellCenter(i - 1,j);
			int idx_i_plus1 = 	twoDToLinearCellCenter(i + 1,j);
			int idx_j_minus1 = 	twoDToLinearCellCenter(i,j - 1);
			int idx_j_plus1 = 	twoDToLinearCellCenter(i,j + 1);

			// Set values in A. Check in all dimensions
			int n_non_solid_neighbors = 0;
			if (cellType(i - 1, j) != SOLID) {
				A.insert(idx_i_minus1, idx) = 1;
				n_non_solid_neighbors++;
			}	
			if (cellType(i + 1, j) != SOLID) {
				A.insert(idx_i_plus1, idx) = 1;
				n_non_solid_neighbors++;
			}
			if (cellType(i, j - 1) != SOLID) {
				A.insert(idx_j_minus1, idx) = 1;
				n_non_solid_neighbors++;
			}
			if (cellType(i, j + 1) != SOLID) {
				A.insert(idx_j_plus1, idx) = 1;
				n_non_solid_neighbors++;
			}

			// Set diagonal value
			A.insert(idx, idx) = - n_non_solid_neighbors;
		}
	}
}

MacGrid::~MacGrid()
{

}

void MacGrid::advect(double dt)
{
	// Advection is done separately since the attributes vel_x, vel_y and
	// color are strored in different grids (MAC grid)
	// Advect color through fluid
	_advectColor(dt);
	// Self advection of velocity components
	_advectVelX(dt);
	_advectVelY(dt);

	_swapBuffers();
}

void MacGrid::addExternalForce(double dt, double F_x, double F_y)
{
	/*
	for (int j = _SIZE_Y * 3 / 7; j < _SIZE_Y * 4 / 7; ++j)
	{
		for (int i = _SIZE_X * 3 / 7; i < _SIZE_X * 4 / 7; ++i)
		{
			if (cellTypeXHalfIndexed(i, j) == LIQUID)
			{ // Only add force to the liquid cells
				// Euler integration (here write directly to front buffer for now)
				_vel_x_front_buffer(i,j) = 
					_vel_x_front_buffer(i,j) + F_x * dt;
				_vel_y_front_buffer(i,j) = 
					_vel_y_front_buffer(i,j) + F_y * dt;
			}
		}
	}
*/
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i =0 ; i < _SIZE_X; ++i)
		{
			if (cellType(i, j) == LIQUID)
			{ // Only add force to the liquid cells
				// Euler integration (here write directly to front buffer for now)
				double x_pos = (i + 0.5) * _DELTA_X;
				double y_pos = (j + 0.5) * _DELTA_Y;
				double vel_y = velYInterpolated(x_pos, y_pos);
				addToVelYInterpolated(x_pos, y_pos, F_y * dt);
			}
		}
	}
	_swapBuffers();
}

void MacGrid::pressureSolve(double dt)
{
    int n_elements = _SIZE_X * _SIZE_Y;
    // Sparse matrix containing connectivity information
    // for the laplace operator
    //Eigen::SparseMatrix<double> A(n_elements, n_elements);

    // Vector containing negative divergences in each cell
	Eigen::VectorXd b(n_elements);    

	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			int idx = twoDToLinearCellCenter(i,j);

			// Calculate divergence and store in b
			b[idx] = - (divVelX(i, j) + divVelY(i, j));
		}
	}

	// Solver of linear system
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
	solver.setMaxIterations(5);

	// Vector containing pressures for each cell
	Eigen::VectorXd x(n_elements);
	solver.compute(A);
	x = solver.solve(b);

	double density = 0.1;
	// X vel
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 1; i < _SIZE_X; ++i)
		{
			// Calculate indices
			int idx_center = 			twoDToLinearCellCenter(i, j);
			int idx_center_i_minus1 = 	twoDToLinearCellCenter(i - 1, j);
			// Get the current velocity and pressure difference
			double vel_x = _vel_x_front_buffer(i,j);
			double pressure_diff = x(idx_center) - x(idx_center_i_minus1);
			// Calculate new velocity
			double new_vel_x = vel_x + dt * 1 / density * pressure_diff / _DELTA_X;
			// Write data
			_vel_x_back_buffer(i,j) = new_vel_x;
		}
	}
	// Y vel
	for (int j = 1; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			// Calculate indices
			int idx_center = 			twoDToLinearCellCenter(i, j);
			int idx_center_j_minus1 = 	twoDToLinearCellCenter(i, j - 1);
			// Get the current velocity and pressure difference
			double vel_y = _vel_y_front_buffer(i,j);
			double pressure_diff = x(idx_center) - x(idx_center_j_minus1);
			// Calculate new velocity
			double new_vel_y = vel_y + dt * 1 / density * pressure_diff / _DELTA_Y;
			// Write data
			_vel_y_back_buffer(i,j) = new_vel_y;
		}
	}

	// Color (does not change)
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			double new_color = _color_front_buffer(i,j);
			_color_back_buffer(i,j) = new_color;
		}
	}

	_swapBuffers();
}

void MacGrid::enforceDirichlet()
{
	// X vel
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X + 1; ++i)
		{
			if(	(cellType(i - 1, j) == SOLID && velXHalfIndexed(i,j) < 0) ||
				(cellType(i, j)  == SOLID && velXHalfIndexed(i,j) > 0))
				setVelX(i, j, 0);
		}
	}
	// Y vel
	for (int j = 0; j < _SIZE_Y + 1; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			if(	(cellType(i, j - 1) == SOLID && velYHalfIndexed(i,j) < 0) ||
				(cellType(i, j) == SOLID && velYHalfIndexed(i,j) > 0))
				setVelY(i, j, 0);
		}
	}
}

void MacGrid::clearCellTypeBuffer()
{
	// Set all to AIR
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
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

double MacGrid::color(int i, int j) const
{
	return _color_front_buffer.value(i,j);
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

double MacGrid::colorInterpolated(double x, double y) const
{
	double color = _vel_y_front_buffer.valueInterpolated(
		x - _DELTA_X * 0.5,
		y - _DELTA_Y * 0.5);
	// -0.5 Due to the MAC grid structure
	return color;
}

CellType MacGrid::cellType(int i, int j) const
{
	// Currently only boundaries are solid
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

/**
	Sampling on borders yields LIQUID if the point neighbors a liquid cell
*/
CellType MacGrid::cellTypeXHalfIndexed(int i, int j) const
{
	// Look at the cells on either side of the point
	// (which is on the border of two cells)
	CellType left_type = cellType(i, j);
	CellType right_type = cellType(i - 1, j);

	return
		left_type == SOLID || right_type == SOLID ? SOLID :
		(left_type == LIQUID || right_type == LIQUID ? LIQUID :
		AIR);
}

/**
	Sampling on borders yields LIQUID if the point neighbors a liquid cell
*/
CellType MacGrid::cellTypeYHalfIndexed(int i, int j) const
{
	// Look at the cells on either side of the point
	// (which is on the border of two cells)
	CellType upper_type = cellType(i, j);
	CellType lower_type = cellType(i, j - 1);

	return
		upper_type == SOLID || lower_type == SOLID ? SOLID :
		(upper_type == LIQUID || lower_type == LIQUID ? LIQUID :
		AIR);
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
void MacGrid::setVelX(int i, int j, double vel_x)
{
	_vel_x_front_buffer(i, j) = vel_x;
}

/**
*/
void MacGrid::setVelY(int i, int j, double vel_y)
{
	_vel_y_front_buffer(i, j) = vel_y;
}

void MacGrid::setColor(int i, int j, double color)
{
	_color_front_buffer(i, j) = color;
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

/**
	Writes to the four closest grid points, writes to back buffer.
*/
void MacGrid::addToColorInterpolated(double x, double y, double color)
{
	_color_back_buffer.addToValueInterpolated(
		x - 0.5 * _DELTA_X,
		y - 0.5 * _DELTA_Y,
		color);
	// -0.5 Due to the MAC grid structure
}

void MacGrid::_getAdvectedPositionRK3(
		double x_pos,
		double y_pos,
		double dt,
		double* x,
		double* y)
{
	double k_1x = velXInterpolated(x_pos, y_pos);
	double k_1y = velYInterpolated(x_pos, y_pos);
	
	double k_2x = velXInterpolated(
		x_pos + 1.0/2 * dt * k_1x,
		y_pos + 1.0/2 * dt * k_1y);
	double k_2y = velYInterpolated(
		x_pos + 1.0/2 * dt * k_1x,
		y_pos + 1.0/2 * dt * k_1y);

	double k_3x = velXInterpolated(
		x_pos + 3.0/4 * dt * k_2x,
		y_pos + 3.0/4 * dt * k_2y);
	double k_3y = velYInterpolated(
		x_pos + 3.0/4 * dt * k_2x,
		y_pos + 3.0/4 * dt * k_2y);

	double x_pos_prev = x_pos -
	dt * 1.0/9 * (
		2 * k_1x +
		3 * k_2x +
		4 * k_3x);
	double y_pos_prev = y_pos -
	dt * 1.0/6 * (
		2 * k_1y +
		3 * k_2y +
		4 * k_3y);

	*x = x_pos_prev;
	*y = y_pos_prev;
}

void MacGrid::_getAdvectedPositionForwardEuler(
		double x_pos,
		double y_pos,
		double dt,
		double* x,
		double* y)
{
	// Velcoity
	double v_x = velXInterpolated(x_pos, y_pos);
	double v_y = velYInterpolated(x_pos, y_pos);

	// Previous particle position
	double x_pos_prev = x_pos - v_x * dt;
	double y_pos_prev = y_pos - v_y * dt;

	*x = x_pos_prev;
	*y = y_pos_prev;
}

void MacGrid::_advectVelX(double dt)
{
	// Set all to 0
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X + 1; ++i)
		{
			// Write data
			_vel_x_back_buffer(i, j) = 0;
		}
	}
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X + 1; ++i)
		{
			if (cellTypeXHalfIndexed(i, j) == LIQUID)
			{ // Only advect the liquid cells
				double x_pos = i * _DELTA_X;
				double y_pos = (j + 0.5) * _DELTA_Y;
				double x_pos_prev, y_pos_prev;
//#define USE_EULER
#ifdef USE_EULER
				_getAdvectedPositionForwardEuler(
					x_pos,
					y_pos,
					dt,
					&x_pos_prev,
					&y_pos_prev);
#else
				_getAdvectedPositionRK3(
					x_pos,
					y_pos,
					dt,
					&x_pos_prev,
					&y_pos_prev);
#endif
				// Velocity to be stored
				double v_x = velXInterpolated(x_pos, y_pos);
				addToVelXInterpolated(x_pos_prev, y_pos_prev, v_x);
			}
		}
	}
}

void MacGrid::_advectVelY(double dt)
{
	// Set all to 0
	for (int j = 0; j < _SIZE_Y + 1; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			// Write data
			_vel_y_back_buffer(i, j) = 0;
		}
	}
	for (int j = 0; j < _SIZE_Y + 1; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			if (cellTypeYHalfIndexed(i, j) == LIQUID)
			{ // Only advect the liquid cells
				// Particle position
				double x_pos = (i + 0.5) * _DELTA_X;
				double y_pos = j * _DELTA_Y;
				double x_pos_prev, y_pos_prev;
#ifdef USE_EULER
				_getAdvectedPositionForwardEuler(
					x_pos,
					y_pos,
					dt,
					&x_pos_prev,
					&y_pos_prev);
#else
				_getAdvectedPositionRK3(
					x_pos,
					y_pos,
					dt,
					&x_pos_prev,
					&y_pos_prev);
#endif
				// Velocity to be stored
				double v_y = velYInterpolated(x_pos, y_pos);
				addToVelYInterpolated(x_pos_prev, y_pos_prev, v_y);
			}
		}
	}
}

void MacGrid::_advectColor(double dt)
{
	// Set all to 0
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			// Write data
			_color_back_buffer(i, j) = 0;
		}
	}
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			if (cellType(i, j) == LIQUID)
			{ // Only advect the liquid cells
				// Particle position
				double x_pos = (i + 0.5) * _DELTA_X;
				double y_pos = (j + 0.5) * _DELTA_Y;
				double x_pos_prev, y_pos_prev;
#ifdef USE_EULER
				_getAdvectedPositionForwardEuler(
					x_pos,
					y_pos,
					dt,
					&x_pos_prev,
					&y_pos_prev);
#else
				_getAdvectedPositionRK3(
					x_pos,
					y_pos,
					dt,
					&x_pos_prev,
					&y_pos_prev);
#endif
				// Color to be stored
				double color = colorInterpolated(x_pos, y_pos);
				addToColorInterpolated(x_pos_prev, y_pos_prev, color);
			}
		}
	}
}

void MacGrid::_swapBuffers()
{
	SizedGrid<double> vel_x_tmp = std::move(_vel_x_front_buffer);
	SizedGrid<double> vel_y_tmp = std::move(_vel_y_front_buffer);
	SizedGrid<double> color_tmp = std::move(_color_front_buffer);

	_vel_x_front_buffer = std::move(_vel_x_back_buffer);
	_vel_y_front_buffer = std::move(_vel_y_back_buffer);
	_color_front_buffer = std::move(_color_back_buffer);

	_vel_x_back_buffer = std::move(vel_x_tmp);
	_vel_y_back_buffer = std::move(vel_y_tmp);
	_color_back_buffer = std::move(color_tmp);
}
