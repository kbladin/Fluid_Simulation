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
	_DELTA_Y(length_y / size_y)
{
	_vel_x_front_buffer = new double[(size_x + 1) * size_y];
	_vel_y_front_buffer = new double[size_x * (size_y + 1)];
	_color_front_buffer = new double[size_x * size_y];

	_vel_x_back_buffer = new double[(size_x + 1) * size_y];
	_vel_y_back_buffer = new double[size_x * (size_y + 1)];
	_color_back_buffer = new double[size_x * size_y];

	_cell_type_buffer = new CellType[size_x * size_y];

	// Reset data
	for (int i = 0; i < size_x * size_y; ++i)
	{
		_color_front_buffer[i] = 0;
		_color_back_buffer[i] = 0;
		_cell_type_buffer[i] = LIQUID;
	}
	for (int i = 0; i < (size_x + 1) * size_y; ++i)
	{
		_vel_x_front_buffer[i] = 0;
		_vel_x_back_buffer[i] = 0;
	}
	for (int i = 0; i < size_x * (size_y + 1); ++i)
	{
		_vel_y_front_buffer[i] = 0;
		_vel_y_back_buffer[i] = 0;
	}

	A = Eigen::SparseMatrix<double>(_SIZE_X * _SIZE_Y, _SIZE_X * _SIZE_Y);

	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			int idx = i + j * _SIZE_X;
			int idx_i_minus1 = (i - 1) + j * _SIZE_X;
			int idx_i_plus1 = (i + 1) + j * _SIZE_X;
			int idx_j_minus1 = i + (j - 1) * _SIZE_X;
			int idx_j_plus1 = i + (j + 1) * _SIZE_X;

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
	delete[] _vel_x_front_buffer;
	delete[] _vel_y_front_buffer;
	delete[] _color_front_buffer;
	delete[] _cell_type_buffer;

	delete[] _vel_x_back_buffer;
	delete[] _vel_y_back_buffer;
	delete[] _color_back_buffer;
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
	for (int j = _SIZE_Y * 1 / 7; j < _SIZE_Y * 2 / 7; ++j)
	{
		for (int i = _SIZE_X * 1 / 7; i < _SIZE_X * 2 / 7; ++i)
		{
			if (cellTypeXHalfIndexed(i, j) == LIQUID)
			{ // Only add force to the liquid cells
				// Euler integration (here write directly to front buffer for now)
				_vel_x_front_buffer[i + j * (_SIZE_X + 1)] = 
					_vel_x_front_buffer[i + j * (_SIZE_X + 1)] + F_x * dt;
				_vel_y_front_buffer[i + j * _SIZE_X] = 
					_vel_y_front_buffer[i + j * _SIZE_X] + F_y * dt;
			}
		}
	}

	//for (int j = 0; j < _SIZE_Y; ++j)
	//{
	//	for (int i =0 ; i < _SIZE_X; ++i)
	//	{
	//		if (cellType(i, j) == LIQUID)
	//		{ // Only add force to the liquid cells
	//			// Euler integration (here write directly to front buffer for now)
	//			double x_pos = (i + 0.5) * _DELTA_X;
	//			double y_pos = (j + 0.5) * _DELTA_Y;
	//			double vel_y = velYInterpolated(x_pos, y_pos);
	//			addToVelYInterpolated(x_pos, y_pos, F_y * dt);
	//		}
	//	}
	//}
	//_swapBuffers();
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
			int idx = i + j * _SIZE_X;

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
			double vel_x = _vel_x_front_buffer[i + j * (_SIZE_X + 1)];
			double pressure_diff =
				x(i + j * _SIZE_X)
				- x(i - 1 + j * _SIZE_X);
			double new_vel_x = vel_x + dt * 1 / density * pressure_diff / _DELTA_X;
			_vel_x_back_buffer[i + j * (_SIZE_X + 1)] = new_vel_x;
		}
	}
	// Y vel
	for (int j = 1; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			double vel_y = _vel_y_front_buffer[i + j * _SIZE_X];
			double pressure_diff =
				x(i + j * _SIZE_X)
				- x(i + (j - 1) * _SIZE_X);
			double new_vel_y = vel_y + dt * 1 / density * pressure_diff / _DELTA_Y;
			_vel_y_back_buffer[i + j * _SIZE_X] = new_vel_y;
		}
	}

	// Color (does not change)
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			double new_color = _color_front_buffer[i + j * _SIZE_X];
			_color_back_buffer[i + j * _SIZE_X] = new_color;
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
			_cell_type_buffer[i + j * _SIZE_X] = AIR;
		}
	}
}

// Getters
/**
	Returns an interpolated value due to the mac grid structure.
*/
double MacGrid::velX(int x, int y) const
{
	return (_vel_x_front_buffer[x + y * (_SIZE_X + 1)] +
		_vel_x_front_buffer[x + 1 + y * (_SIZE_X + 1)]) / 2;
}

/**
	Returns an interpolated value due to the mac grid structure.
*/
double MacGrid::velY(int x, int y) const
{
	return (_vel_y_front_buffer[x + y * _SIZE_X] +
		_vel_y_front_buffer[x + (y + 1) * _SIZE_X]) / 2;
}

/**
	Returns u[x - 1/2][y], velocities are stored in the borders of the cells.
*/
double MacGrid::velXHalfIndexed(int x, int y) const
{
	return _vel_x_front_buffer[x + y * (_SIZE_X + 1)];
}

/**
	Returns u[x][y - 1/2], velocities are stored in the borders of the cells.
*/
double MacGrid::velYHalfIndexed(int x, int y) const
{
	return _vel_y_front_buffer[x + y * _SIZE_X];
}

double MacGrid::color(int x, int y) const
{
	return _color_front_buffer[x + y * _SIZE_X];
}

/**
	Input is in world coordinates, not necessary on cell borders.
	Currently simple linear interpolation.
*/
double MacGrid::velXInterpolated(double x, double y) const
{
	// Calculate indices
	int i = x / _DELTA_X;
	int j = y / _DELTA_Y - 0.5; // -0.5 Due to the MAC grid structure
	int i_plus1 = i + 1;
	int j_plus1 = j + 1;
	double i_frac = x / _DELTA_X - i;
	double j_frac = y / _DELTA_Y - 0.5 - j;

	assert(i_frac <= 1);
	assert(j_frac <= 1);

	// Border cases
	i = CLAMP(i, 0, _SIZE_X + 1 - 1);
	j = CLAMP(j, 0, _SIZE_Y - 1);
	i_plus1 = CLAMP(i_plus1, 0, _SIZE_X + 1 - 1);
	j_plus1 = CLAMP(j_plus1, 0, _SIZE_Y - 1);
	
	// First interpolate in x, then in y
	double v_x_00 = _vel_x_front_buffer[i + j * (_SIZE_X + 1)];
	double v_x_10 = _vel_x_front_buffer[i_plus1 + j * (_SIZE_X + 1)];
	double v_x_01 = _vel_x_front_buffer[i + j_plus1 * (_SIZE_X + 1)];
	double v_x_11 = _vel_x_front_buffer[i_plus1 + j_plus1 * (_SIZE_X + 1)];
	
	// Interpolate x
	double v_x_0 = (1 - i_frac) * v_x_00 + i_frac * v_x_10;
	double v_x_1 = (1 - i_frac) * v_x_01 + i_frac * v_x_11; 

	// Interpolate y
	double v_x = (1 - j_frac) * v_x_0 + j_frac * v_x_1;
	return v_x;
}

/**
	Input is in world coordinates, not necessary on cell borders.
	Currently simple linear interpolation.
*/
double MacGrid::velYInterpolated(double x, double y) const
{
	// Calculate indices
	int i = x / _DELTA_X - 0.5; // -0.5 Due to the MAC grid structure
	int j = y / _DELTA_Y;
	int i_plus1 = i + 1;
	int j_plus1 = j + 1;
	double i_frac = x / _DELTA_X - 0.5 - i;
	double j_frac = y / _DELTA_Y - j;

	assert(i_frac <= 1);
	assert(j_frac <= 1);

	// Border cases
	i = CLAMP(i, 0, _SIZE_X - 1);
	j = CLAMP(j, 0, _SIZE_Y + 1 - 1);
	i_plus1 = CLAMP(i_plus1, 0, _SIZE_X - 1);
	j_plus1 = CLAMP(j_plus1, 0, _SIZE_Y + 1 - 1);
	
	// First interpolate in x, then in y
	double v_y_00 = _vel_y_front_buffer[i + j * _SIZE_X];
	double v_y_10 = _vel_y_front_buffer[i_plus1 + j * _SIZE_X];
	double v_y_01 = _vel_y_front_buffer[i + j_plus1 * _SIZE_X];
	double v_y_11 = _vel_y_front_buffer[i_plus1 + j_plus1 * _SIZE_X];
	
	// Interpolate x
	double v_y_0 = (1 - i_frac) * v_y_00 + i_frac * v_y_10;
	double v_y_1 = (1 - i_frac) * v_y_01 + i_frac * v_y_11; 

	// Interpolate y
	double v_y = (1 - j_frac) * v_y_0 + j_frac * v_y_1;
	return v_y;
}

double MacGrid::colorInterpolated(double x, double y) const
{
	// Calculate indices
	int i = x / _DELTA_X - 0.5;
	int j = y / _DELTA_Y - 0.5;
	int i_plus1 = i + 1;
	int j_plus1 = j + 1;
	double i_frac = x / _DELTA_X - 0.5 - i;
	double j_frac = y / _DELTA_Y - 0.5 - j;

	assert(i_frac <= 1);
	assert(j_frac <= 1);

	// Border cases
	i = CLAMP(i, 0, _SIZE_X - 1);
	j = CLAMP(j, 0, _SIZE_Y - 1);
	i_plus1 = CLAMP(i_plus1, 0, _SIZE_X - 1);
	j_plus1 = CLAMP(j_plus1, 0, _SIZE_Y - 1);
	
	// First interpolate in x, then in y
	double color_00 = _color_front_buffer[i + j * _SIZE_X];
	double color_10 = _color_front_buffer[i_plus1 + j * _SIZE_X];
	double color_01 = _color_front_buffer[i + j_plus1 * _SIZE_X];
	double color_11 = _color_front_buffer[i_plus1 + j_plus1 * _SIZE_X];
	
	// Interpolate x
	double color_0 = (1 - i_frac) * color_00 + i_frac * color_10;
	double color_1 = (1 - i_frac) * color_01 + i_frac * color_11; 

	// Interpolate y
	double color = (1 - j_frac) * color_0 + j_frac * color_1;
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
		return _cell_type_buffer[i + j * _SIZE_X];
	}
}

/**
	Sampling on borders yields LIQUID if the point neighbors a liquid cell
*/
CellType MacGrid::cellTypeXHalfIndexed(int x, int y) const
{
	// Look at the cells on either side of the point
	// (which is on the border of two cells)
	CellType left_type = cellType(x, y);
	CellType right_type = cellType(x - 1, y);

	return
		left_type == SOLID || right_type == SOLID ? SOLID :
		(left_type == LIQUID || right_type == LIQUID ? LIQUID :
		AIR);
}

/**
	Sampling on borders yields LIQUID if the point neighbors a liquid cell
*/
CellType MacGrid::cellTypeYHalfIndexed(int x, int y) const
{
	// Look at the cells on either side of the point
	// (which is on the border of two cells)
	CellType upper_type = cellType(x, y);
	CellType lower_type = cellType(x, y - 1);

	return
		upper_type == SOLID || lower_type == SOLID ? SOLID :
		(upper_type == LIQUID || lower_type == LIQUID ? LIQUID :
		AIR);
}

double MacGrid::divVelX(int x, int y) const
{
	return 
		(_vel_x_front_buffer[x + 1 + y * (_SIZE_X + 1)]
		- _vel_x_front_buffer[x + y * (_SIZE_X + 1)])
		/ _DELTA_X;
}

double MacGrid::divVelY(int x, int y) const
{
	return
		(_vel_y_front_buffer[x + (y + 1) * _SIZE_X]
		- _vel_y_front_buffer[x + y * _SIZE_X])
		/ _DELTA_Y;
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
void MacGrid::setVelX(int x, int y, double vel_x)
{
	_vel_x_front_buffer[x + y * (_SIZE_X + 1)] = vel_x;
}

/**
*/
void MacGrid::setVelY(int x, int y, double vel_y)
{
	_vel_y_front_buffer[x + y * _SIZE_X] = vel_y;
}

void MacGrid::setColor(int x, int y, double color)
{
	_color_front_buffer[x + y * _SIZE_X] = color;
}

void MacGrid::setCellType(int x, int y, CellType cell_type)
{
	_cell_type_buffer[x + y * _SIZE_X] = cell_type;
}

/**
	Writes to the four closest grid points, writes to back buffer.
*/
void MacGrid::addToVelXInterpolated(double x, double y, double vel_x)
{
	// Calculate indices
	int i = x / _DELTA_X;
	int j = y / _DELTA_Y - 0.5; // -0.5 Due to the MAC grid structure
	int i_plus1 = i + 1;
	int j_plus1 = j + 1;
	double i_frac = x / _DELTA_X - i;
	double j_frac = y / _DELTA_Y - 0.5 - j;

	assert(i_frac <= 1);
	assert(j_frac <= 1);

	// Border cases
	i = CLAMP(i, 0, _SIZE_X + 1 - 1);
	j = CLAMP(j, 0, _SIZE_Y - 1);
	i_plus1 = CLAMP(i_plus1, 0, _SIZE_X + 1 - 1);
	j_plus1 = CLAMP(j_plus1, 0, _SIZE_Y - 1);
	
	// Spread in y
	double v_x_0 = (1 - j_frac) * vel_x;
	double v_x_1 = j_frac * vel_x;

	// Spread in X
	double v_x_00 = (1 - i_frac) * v_x_0;
	double v_x_10 = i_frac * v_x_0;
	double v_x_01 = (1 - i_frac) * v_x_1;
	double v_x_11 = i_frac * v_x_1;

	// Write data
	_vel_x_back_buffer[i + j * (_SIZE_X + 1)] += v_x_00;
	_vel_x_back_buffer[i_plus1 + j * (_SIZE_X + 1)] += v_x_10;
	_vel_x_back_buffer[i + j_plus1 * (_SIZE_X + 1)] += v_x_01;
	_vel_x_back_buffer[i_plus1 + j_plus1 * (_SIZE_X + 1)] += v_x_11;
}

/**
	Writes to the four closest grid points, writes to back buffer.
*/
void MacGrid::addToVelYInterpolated(double x, double y, double vel_y)
{
	// Calculate indices
	int i = x / _DELTA_X - 0.5; // -0.5 Due to the MAC grid structure
	int j = y / _DELTA_Y;
	int i_plus1 = i + 1;
	int j_plus1 = j + 1;
	double i_frac = x / _DELTA_X - 0.5 - i;
	double j_frac = y / _DELTA_Y - j;

	assert(i_frac <= 1);
	assert(j_frac <= 1);

	// Border cases
	i = CLAMP(i, 0, _SIZE_X - 1);
	j = CLAMP(j, 0, _SIZE_Y + 1 - 1);
	i_plus1 = CLAMP(i_plus1, 0, _SIZE_X - 1);
	j_plus1 = CLAMP(j_plus1, 0, _SIZE_Y + 1 - 1);
	
	// Spread in y
	double v_y_0 = (1 - j_frac) * vel_y;
	double v_y_1 = j_frac * vel_y;

	// Spread in x
	double v_y_00 = (1 - i_frac) * v_y_0;
	double v_y_10 = i_frac * v_y_0;
	double v_y_01 = (1 - i_frac) * v_y_1;
	double v_y_11 = i_frac * v_y_1;

	// Write data
	_vel_y_back_buffer[i + j * _SIZE_X] += v_y_00;
	_vel_y_back_buffer[i_plus1 + j * _SIZE_X] += v_y_10;
	_vel_y_back_buffer[i + j_plus1 * _SIZE_X] += v_y_01;
	_vel_y_back_buffer[i_plus1 + j_plus1 * _SIZE_X] += v_y_11;
}

/**
	Writes to the four closest grid points, writes to back buffer.
*/
void MacGrid::addToColorInterpolated(double x, double y, double color)
{
	// Calculate indices
	int i = x / _DELTA_X - 0.5;
	int j = y / _DELTA_Y - 0.5;
	int i_plus1 = i + 1;
	int j_plus1 = j + 1;
	double i_frac = x / _DELTA_X - 0.5 - i;
	double j_frac = y / _DELTA_Y - 0.5 - j;

	assert(i_frac <= 1);
	assert(j_frac <= 1);

	// Border cases
	i = CLAMP(i, 0, _SIZE_X - 1);
	j = CLAMP(j, 0, _SIZE_Y - 1);
	i_plus1 = CLAMP(i_plus1, 0, _SIZE_X - 1);
	j_plus1 = CLAMP(j_plus1, 0, _SIZE_Y - 1);
		
	// Spread in y
	double color_0 = (1 - j_frac) * color;
	double color_1 = j_frac * color;

	// Spread in x
	double color_00 = (1 - i_frac) * color_0;
	double color_10 = i_frac * color_0;
	double color_01 = (1 - i_frac) * color_1;
	double color_11 = i_frac * color_1;
	
	// Write data
	_color_back_buffer[i + j * _SIZE_X] += color_00;
	_color_back_buffer[i_plus1 + j * _SIZE_X] += color_10;
	_color_back_buffer[i + j_plus1 * _SIZE_X] += color_01;
	_color_back_buffer[i_plus1 + j_plus1 * _SIZE_X] += color_11;
}


void MacGrid::_advectVelX(double dt)
{
	// Set all to 0
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X + 1; ++i)
		{
			// Write data
			_vel_x_back_buffer[i + j * (_SIZE_X + 1)] = 0;
		}
	}
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X + 1; ++i)
		{
			if (cellTypeXHalfIndexed(i, j) == LIQUID)
			{ // Only advect the liquid cells
				// Particle position
				double x_pos = i * _DELTA_X;
				double y_pos = (j + 0.5) * _DELTA_Y;

				// Velcoity
				double v_x = velXInterpolated(x_pos, y_pos);
				double v_y = velYInterpolated(x_pos, y_pos);

				// Previous particle position
				double x_pos_prev = x_pos - v_x * dt;
				double y_pos_prev = y_pos - v_y * dt;

				// Write data (in back buffer)
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
			_vel_y_back_buffer[i + j * _SIZE_X] = 0;
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

				// Velcoity
				double v_x = velXInterpolated(x_pos, y_pos);
				double v_y = velYInterpolated(x_pos, y_pos);

				// Previous particle position
				double x_pos_prev = x_pos - v_x * dt;
				double y_pos_prev = y_pos - v_y * dt;

				// Write data (in back buffer)
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
			_color_back_buffer[i + j * _SIZE_X] = 0;
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

				// Velcoity
				double v_x = velXInterpolated(x_pos, y_pos);
				double v_y = velYInterpolated(x_pos, y_pos);
				double color = colorInterpolated(x_pos, y_pos);

				// Previous particle position
				double x_pos_prev = x_pos - v_x * dt;
				double y_pos_prev = y_pos - v_y * dt;

				// Write data (in back buffer)
				addToColorInterpolated(x_pos_prev, y_pos_prev, color);
			}
		}
	}
}

void MacGrid::_swapBuffers()
{
	double* vel_x_tmp = _vel_x_front_buffer;
	double* vel_y_tmp = _vel_y_front_buffer;
	double* color_tmp = _color_front_buffer;

	_vel_x_front_buffer = _vel_x_back_buffer;
	_vel_y_front_buffer = _vel_y_back_buffer;
	_color_front_buffer = _color_back_buffer;

	_vel_x_back_buffer = vel_x_tmp;
	_vel_y_back_buffer = vel_y_tmp;
	_color_back_buffer = color_tmp;
}
