#include <FluidSolver.h>

static const int N_PARTICLES = 10000;

FluidSolver::FluidSolver()
{
	_cg_solver.setMaxIterations(5);
}

FluidSolver::~FluidSolver()
{

}

void FluidSolver::step(FluidDomain& fluid_domain, double dt)
{	
	// Self advection
	advectVelocity(fluid_domain.macGrid(), dt);
	// Solve boundary condition
	enforceDirichlet(fluid_domain.macGrid());
	// Ensure divergence free
	pressureSolve(fluid_domain.macGrid(), dt);
	// Solve boundary condition again due to numerical errors in previous step
	enforceDirichlet(fluid_domain.macGrid());
	// Extend velocities outside of liquid cells so that liquid can flow
	extendVelocity(fluid_domain.macGrid());
}

void FluidSolver::enforceDirichlet(MacGrid& mac_grid)
{
	// X vel
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			// X velocity
			if(	(mac_grid.cellType(i - 1, j) == SOLID && mac_grid.velXHalfIndexed(i,j) < 0) ||
				(mac_grid.cellType(i, j)  == SOLID && mac_grid.velXHalfIndexed(i,j) > 0))
				mac_grid.setVelXHalfIndexed(i, j, 0);
			
			// Y velocity
			if(	(mac_grid.cellType(i, j - 1) == SOLID && mac_grid.velYHalfIndexed(i,j) < 0) ||
				(mac_grid.cellType(i, j) == SOLID && mac_grid.velYHalfIndexed(i,j) > 0))
				mac_grid.setVelYHalfIndexed(i, j, 0);
		}
	}
}


void FluidSolver::pressureSolve(MacGrid& mac_grid, double dt)
{
	// Find which cells are liquid by using a grid with indices
	Grid<int> fluid_indices(mac_grid.sizeX(), mac_grid.sizeY());
	int n_fluid_cells = 0;
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			if (mac_grid.cellType(i, j) == LIQUID)
			{
				fluid_indices(i,j) = n_fluid_cells;
				n_fluid_cells++;
			}
			else
			{
				fluid_indices(i,j) = -1;
			}
		}
	}

	// Allocate matrix in which to store connectivity information
    Eigen::SparseMatrix<double> A(n_fluid_cells, n_fluid_cells);
	// Vector containing negative divergence
	Eigen::VectorXd b(n_fluid_cells);

	// Loop through all cells
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			if (mac_grid.cellType(i, j) == LIQUID)
			{
				int idx = fluid_indices(i, j);
				int n_non_solid_neighbors = 0;
				if (mac_grid.cellType(i - 1, j) != SOLID) {
					if (mac_grid.cellType(i - 1, j) == LIQUID) {
						A.insert(fluid_indices(i - 1, j), idx) = 1;
					}
					n_non_solid_neighbors++;
				}
				if (mac_grid.cellType(i + 1, j) != SOLID) {
					if (mac_grid.cellType(i + 1, j) == LIQUID) {
						A.insert(fluid_indices(i + 1, j), idx) = 1;
					}
					n_non_solid_neighbors++;
				}
				if (mac_grid.cellType(i, j - 1) != SOLID) {
					if (mac_grid.cellType(i, j - 1) == LIQUID) {
						A.insert(fluid_indices(i, j - 1), idx) = 1;
					}
					n_non_solid_neighbors++;
				}
				if (mac_grid.cellType(i, j + 1) != SOLID) {
					if (mac_grid.cellType(i, j + 1) == LIQUID) {
						A.insert(fluid_indices(i, j + 1), idx) = 1;
					}
					n_non_solid_neighbors++;
				}

				// Set diagonal value of matrix
				A.insert(idx, idx) = - n_non_solid_neighbors;

				// Calculate negative divergence and store in b
				b[idx] = - (mac_grid.divVelX(i, j) + mac_grid.divVelY(i, j));
			}
		}
	}

	// Vector containing pressures for each cell
	Eigen::VectorXd x(n_fluid_cells);
	_cg_solver.compute(A);
	x = _cg_solver.solve(b);

	double density = 0.1;

	// Loop throgh all cells to set new velocities from pressures
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			// Calculate indices
			int idx = fluid_indices(i, j);
			int idx_i_minus1 = fluid_indices(i - 1, j);
			int idx_j_minus1 = fluid_indices(i, j - 1);
			
			if (idx > 0)
			{
				double p = x(idx); // Pressure
				// X dimension
				if (idx_i_minus1 > 0)
				{
					// When looking at x velocities at the borders of the cells,
					// the pressure 'p' is to the right and p_i_minus1 is to
					// the left of the border
					double p_i_minus1 = idx_i_minus1 > 0 ? x(idx_i_minus1) : 0;

					// Get Pressure difference in x and y dimension
					double pressure_diff_x = p - p_i_minus1;
					
					// Current velocities
					double vel_x = mac_grid.velXHalfIndexed(i,j);
					
					// Calculate new velocity
					double new_vel_x = vel_x + dt * 1 / density * pressure_diff_x * mac_grid.deltaX();
					// Write data
					mac_grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);	
				}
				// Y dimension
				if (idx_j_minus1 > 0)
				{
					// When looking at y velocities at the borders of the cells,
					// the pressure 'p' above and p_i_minus1 is below the border
					double p_j_minus1 = idx_j_minus1 > 0 ? x(idx_j_minus1) : 0;

					// Get Pressure difference in x and y dimension
					double pressure_diff_y = p - p_j_minus1;

					// Current velocities
					double vel_y = mac_grid.velYHalfIndexed(i,j);

					// Calculate new velocity
					double new_vel_y = vel_y + dt * 1 / density * pressure_diff_y * mac_grid.deltaY();
					// Write data
					mac_grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);	
				}
			}
		}
	}
	mac_grid.swapBuffers();
}

/**
	Extends the velocity of the mac grid from the interface of the level set.
*/
void FluidSolver::extendVelocity(MacGrid& mac_grid)
{
	Grid<int> valid_mask(mac_grid.sizeX(), mac_grid.sizeY());
    for (int j = 0; j < mac_grid.sizeY(); ++j)
    {
        for (int i = 0; i < mac_grid.sizeX(); ++i)
        {
        	if(mac_grid.cellType(i, j) == LIQUID)
        	{
        		valid_mask(i,j) = 1;
				mac_grid.setVelXBackBufferHalfIndexed(i, j, mac_grid.velXHalfIndexed(i, j));
	            mac_grid.setVelYBackBufferHalfIndexed(i, j, mac_grid.velYHalfIndexed(i, j));
        	}
        }
    }

    int iterations = 5;
    for (int iter = 0; iter < iterations; ++iter)
    {
        for (int j = 0; j < mac_grid.sizeY(); ++j)
        {
            for (int i = 0; i < mac_grid.sizeX(); ++i)
            {
            	if (valid_mask.value(i, j) == 0)
            	{
	            	double new_vel_x = 0;
	            	double new_vel_y = 0;
	            	int n_valid_neighbors = 0;

	            	// Get values of all neighbors
	            	if(valid_mask.value(i-1, j) == 1)
	            	{
	            		new_vel_x += mac_grid.velX(i-1, j);
	            		new_vel_y += mac_grid.velY(i-1, j);
	            		n_valid_neighbors++;
	            	}
	            	if(valid_mask.value(i+1, j) == 1)
	            	{
	            		new_vel_x += mac_grid.velX(i+1, j);
	            		new_vel_y += mac_grid.velY(i+1, j);
	            		n_valid_neighbors++;
	            	}
	            	if(valid_mask.value(i, j-1) == 1)
	            	{
	            		new_vel_x += mac_grid.velX(i, j-1);
	            		new_vel_y += mac_grid.velY(i, j-1);
	            		n_valid_neighbors++;
	            	}
	            	if(valid_mask.value(i, j+1) == 1)
	            	{
	            		new_vel_x += mac_grid.velX(i, j+1);
	            		new_vel_y += mac_grid.velY(i, j+1);
	            		n_valid_neighbors++;
	            	}

	            	// Average the value for the current cell
	            	if (n_valid_neighbors > 0)
	            	{
	            		new_vel_x /= n_valid_neighbors;
	            		new_vel_y /= n_valid_neighbors;

	            		mac_grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);
						mac_grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);
	            		valid_mask(i, j) = 1;
	            	}
            	}
            }
        }
    }
    mac_grid.swapBuffers();
}

void FluidSolver::advectVelocity(MacGrid& mac_grid, double dt)
{
	// Set all to 0
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			mac_grid.setVelXBackBuffer(i, j, 0);
			mac_grid.setVelYBackBuffer(i, j, 0);
		}
	}


	// Merge these loops later




	// X
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			if (mac_grid.cellType(i, j) == LIQUID ||
				mac_grid.cellType(i - 1, j) == LIQUID)
			{ // Only advect the liquid cells
				// Position in world coordinates
				double x_pos = i * mac_grid.deltaX();
				double y_pos = (j + 0.5) * mac_grid.deltaY();
				double x_pos_prev, y_pos_prev;
//#define USE_EULER
#ifdef USE_EULER
				getAdvectedPositionForwardEuler(
					mac_grid,
					x_pos,
					y_pos,
					dt,
					&x_pos_prev,
					&y_pos_prev);
#else
				getAdvectedPositionRK3(
					mac_grid,
					x_pos,
					y_pos,
					dt,
					&x_pos_prev,
					&y_pos_prev);
#endif
				// Velocity to be stored
				double v_x = mac_grid.velXInterpolated(x_pos, y_pos);
				mac_grid.addToVelXInterpolated(x_pos_prev, y_pos_prev, v_x);
			}
		}
	}









	// Y
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			if (mac_grid.cellType(i, j) == LIQUID ||
				mac_grid.cellType(i, j - 1) == LIQUID)
			{ // Only advect the liquid cells
				// Position in world coordinates
				double x_pos = (i + 0.5) * mac_grid.deltaX();
				double y_pos = j * mac_grid.deltaY();
				double x_pos_prev, y_pos_prev;
#ifdef USE_EULER
				getAdvectedPositionForwardEuler(
					mac_grid,
					x_pos,
					y_pos,
					dt,
					&x_pos_prev,
					&y_pos_prev);
#else
				getAdvectedPositionRK3(
					mac_grid,
					x_pos,
					y_pos,
					dt,
					&x_pos_prev,
					&y_pos_prev);
#endif
				// Velocity to be stored
				double v_y = mac_grid.velYInterpolated(x_pos, y_pos);
				mac_grid.addToVelYInterpolated(x_pos_prev, y_pos_prev, v_y);
			}
		}
	}
}

void FluidSolver::getAdvectedPositionRK3(
	MacGrid& mac_grid,
	double x_pos,
	double y_pos,
	double dt,
	double* x,
	double* y)
{
	double k_1x = mac_grid.velXInterpolated(x_pos, y_pos);
	double k_1y = mac_grid.velYInterpolated(x_pos, y_pos);
	
	double k_2x = mac_grid.velXInterpolated(
		x_pos + 1.0/2 * dt * k_1x,
		y_pos + 1.0/2 * dt * k_1y);
	double k_2y = mac_grid.velYInterpolated(
		x_pos + 1.0/2 * dt * k_1x,
		y_pos + 1.0/2 * dt * k_1y);

	double k_3x = mac_grid.velXInterpolated(
		x_pos + 3.0/4 * dt * k_2x,
		y_pos + 3.0/4 * dt * k_2y);
	double k_3y = mac_grid.velYInterpolated(
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

void FluidSolver::getAdvectedPositionForwardEuler(
	MacGrid& mac_grid,
	double x_pos,
	double y_pos,
	double dt,
	double* x,
	double* y)
{
	// Velcoity
	double v_x = mac_grid.velXInterpolated(x_pos, y_pos);
	double v_y = mac_grid.velYInterpolated(x_pos, y_pos);

	// Previous particle position
	double x_pos_prev = x_pos - v_x * dt;
	double y_pos_prev = y_pos - v_y * dt;

	*x = x_pos_prev;
	*y = y_pos_prev;
}
