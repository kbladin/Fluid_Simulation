#include <FluidSolver.h>

FluidSolver::FluidSolver(int size_x, int size_y)
: _fluid_indices(size_x, size_y)
, _valid_mask(size_x, size_y)
, _valid_mask_back_buffer(size_x, size_y)
, _SIZE_X(size_x)
, _SIZE_Y(size_y)
{
	_cg_solver.setMaxIterations(15);
}

FluidSolver::~FluidSolver()
{

}

void FluidSolver::step(FluidDomain& fluid_domain, MyFloat dt)
{
	assert(fluid_domain.macGrid().sizeX() == _SIZE_X);
	assert(fluid_domain.macGrid().sizeY() == _SIZE_Y);
	// Self advection
	// (should be done on a divergence free field, therefore done first)
	advectVelocitySemiLagrangian(fluid_domain.macGrid(), dt);
	// Add external force
	addExternalAcceleration(fluid_domain.macGrid(), 0, 9.82, dt);
    // Solve boundary condition
	enforceDirichlet(fluid_domain.macGrid());
	// Ensure divergence free
	pressureSolve(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet(),
		fluid_domain.density());
    // Solve boundary condition again due to numerical errors in previous step
    enforceDirichlet(fluid_domain.macGrid());
    // Extend velocities outside of liquid cells so that liquid can flow
    extendVelocity(fluid_domain.macGrid());
}

void FluidSolver::addExternalForce(FluidDomain& fluid_domain, MyFloat F_x, MyFloat F_y, MyFloat dt)
{
	for (int j = 0; j < fluid_domain.macGrid().sizeY(); ++j)
	{
		for (int i = 0; i < fluid_domain.macGrid().sizeY(); ++i)
		{
			if (fluid_domain.macGrid().cellType(i, j) == LIQUID)
			{ // Only add force to the liquid cells
                fluid_domain.macGrid().setVelXHalfIndexed(i,j,
                	fluid_domain.macGrid().velXHalfIndexed(i,j)
                	+ F_x / fluid_domain.density() * dt);
                fluid_domain.macGrid().setVelYHalfIndexed(i,j,
                	fluid_domain.macGrid().velYHalfIndexed(i,j)
                	+ F_y / fluid_domain.density() * dt);
			}
		}
	}
}

void FluidSolver::addExternalAcceleration(MacGrid& mac_grid, MyFloat a_x, MyFloat a_y, MyFloat dt)
{
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeY(); ++i)
		{
			if (mac_grid.cellType(i, j) == LIQUID)
			{ // Only add force to the liquid cells
                mac_grid.setVelXHalfIndexed(i,j,
                	mac_grid.velXHalfIndexed(i,j) + a_x * dt);
                mac_grid.setVelYHalfIndexed(i,j,
                	mac_grid.velYHalfIndexed(i,j) + a_y * dt);
			}
		}
	}
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

void FluidSolver::pressureSolve(
	MacGrid& mac_grid,
	MarkerParticleSet& particle_set,
	MyFloat density)
{
	int n_fluid_cells = 0;
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			if (mac_grid.cellType(i, j) == LIQUID)
			{
				_fluid_indices(i,j) = n_fluid_cells;
				n_fluid_cells++;
			}
			else
			{
				_fluid_indices(i,j) = -1;
			}
		}
	}
    if (n_fluid_cells == 0) {
        return;
    }
/*
	// Loop through all particles
	for (auto it = particle_set.begin(); it != particle_set.end(); it++)
	{
		// Find the particles indices in the grid
		int x = (it->posX() / mac_grid.lengthX()) * mac_grid.sizeX();
		int y = (it->posY() / mac_grid.lengthY()) * mac_grid.sizeY();

		_n_particles(x, y)++;
	}
*/
	// Allocate matrix in which to store discrete laplacian
    A.resize(n_fluid_cells, n_fluid_cells);
    A.reserve(Eigen::VectorXi::Constant(5, n_fluid_cells));
    // Vector containing divergence
	VectorX b(n_fluid_cells);

	// Loop through all cells
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			if (_fluid_indices(i, j) != -1)
			{
				int idx = _fluid_indices(i, j);
				int n_non_solid_neighbors = 0;
				if (mac_grid.cellType(i - 1, j) != SOLID) {
					if (mac_grid.cellType(i - 1, j) == LIQUID) {
						A.insert(_fluid_indices(i - 1, j), idx) = 1 / pow(mac_grid.deltaX(), 2);
					}
					n_non_solid_neighbors++;
				}
				if (mac_grid.cellType(i + 1, j) != SOLID) {
					if (mac_grid.cellType(i + 1, j) == LIQUID) {
						A.insert(_fluid_indices(i + 1, j), idx) = 1 / pow(mac_grid.deltaX(), 2);
					}
					n_non_solid_neighbors++;
				}
				if (mac_grid.cellType(i, j - 1) != SOLID) {
					if (mac_grid.cellType(i, j - 1) == LIQUID) {
						A.insert(_fluid_indices(i, j - 1), idx) = 1 / pow(mac_grid.deltaX(), 2);
					}
					n_non_solid_neighbors++;
				}
				if (mac_grid.cellType(i, j + 1) != SOLID) {
					if (mac_grid.cellType(i, j + 1) == LIQUID) {
						A.insert(_fluid_indices(i, j + 1), idx) = 1 / pow(mac_grid.deltaX(), 2);
					}
					n_non_solid_neighbors++;
				}

				// Set diagonal value of matrix
				A.insert(idx, idx) = - n_non_solid_neighbors / pow(mac_grid.deltaX(), 2);

				// Calculate negative divergence and store in b
				b[idx] = mac_grid.divVelX(i, j) + mac_grid.divVelY(i, j);
			}
		}
	}

	// Vector containing pressures for each cell
	VectorX x(n_fluid_cells);
	_cg_solver.compute(A);
	x = _cg_solver.solve(b);

	// Loop throgh all cells to set new velocities from pressures
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			// Calculate indices
			int idx = _fluid_indices(i, j);
			int idx_i_minus1 = _fluid_indices(i - 1, j);
			int idx_j_minus1 = _fluid_indices(i, j - 1);
			
			if (idx >= 0 || idx_i_minus1 >= 0 || idx_j_minus1 >= 0)
			{
				// UGLY SOLUTION TO VOLUME LOSS HERE!!
				// MAKE PRESSURE PROPORTIONAL TO NUMBER OF PARTICLES
				//MyFloat k = 0;
                
				MyFloat p = idx >= 0 ? x(idx) : 0;

				MyFloat p_i_minus1 = idx_i_minus1 >= 0 ? x(idx_i_minus1) : 0;
				MyFloat p_j_minus1 = idx_j_minus1 >= 0 ? x(idx_j_minus1) : 0;

				// Get Pressure difference in x and y dimension
				MyFloat pressure_diff_x = p - p_i_minus1;
				MyFloat pressure_diff_y = p - p_j_minus1;

				// Current velocities
				MyFloat vel_x = mac_grid.velXHalfIndexed(i,j);
				MyFloat vel_y = mac_grid.velYHalfIndexed(i,j);

				// Calculate new velocity
				MyFloat new_vel_x = vel_x - 1.0 / density * pressure_diff_x / mac_grid.deltaX();
				MyFloat new_vel_y = vel_y - 1.0 / density * pressure_diff_y / mac_grid.deltaY();

				// Write data
				mac_grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);	
				mac_grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);
			}
		}
	}
	mac_grid.swapBuffers();
}

void FluidSolver::extendVelocity(MacGrid& mac_grid)
{
    for (int j = 0; j < mac_grid.sizeY(); ++j)
    {
        for (int i = 0; i < mac_grid.sizeX(); ++i)
        {
        	if(mac_grid.cellType(i, j) == LIQUID)
        	{
                _valid_mask(i,j) = 1;
                _valid_mask_back_buffer(i,j) = 1;
        	}
        	else
        	{
        		_valid_mask(i,j) = 0;
				_valid_mask_back_buffer(i,j) = 0;
        	}
            mac_grid.setVelXBackBufferHalfIndexed(i, j, mac_grid.velXHalfIndexed(i, j));
            mac_grid.setVelYBackBufferHalfIndexed(i, j, mac_grid.velYHalfIndexed(i, j));
        }
    }

    int iterations = 3;
    for (int iter = 0; iter < iterations; ++iter)
    {
        for (int j = 0; j < mac_grid.sizeY(); ++j)
        {
            for (int i = 0; i < mac_grid.sizeX(); ++i)
            {
            	if (_valid_mask.value(i, j) == 0 && mac_grid.cellType(i, j) != SOLID)
            	{
	            	MyFloat new_vel_x = 0;
	            	MyFloat new_vel_y = 0;
	            	int n_valid_neighbors = 0;

	            	// Get values of all neighbors
                    if(_valid_mask.value(i-1, j) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBuffer(i-1, j);
                        new_vel_y += mac_grid.velYBackBuffer(i-1, j);
                        n_valid_neighbors++;
                    }
                    if(_valid_mask.value(i, j-1) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBuffer(i, j-1);
                        new_vel_y += mac_grid.velYBackBuffer(i, j-1);
                        n_valid_neighbors++;
                    }
                    if(_valid_mask.value(i, j+1) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBuffer(i, j+1);
                        new_vel_y += mac_grid.velYBackBuffer(i, j+1);
                        n_valid_neighbors++;
                    }
                    if(_valid_mask.value(i+1, j) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBuffer(i+1, j);
                        new_vel_y += mac_grid.velYBackBuffer(i+1, j);
                        n_valid_neighbors++;
                    }

	            	// Average the value for the current cell
	            	if (n_valid_neighbors > 0)
	            	{
	            		new_vel_x /= n_valid_neighbors;
	            		new_vel_y /= n_valid_neighbors;

	            		mac_grid.setVelXBackBuffer(i, j, new_vel_x);
						mac_grid.setVelYBackBuffer(i, j, new_vel_y);
	            		_valid_mask_back_buffer(i, j) = 1;
	            	}
            	}
                /*
                if (_valid_mask.value(i, j) == 0 && mac_grid.cellType(i, j) != SOLID)
                {
                    MyFloat new_vel_x = 0;
                    MyFloat new_vel_y = 0;
                    int n_valid_neighbors = 0;
                    
                    // Get values of all neighbors
                    if(_valid_mask.value(i-1, j) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i-1, j);
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i-1, j);
                        n_valid_neighbors++;
                    }
                    if(_valid_mask.value(i+1, j) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i+1, j);
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i+1, j);
                        n_valid_neighbors++;
                    }
                    if(_valid_mask.value(i, j-1) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i, j-1);
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i, j-1);
                        n_valid_neighbors++;
                    }
                    if(_valid_mask.value(i, j+1) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i, j+1);
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i, j+1);
                        n_valid_neighbors++;
                    }
                    
                    // Average the value for the current cell
                    if (n_valid_neighbors > 0)
                    {
                        new_vel_x /= n_valid_neighbors;
                        new_vel_y /= n_valid_neighbors;
                        
                        mac_grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);
                        mac_grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);
                        _valid_mask_back_buffer(i, j) = 1;
                    }
                }*/

            }
        }
        // Swap valid mask
        Grid<char>* tmp = &_valid_mask;
        _valid_mask = std::move(_valid_mask_back_buffer);
        _valid_mask_back_buffer = std::move(*tmp);
    }
    mac_grid.swapBuffers();
}

void FluidSolver::advectVelocitySemiLagrangian(MacGrid& mac_grid, MyFloat dt)
{
	// Set all to 0
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			mac_grid.setVelXBackBufferHalfIndexed(i, j, 0);
			mac_grid.setVelYBackBufferHalfIndexed(i, j, 0);
		}
	}

	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			// X dimension
			if (mac_grid.cellType(i, j) == LIQUID ||
				mac_grid.cellType(i - 1, j) == LIQUID)
			{ // Only advect the liquid cells
				// Position in world coordinates
				MyFloat x_pos = i * mac_grid.deltaX();
				MyFloat y_pos = (j + 0.5) * mac_grid.deltaY();
				MyFloat x_pos_prev, y_pos_prev;
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
				MyFloat v_x = mac_grid.velXInterpolated(x_pos, y_pos);
				mac_grid.addToVelXInterpolated(x_pos_prev, y_pos_prev, v_x);
			}

			// Y dimension
			if (mac_grid.cellType(i, j) == LIQUID ||
				mac_grid.cellType(i, j - 1) == LIQUID)
			{ // Only advect the liquid cells
				// Position in world coordinates
				MyFloat x_pos = (i + 0.5) * mac_grid.deltaX();
				MyFloat y_pos = j * mac_grid.deltaY();
				MyFloat x_pos_prev, y_pos_prev;
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
				MyFloat v_y = mac_grid.velYInterpolated(x_pos, y_pos);
				mac_grid.addToVelYInterpolated(x_pos_prev, y_pos_prev, v_y);
			}
		}
	}
}

void FluidSolver::getAdvectedPositionRK3(
	MacGrid& mac_grid,
	MyFloat x_pos,
	MyFloat y_pos,
	MyFloat dt,
	MyFloat* x,
	MyFloat* y)
{
	MyFloat k_1x = mac_grid.velXInterpolated(x_pos, y_pos);
	MyFloat k_1y = mac_grid.velYInterpolated(x_pos, y_pos);
	
	MyFloat k_2x = mac_grid.velXInterpolated(
		x_pos + 1.0/2 * dt * k_1x,
		y_pos + 1.0/2 * dt * k_1y);
	MyFloat k_2y = mac_grid.velYInterpolated(
		x_pos + 1.0/2 * dt * k_1x,
		y_pos + 1.0/2 * dt * k_1y);

	MyFloat k_3x = mac_grid.velXInterpolated(
		x_pos + 3.0/4 * dt * k_2x,
		y_pos + 3.0/4 * dt * k_2y);
	MyFloat k_3y = mac_grid.velYInterpolated(
		x_pos + 3.0/4 * dt * k_2x,
		y_pos + 3.0/4 * dt * k_2y);

	MyFloat x_pos_prev = x_pos -
	dt * 1.0/9 * (
		2 * k_1x +
		3 * k_2x +
		4 * k_3x);
	MyFloat y_pos_prev = y_pos -
	dt * 1.0/6 * (
		2 * k_1y +
		3 * k_2y +
		4 * k_3y);

	*x = x_pos_prev;
	*y = y_pos_prev;
}

void FluidSolver::getAdvectedPositionForwardEuler(
	MacGrid& mac_grid,
	MyFloat x_pos,
	MyFloat y_pos,
	MyFloat dt,
	MyFloat* x,
	MyFloat* y)
{
	// Velcoity
	MyFloat v_x = mac_grid.velXInterpolated(x_pos, y_pos);
	MyFloat v_y = mac_grid.velYInterpolated(x_pos, y_pos);

	// Previous particle position
	MyFloat x_pos_prev = x_pos - v_x * dt;
	MyFloat y_pos_prev = y_pos - v_y * dt;

	*x = x_pos_prev;
	*y = y_pos_prev;
}
