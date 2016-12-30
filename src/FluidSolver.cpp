#include <FluidSolver.h>

FluidSolver::FluidSolver(int size_x, int size_y, MyFloat delta_x, MyFloat delta_y)
: _fluid_indices(size_x, size_y)
, _valid_mask_x(size_x, size_y)
, _valid_mask_x_back_buffer(size_x, size_y)
, _valid_mask_y(size_x, size_y)
, _valid_mask_y_back_buffer(size_x, size_y)
, _vel_x_sum(size_x, size_y, delta_x, delta_y)
, _vel_y_sum(size_x, size_y, delta_x, delta_y)
, _weight_vel_x_sum(size_x, size_y, delta_x, delta_y)
, _weight_vel_y_sum(size_x, size_y, delta_x, delta_y)
, _SIZE_X(size_x)
, _SIZE_Y(size_y)
, _DELTA_X(delta_x)
, _DELTA_Y(delta_y)
{
	_cg_solver.setMaxIterations(10);
}

FluidSolver::~FluidSolver()
{

}

void FluidSolver::stepSemiLagrangian(FluidDomain& fluid_domain, MyFloat dt)
{
	assert(fluid_domain.macGrid().sizeX() == _SIZE_X);
	assert(fluid_domain.macGrid().sizeY() == _SIZE_Y);
	// Classify the cells of the domain (AIR, LIQUID or SOLID)
	fluid_domain.classifyCells(fluid_domain.markerParticleSet());
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
    extendVelocity(fluid_domain.macGrid(), 2);
	
	// Advect the fluid through the newly updated field
	fluid_domain.advectParticlesWithGrid(dt);
}

void FluidSolver::stepPIC(FluidDomain& fluid_domain, MyFloat dt)
{
	assert(fluid_domain.macGrid().sizeX() == _SIZE_X);
	assert(fluid_domain.macGrid().sizeY() == _SIZE_Y);
	assert(fluid_domain.macGrid().deltaX() == _DELTA_X);
	assert(fluid_domain.macGrid().deltaY() == _DELTA_Y);

	// Classify the cells of the domain (AIR, LIQUID or SOLID)
	fluid_domain.classifyCells(fluid_domain.markerParticleSet());

	// Transfer to grid
	transferVelocityToGridSpread(fluid_domain.markerParticleSet(), fluid_domain.macGrid());
	// Resolve forces on grid
    addExternalAcceleration(fluid_domain.macGrid(), 0, 9.82, dt);
    enforceDirichlet(fluid_domain.macGrid());
	pressureSolve(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet(),
		fluid_domain.density());
    
    enforceDirichlet(fluid_domain.macGrid());
    extendVelocity(fluid_domain.macGrid(), 2);

    // Transfer to particles
	transferVelocityToParticlesPIC(fluid_domain.macGrid(), fluid_domain.markerParticleSet());

    // Advect
	fluid_domain.markerParticleSet().advect(dt);
}


void FluidSolver::stepFLIP(FluidDomain& fluid_domain, MyFloat dt)
{
	assert(fluid_domain.macGrid().sizeX() == _SIZE_X);
	assert(fluid_domain.macGrid().sizeY() == _SIZE_Y);
	assert(fluid_domain.macGrid().deltaX() == _DELTA_X);
	assert(fluid_domain.macGrid().deltaY() == _DELTA_Y);

	// Classify the cells of the domain (AIR, LIQUID or SOLID)
	fluid_domain.classifyCells(fluid_domain.markerParticleSet());

	// Transfer to grid
	transferVelocityToGridSpread(fluid_domain.markerParticleSet(), fluid_domain.macGrid());
	fluid_domain.macGrid().updatePreviousVelocityBuffer();
	// Resolve forces on grid
    addExternalAcceleration(fluid_domain.macGrid(), 0, 9.82, dt);
    enforceDirichlet(fluid_domain.macGrid());
	pressureSolve(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet(),
		fluid_domain.density());
    
    enforceDirichlet(fluid_domain.macGrid());
    //extendVelocity(fluid_domain.macGrid(), 2);

    fluid_domain.macGrid().updateVelocityDiffBuffer();

    // Transfer to particles
	transferVelocityToParticlesFLIP(fluid_domain.macGrid(), fluid_domain.markerParticleSet());

    // Advect
	fluid_domain.markerParticleSet().advect(dt);
}

void FluidSolver::stepPICFLIP(FluidDomain& fluid_domain, MyFloat dt, MyFloat pic_ratio)
{
	assert(fluid_domain.macGrid().sizeX() == _SIZE_X);
	assert(fluid_domain.macGrid().sizeY() == _SIZE_Y);
	assert(fluid_domain.macGrid().deltaX() == _DELTA_X);
	assert(fluid_domain.macGrid().deltaY() == _DELTA_Y);

	// Classify the cells of the domain (AIR, LIQUID or SOLID)
	fluid_domain.classifyCells(fluid_domain.markerParticleSet());

	// Transfer to grid
	transferVelocityToGridSpread(fluid_domain.markerParticleSet(), fluid_domain.macGrid());
	fluid_domain.macGrid().updatePreviousVelocityBuffer();
	// Resolve forces on grid
	addExternalAcceleration(fluid_domain.macGrid(), 0, 9.82, dt);
    enforceDirichlet(fluid_domain.macGrid());

    extendVelocity(fluid_domain.macGrid(), 2);
	pressureSolve(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet(),
		fluid_domain.density());
    enforceDirichlet(fluid_domain.macGrid());

    fluid_domain.macGrid().updateVelocityDiffBuffer();

    // Transfer to particles
	transferVelocityToParticlesPICFLIP(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet(),
		pic_ratio);

    // Advect
	fluid_domain.markerParticleSet().advect(dt);
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
            int i_minus1 = CLAMP(i - 1, 0, mac_grid.sizeX());
            int j_minus1 = CLAMP(j - 1, 0, mac_grid.sizeY());
			// X velocity
			if(	(mac_grid.cellType(i_minus1, j) == SOLID && mac_grid.velXHalfIndexed(i,j) < 0) ||
				(mac_grid.cellType(i, j)  == SOLID && mac_grid.velXHalfIndexed(i,j) > 0))
				mac_grid.setVelXHalfIndexed(i, j, 0);
			
			// Y velocity
			if(	(mac_grid.cellType(i, j_minus1) == SOLID && mac_grid.velYHalfIndexed(i,j) < 0) ||
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
    //A.reserve(Eigen::VectorXi::Constant(5, n_fluid_cells));
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
			int idx_i_minus1 = _fluid_indices(CLAMP(i - 1, 0, mac_grid.sizeX() - 1), j);
			int idx_j_minus1 = _fluid_indices(i, CLAMP(j - 1, 0, mac_grid.sizeY() - 1));
			
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

void FluidSolver::extendVelocity(MacGrid& mac_grid, int n_iterations)
{
    for (int j = 0; j < mac_grid.sizeY(); ++j)
    {
        for (int i = 0; i < mac_grid.sizeX(); ++i)
        {
        	if(mac_grid.cellType(i, j) == LIQUID || mac_grid.cellType(i - 1, j) == LIQUID)
        	{
                _valid_mask_x(i,j) = 1;
                _valid_mask_x_back_buffer(i,j) = 1;
	            mac_grid.setVelXBackBufferHalfIndexed(i, j, mac_grid.velXHalfIndexed(i, j));
	            mac_grid.setVelXHalfIndexed(i, j, mac_grid.velXHalfIndexed(i, j));
        	}
        	else
        	{
        		_valid_mask_x(i,j) = 0;
				_valid_mask_x_back_buffer(i,j) = 0;
	            mac_grid.setVelXBackBufferHalfIndexed(i, j, 0);
	            mac_grid.setVelXHalfIndexed(i, j, 0);
        	}

        	if(mac_grid.cellType(i, j) == LIQUID || mac_grid.cellType(i, j - 1) == LIQUID)
        	{
                _valid_mask_y(i,j) = 1;
                _valid_mask_y_back_buffer(i,j) = 1;
	            mac_grid.setVelYBackBufferHalfIndexed(i, j, mac_grid.velYHalfIndexed(i, j));
	            mac_grid.setVelYHalfIndexed(i, j, mac_grid.velYHalfIndexed(i, j));
        	}
        	else
        	{
        		_valid_mask_y(i,j) = 0;
				_valid_mask_y_back_buffer(i,j) = 0;
	            mac_grid.setVelYBackBufferHalfIndexed(i, j, 0);
	            mac_grid.setVelXHalfIndexed(i, j, 0);
        	}
        }
    }

    for (int iter = 0; iter < n_iterations; ++iter)
    {
        for (int j = 0; j < mac_grid.sizeY(); ++j)
        {
            for (int i = 0; i < mac_grid.sizeX(); ++i)
            {
            	if (_valid_mask_x.value(i, j) == 0 && mac_grid.cellType(i, j) != SOLID && mac_grid.cellType(i - 1, j) != SOLID)
            	{
	            	MyFloat new_vel_x = 0;
	            	int n_valid_neighbors_x = 0;

	            	// Get values of all neighbors
                    if(_valid_mask_x.value(i-1, j) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i-1, j);
                        n_valid_neighbors_x++;
                    }
                    if(_valid_mask_x.value(i, j-1) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i, j-1);
                        n_valid_neighbors_x++;
                    }
                    if(_valid_mask_x.value(i, j+1) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i, j+1);
                        n_valid_neighbors_x++;
                    }
                    if(_valid_mask_x.value(i+1, j) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i+1, j);
                        n_valid_neighbors_x++;
                    }

	            	// Average the value for the current cell
	            	if (n_valid_neighbors_x > 0)
	            	{
	            		new_vel_x /= n_valid_neighbors_x;
	            	
	            		mac_grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);
						_valid_mask_x_back_buffer(i, j) = 1;
	            	}
            	}

            	if (_valid_mask_y.value(i, j) == 0 && mac_grid.cellType(i, j) != SOLID && mac_grid.cellType(i, j - 1) != SOLID)
            	{
	            	MyFloat new_vel_y = 0;
	            	int n_valid_neighbors_y = 0;

	            	// Get values of all neighbors
                    if(_valid_mask_y.value(i-1, j) == 1)
                    {
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i-1, j);
                        n_valid_neighbors_y++;
                    }
                    if(_valid_mask_y.value(i, j-1) == 1)
                    {
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i, j-1);
                        n_valid_neighbors_y++;
                    }
                    if(_valid_mask_y.value(i, j+1) == 1)
                    {
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i, j+1);
                        n_valid_neighbors_y++;
                    }
                    if(_valid_mask_y.value(i+1, j) == 1)
                    {
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i+1, j);
                        n_valid_neighbors_y++;
                    }

	            	// Average the value for the current cell
	            	if (n_valid_neighbors_y > 0)
	            	{
	            		new_vel_y /= n_valid_neighbors_y;
	            	
	            		mac_grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);
						_valid_mask_y_back_buffer(i, j) = 1;
	            	}
            	}
            }
        }
        // Swap valid mask
        Grid<char>* tmp = &_valid_mask_x;
        _valid_mask_x = std::move(_valid_mask_x_back_buffer);
        _valid_mask_x_back_buffer = std::move(*tmp);

        Grid<char>* tmp2 = &_valid_mask_y;
        _valid_mask_y = std::move(_valid_mask_y_back_buffer);
        _valid_mask_y_back_buffer = std::move(*tmp2);
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

void FluidSolver::transferVelocityToGridGather(
	MarkerParticleSet& particle_set,
	MacGrid& mac_grid)
{
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			MyFloat x_x_vel_grid = i * mac_grid.deltaX();
			MyFloat y_x_vel_grid = (j + 0.5) * mac_grid.deltaY();
			MyFloat x_y_vel_grid = (i + 0.5) * mac_grid.deltaX();
			MyFloat y_y_vel_grid = j * mac_grid.deltaY();
			
			MyFloat weight_vel_x_sum = 0;
			MyFloat vel_x_sum = 0;
			MyFloat weight_vel_y_sum = 0;
			MyFloat vel_y_sum = 0;
			for (auto p = particle_set.begin(); p != particle_set.end(); p++)
			{
				// x velocity
				MyFloat abs_x_diff = std::abs(p->posX() - x_x_vel_grid);
				MyFloat abs_y_diff = std::abs(p->posY() - y_x_vel_grid);
				MyFloat weight_vel_x =
					1 - (abs_x_diff / mac_grid.deltaX() +
					abs_y_diff / mac_grid.deltaY());
				if (weight_vel_x >= 1)
				{
					weight_vel_x_sum += weight_vel_x;
					vel_x_sum += p->velX();
				}
				// y velocity
				abs_x_diff = std::abs(p->posX() - x_y_vel_grid);
				abs_y_diff = std::abs(p->posY() - y_y_vel_grid);
				MyFloat weight_vel_y =
					1 - (abs_x_diff / mac_grid.deltaX() +
					abs_y_diff / mac_grid.deltaY());
				if (weight_vel_y >= 1)
				{
					weight_vel_y_sum += weight_vel_y;
					vel_y_sum += p->velY();
				}
			}
			if (weight_vel_x_sum)
			{
				vel_x_sum /= weight_vel_x_sum;
				mac_grid.setVelXBackBufferHalfIndexed(i, j, vel_x_sum);
			}
			if (weight_vel_y_sum)
			{
				vel_y_sum /= weight_vel_y_sum;
				mac_grid.setVelYBackBufferHalfIndexed(i, j, vel_y_sum);
			}
		}
	}
	mac_grid.swapBuffers();
}


void FluidSolver::transferVelocityToGridSpread(
	MarkerParticleSet& particle_set,
	MacGrid& mac_grid)
{
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			_vel_x_sum(i, j) = 0;
			_vel_y_sum(i, j) = 0;
			_weight_vel_x_sum(i, j) = 0;
			_weight_vel_y_sum(i, j) = 0;
		}
	}
	for (auto p = particle_set.begin(); p != particle_set.end(); p++)
	{
		_vel_x_sum.addToValueInterpolated(
			p->posX(), p->posY() - 0.5 * mac_grid.deltaY(), p->velX());
		_vel_y_sum.addToValueInterpolated(
			p->posX() - 0.5 * mac_grid.deltaX(), p->posY(), p->velY());
		_weight_vel_x_sum.addToValueInterpolated(
			p->posX(), p->posY() - 0.5 * mac_grid.deltaY(), 1.0);
		_weight_vel_y_sum.addToValueInterpolated(
			p->posX() - 0.5 * mac_grid.deltaX(), p->posY(), 1.0);
	}
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			if (_weight_vel_x_sum(i, j) > 0.000001)
			{
				MyFloat new_vel_x = _vel_x_sum(i, j) / _weight_vel_x_sum(i, j);
				mac_grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);
			}
			if (_weight_vel_y_sum(i, j) > 0.000001)
			{
				MyFloat new_vel_y = _vel_y_sum(i, j) / _weight_vel_y_sum(i, j);
				mac_grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);
			}
		}
	}
	mac_grid.swapBuffers();
}

void FluidSolver::transferVelocityToParticlesPIC(
	MacGrid& mac_grid,
	MarkerParticleSet& particle_set)
{
	for (auto p = particle_set.begin(); p != particle_set.end(); p++)
	{
		p->setVelocity(
			mac_grid.velXInterpolated(p->posX(), p->posY()),
			mac_grid.velYInterpolated(p->posX(), p->posY()));
	}
}


void FluidSolver::transferVelocityToParticlesFLIP(
	MacGrid& mac_grid,
	MarkerParticleSet& particle_set)
{
	for (auto p = particle_set.begin(); p != particle_set.end(); p++)
	{
        p->setVelocity(
			p->velX() + mac_grid.velXDiffInterpolated(p->posX(), p->posY()),
			p->velY() + mac_grid.velYDiffInterpolated(p->posX(), p->posY()));
	}
}

void FluidSolver::transferVelocityToParticlesPICFLIP(
	MacGrid& mac_grid,
	MarkerParticleSet& particle_set,
	MyFloat pic_ratio)
{
	for (auto p = particle_set.begin(); p != particle_set.end(); p++)
	{
        MyFloat pic_vel_x = mac_grid.velXInterpolated(p->posX(), p->posY());
        MyFloat pic_vel_y = mac_grid.velYInterpolated(p->posX(), p->posY());
        MyFloat flip_vel_x = p->velX() + mac_grid.velXDiffInterpolated(p->posX(), p->posY());
        MyFloat flip_vel_y = p->velY() + mac_grid.velYDiffInterpolated(p->posX(), p->posY());

        p->setVelocity(
        	pic_vel_x * pic_ratio + flip_vel_x * (1 - pic_ratio),
        	pic_vel_y * pic_ratio + flip_vel_y * (1 - pic_ratio));
	}
}
