#include <FluidSolver.h>

/*
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/sequence.h>

#include <cusparse.h>
#include <cublas_v2.h>

template<typename T>
struct Fun
{
    __device__ T operator()(T t1, T t2)
    {
        auto result = t1+t2;
        return result;
    }
};

int run()
{
    const int N = 100;
    thrust::device_vector<int> vec(N);
    thrust::sequence(vec.begin(),vec.end());
    auto op = Fun<int>();
    return thrust::reduce(vec.begin(),vec.end(),0,op);
}
*/

FluidSolverMemoryPool::FluidSolverMemoryPool(
	int size_x,
	int size_y,
	MyFloat delta_x,
	MyFloat delta_y) :
	GridInterface(size_x, size_y, delta_x, delta_y),
	fluid_indices(size_x, size_y),
	n_particles(size_x, size_y),
	vel_x_sum(size_x, size_y, delta_x, delta_y),
	vel_y_sum(size_x, size_y, delta_x, delta_y),
	weight_vel_x_sum(size_x, size_y, delta_x, delta_y),
	weight_vel_y_sum(size_x, size_y, delta_x, delta_y)
{
	valid_mask_x_front_buffer =
		std::make_unique<Grid <unsigned char> >(size_x, size_y);
	valid_mask_x_back_buffer =
		std::make_unique<Grid <unsigned char> >(size_x, size_y);
	valid_mask_y_front_buffer =
		std::make_unique<Grid <unsigned char> >(size_x, size_y);
	valid_mask_y_back_buffer =
		std::make_unique<Grid <unsigned char> >(size_x, size_y);
}

FluidSolverMemoryPool::FluidSolverMemoryPool(
	const FluidSolverMemoryPool& other) :
	FluidSolverMemoryPool(
		other.sizeX(),
		other.sizeY(),
		other.deltaX(),
		other.deltaX())
{

}

FluidSolverMemoryPool::~FluidSolverMemoryPool()
{

}

void FluidSolverMemoryPool::swapValidMaskBuffer()
{
    valid_mask_x_front_buffer.swap(valid_mask_x_back_buffer);
    valid_mask_y_front_buffer.swap(valid_mask_y_back_buffer);
}

FluidSolver::FluidSolver()
{
	_cg_solver.setMaxIterations(10);
}

FluidSolver::~FluidSolver()
{

}

void FluidSolver::initialize(FluidSolverMemoryPool mem_pool)
{
	_mem_pool = std::make_unique<FluidSolverMemoryPool>(mem_pool);
}

void FluidSolver::stepSemiLagrangian(FluidDomain& fluid_domain, MyFloat dt)
{
	assert(_mem_pool);

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
		fluid_domain.density(),
		dt);
    // Solve boundary condition again due to numerical errors in previous step
    enforceDirichlet(fluid_domain.macGrid());
	// Extend velocities outside of liquid cells so that liquid can flow
    extendVelocityAvarageing(fluid_domain.macGrid(), 2);
	// Advect the fluid through the newly updated field
	fluid_domain.advectParticlesWithGrid(dt);
}

void FluidSolver::stepPIC(FluidDomain& fluid_domain, MyFloat dt)
{
	assert(_mem_pool);

	// Classify the cells of the domain (AIR, LIQUID or SOLID)
	fluid_domain.classifyCells(fluid_domain.markerParticleSet());
	// Transfer to grid
	transferVelocityToGridSpread(
		fluid_domain.markerParticleSet(),
		fluid_domain.macGrid());
	
	// Resolve forces on grid
    addExternalAcceleration(fluid_domain.macGrid(), 0, 9.82, dt);
    enforceDirichlet(fluid_domain.macGrid());
	pressureSolve(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet(),
		fluid_domain.density(),
		dt);
    enforceDirichlet(fluid_domain.macGrid());
    extendVelocityIndividual(fluid_domain.macGrid(), 2);

    // Transfer to particles
	transferVelocityToParticlesPIC(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet());
    // Advect
	fluid_domain.markerParticleSet().advect(dt);
}


void FluidSolver::stepFLIP(FluidDomain& fluid_domain, MyFloat dt)
{
	assert(_mem_pool);

	// Classify the cells of the domain (AIR, LIQUID or SOLID)
	fluid_domain.classifyCells(fluid_domain.markerParticleSet());
	// Transfer to grid
	transferVelocityToGridSpread(
		fluid_domain.markerParticleSet(),
		fluid_domain.macGrid());
	// Save velocity buffer
	fluid_domain.macGrid().updatePreviousVelocityBuffer();
	// Resolve forces on grid
    addExternalAcceleration(fluid_domain.macGrid(), 0, 9.82, dt);
    enforceDirichlet(fluid_domain.macGrid());
	pressureSolve(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet(),
		fluid_domain.density(),
		dt);
    enforceDirichlet(fluid_domain.macGrid());
    extendVelocityIndividual(fluid_domain.macGrid(), 2);
    // Save velocity difference
    fluid_domain.macGrid().updateVelocityDiffBuffer();
    // Transfer to particles
	transferVelocityToParticlesFLIP(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet());
    // Advect
	fluid_domain.markerParticleSet().advect(dt);
}

void FluidSolver::stepPICFLIP(
	FluidDomain& fluid_domain,
	MyFloat dt,
	MyFloat pic_ratio)
{
	assert(_mem_pool);

	// Classify the cells of the domain (AIR, LIQUID or SOLID)
	fluid_domain.classifyCells(fluid_domain.markerParticleSet());
	// Transfer to grid
	transferVelocityToGridSpread(
		fluid_domain.markerParticleSet(),
		fluid_domain.macGrid());
	// Save velocity buffer
	fluid_domain.macGrid().updatePreviousVelocityBuffer();
	// Resolve forces on grid
	addExternalAcceleration(fluid_domain.macGrid(), 0, 9.82, dt);
    enforceDirichlet(fluid_domain.macGrid());
    extendVelocityIndividual(fluid_domain.macGrid(), 2);
	pressureSolve(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet(),
		fluid_domain.density(),
		dt);
    enforceDirichlet(fluid_domain.macGrid());
    // Save velocity difference
    fluid_domain.macGrid().updateVelocityDiffBuffer();
    // Transfer to particles
	transferVelocityToParticlesPICFLIP(
		fluid_domain.macGrid(),
		fluid_domain.markerParticleSet(),
		pic_ratio);
    // Advect
	fluid_domain.markerParticleSet().advect(dt);
}

void FluidSolver::addExternalForce(
	FluidDomain& fluid_domain,
	MyFloat F_x,
	MyFloat F_y,
	MyFloat dt)
{
	for (int j = 0; j < fluid_domain.macGrid().sizeY(); ++j)
	{
		for (int i = 0; i < fluid_domain.macGrid().sizeX(); ++i)
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

void FluidSolver::addExternalAcceleration(
	MacGrid& mac_grid,
	MyFloat a_x,
	MyFloat a_y,
	MyFloat dt)
{
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
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
            int i_minus1 = CLAMP(i - 1, 0, mac_grid.sizeX() - 1);
            int j_minus1 = CLAMP(j - 1, 0, mac_grid.sizeY() - 1);
			// X velocity
			if(	(mac_grid.cellType(i_minus1, j) == SOLID &&
				mac_grid.velXHalfIndexed(i,j) < 0) ||
				(mac_grid.cellType(i, j)  == SOLID &&
				mac_grid.velXHalfIndexed(i,j) > 0))
				mac_grid.setVelXHalfIndexed(i, j, 0);
			
			// Y velocity
			if(	(mac_grid.cellType(i, j_minus1) == SOLID &&
				mac_grid.velYHalfIndexed(i,j) < 0) ||
				(mac_grid.cellType(i, j) == SOLID &&
				mac_grid.velYHalfIndexed(i,j) > 0))
				mac_grid.setVelYHalfIndexed(i, j, 0);
		}
	}
}

void FluidSolver::pressureSolve(
	MacGrid& mac_grid,
	MarkerParticleSet& particle_set,
	MyFloat density,
	MyFloat dt)
{
	int n_fluid_cells = 0;
	// Index all cells with fluid
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			if (mac_grid.cellType(i, j) == LIQUID)
			{
				_mem_pool->fluid_indices(i,j) = n_fluid_cells;
				n_fluid_cells++;
			}
			else
			{
				_mem_pool->fluid_indices(i,j) = -1;
			}
            _mem_pool->n_particles(i, j) = 0;
		}
	}
    if (n_fluid_cells == 0)
    {
        return;
    }

	// Loop through all particles to count the number of particles in each cell
	for (auto it = particle_set.begin(); it != particle_set.end(); it++)
	{
		// Find the particles indices in the grid
		int x = (it->posX() / mac_grid.lengthX()) * mac_grid.sizeX();
		int y = (it->posY() / mac_grid.lengthY()) * mac_grid.sizeY();
        
		_mem_pool->n_particles(x, y)++;
	}

	// Allocate matrix in which to store discrete laplacian
    _A.resize(n_fluid_cells, n_fluid_cells);
    _A.reserve(Eigen::VectorXi::Constant(5, n_fluid_cells));
    // Vector containing divergence
	VectorX b(n_fluid_cells);

	// Loop through all cells
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			if (_mem_pool->fluid_indices(i, j) != -1)
			{
				// This is a cell with fluid
				// The decrete Laplacian is computed for all liquid cells
				int idx = _mem_pool->fluid_indices(i, j);
				int n_non_solid_neighbors = 0;
				if (mac_grid.cellType(i - 1, j) != SOLID) {
					if (mac_grid.cellType(i - 1, j) == LIQUID) {
						_A.insert(_mem_pool->fluid_indices(i - 1, j), idx) =
							1 / pow(mac_grid.deltaX(), 2);
					}
					n_non_solid_neighbors++;
				}
				if (mac_grid.cellType(i + 1, j) != SOLID) {
					if (mac_grid.cellType(i + 1, j) == LIQUID) {
						_A.insert(_mem_pool->fluid_indices(i + 1, j), idx) =
							1 / pow(mac_grid.deltaX(), 2);
					}
					n_non_solid_neighbors++;
				}
				if (mac_grid.cellType(i, j - 1) != SOLID) {
					if (mac_grid.cellType(i, j - 1) == LIQUID) {
						_A.insert(_mem_pool->fluid_indices(i, j - 1), idx) =
							1 / pow(mac_grid.deltaX(), 2);
					}
					n_non_solid_neighbors++;
				}
				if (mac_grid.cellType(i, j + 1) != SOLID) {
					if (mac_grid.cellType(i, j + 1) == LIQUID) {
						_A.insert(_mem_pool->fluid_indices(i, j + 1), idx) =
							1 / pow(mac_grid.deltaX(), 2);
					}
					n_non_solid_neighbors++;
				}

				// Set diagonal value of matrix
				_A.insert(idx, idx) =
					- n_non_solid_neighbors / pow(mac_grid.deltaX(), 2);

				// Calculate divergence and store in b
				b[idx] = mac_grid.divVelX(i, j) + mac_grid.divVelY(i, j);
			}
		}
	}

	// Vector containing pressures for each cell
	VectorX x(n_fluid_cells);
	_cg_solver.compute(_A);
	x = _cg_solver.solve(b);

	// Loop throgh all cells to set new velocities from pressures
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			// Calculate indices
			int i_minus1 = CLAMP(i - 1, 0, mac_grid.sizeX() - 1);
			int j_minus1 = CLAMP(j - 1, 0, mac_grid.sizeY() - 1);
			int idx = _mem_pool->fluid_indices(i, j);
			int idx_i_minus1 = _mem_pool->fluid_indices(i_minus1, j);
			int idx_j_minus1 = _mem_pool->fluid_indices(i, j_minus1);
			
			if (idx >= 0 || idx_i_minus1 >= 0 || idx_j_minus1 >= 0)
			{
				// UGLY SOLUTION TO VOLUME LOSS HERE!!
				// MAKE PRESSURE PROPORTIONAL TO NUMBER OF PARTICLES OVER 9
				MyFloat k = 0.0;
                
                MyFloat particle_pressure_i_j =
                	k * (_mem_pool->n_particles(i, j) > 9 ?
                		_mem_pool->n_particles(i, j) - 8 : 0);
                MyFloat particle_pressure_i_minus1_j =
                	k * (_mem_pool->n_particles(i - 1, j) > 9 ?
                		_mem_pool->n_particles(i - 1, j) - 8 : 0);
                MyFloat particle_pressure_i_j_minus1 =
                	k * (_mem_pool->n_particles(i, j - 1) > 9 ?
                		_mem_pool->n_particles(i, j - 1) - 8 : 0);
                
				MyFloat p = idx >= 0 ? x(idx) + particle_pressure_i_j : 0;

				MyFloat p_i_minus1 = idx_i_minus1 >= 0 ?
					x(idx_i_minus1) + particle_pressure_i_minus1_j : 0;
				MyFloat p_j_minus1 = idx_j_minus1 >= 0 ?
					x(idx_j_minus1) + particle_pressure_i_j_minus1 : 0;

				// Get Pressure difference in x and y dimension
				MyFloat pressure_diff_x = p - p_i_minus1;
				MyFloat pressure_diff_y = p - p_j_minus1;

				// Current velocities
				MyFloat vel_x = mac_grid.velXHalfIndexed(i,j);
				MyFloat vel_y = mac_grid.velYHalfIndexed(i,j);

				// Calculate new velocity
				MyFloat new_vel_x =
					vel_x - dt / density * pressure_diff_x / mac_grid.deltaX();
				MyFloat new_vel_y =
					vel_y - dt / density * pressure_diff_y / mac_grid.deltaY();

				// Write data
				mac_grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);	
				mac_grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);
			}
		}
	}
	mac_grid.swapBuffers();
}

void FluidSolver::extendVelocityIndividual(MacGrid& mac_grid, int n_iterations)
{
	// First decide which velcities are valid. These are the ones adjacent to
	// liquid cells (indices will be different for x and y, hence the need of
	// two "valid mask grids")
    for (int j = 0; j < mac_grid.sizeY(); ++j)
    {
        for (int i = 0; i < mac_grid.sizeX(); ++i)
        {
        	if(	mac_grid.cellType(i, j) == LIQUID ||
        		mac_grid.cellType(i - 1, j) == LIQUID)
        	{
                (*_mem_pool->valid_mask_x_front_buffer)(i,j) = 1;
                (*_mem_pool->valid_mask_x_back_buffer)(i,j) = 1;
	            mac_grid.setVelXBackBufferHalfIndexed(
	            	i, j, mac_grid.velXHalfIndexed(i, j));
	            mac_grid.setVelXHalfIndexed(
	            	i, j, mac_grid.velXHalfIndexed(i, j));
        	}
        	else
        	{
        		(*_mem_pool->valid_mask_x_front_buffer)(i,j) = 0;
				(*_mem_pool->valid_mask_x_back_buffer)(i,j) = 0;
	            mac_grid.setVelXBackBufferHalfIndexed(i, j, 0);
	            mac_grid.setVelXHalfIndexed(i, j, 0);
        	}

        	if(	mac_grid.cellType(i, j) == LIQUID ||
        		mac_grid.cellType(i, j - 1) == LIQUID)
        	{
                (*_mem_pool->valid_mask_y_front_buffer)(i,j) = 1;
                (*_mem_pool->valid_mask_y_back_buffer)(i,j) = 1;
	            mac_grid.setVelYBackBufferHalfIndexed(i, j,
	            	mac_grid.velYHalfIndexed(i, j));
	            mac_grid.setVelYHalfIndexed(i, j,
	            	mac_grid.velYHalfIndexed(i, j));
        	}
        	else
        	{
        		(*_mem_pool->valid_mask_y_front_buffer)(i,j) = 0;
				(*_mem_pool->valid_mask_y_back_buffer)(i,j) = 0;
	            mac_grid.setVelYBackBufferHalfIndexed(i, j, 0);
	            mac_grid.setVelXHalfIndexed(i, j, 0);
        	}
        }
    }

    // Each iteration will extend the velocity outward one cell
    for (int iter = 0; iter < n_iterations; ++iter)
    {
        for (int j = 0; j < mac_grid.sizeY(); ++j)
        {
            for (int i = 0; i < mac_grid.sizeX(); ++i)
            {
            	// Check if this cell will be updated in the x dimension
            	if (_mem_pool->valid_mask_x_front_buffer->value(i, j) == 0 &&
            		mac_grid.cellType(i, j) != SOLID &&
            		mac_grid.cellType(i - 1, j) != SOLID)
            	{
	            	MyFloat new_vel_x = 0;
	            	int n_valid_neighbors_x = 0;

	            	// Get values of all valid neighbors
                    if(_mem_pool->valid_mask_x_front_buffer->value(i-1, j) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i-1, j);
                        n_valid_neighbors_x++;
                    }
                    if(_mem_pool->valid_mask_x_front_buffer->value(i, j-1) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i, j-1);
                        n_valid_neighbors_x++;
                    }
                    if(_mem_pool->valid_mask_x_front_buffer->value(i, j+1) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i, j+1);
                        n_valid_neighbors_x++;
                    }
                    if(_mem_pool->valid_mask_x_front_buffer->value(i+1, j) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBufferHalfIndexed(i+1, j);
                        n_valid_neighbors_x++;
                    }

	            	// Average the value for the current cell
	            	if (n_valid_neighbors_x > 0)
	            	{
	            		new_vel_x /= n_valid_neighbors_x;
	            		mac_grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);
						(*_mem_pool->valid_mask_x_back_buffer)(i, j) = 1;
	            	}
            	}

            	// Check if this cell will be updated in the y dimension
            	if (_mem_pool->valid_mask_y_front_buffer->value(i, j) == 0 &&
            		mac_grid.cellType(i, j) != SOLID &&
            		mac_grid.cellType(i, j - 1) != SOLID)
            	{
	            	MyFloat new_vel_y = 0;
	            	int n_valid_neighbors_y = 0;

	            	// Get values of all valid neighbors
                    if(_mem_pool->valid_mask_y_front_buffer->value(i-1, j) == 1)
                    {
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i-1, j);
                        n_valid_neighbors_y++;
                    }
                    if(_mem_pool->valid_mask_y_front_buffer->value(i, j-1) == 1)
                    {
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i, j-1);
                        n_valid_neighbors_y++;
                    }
                    if(_mem_pool->valid_mask_y_front_buffer->value(i, j+1) == 1)
                    {
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i, j+1);
                        n_valid_neighbors_y++;
                    }
                    if(_mem_pool->valid_mask_y_front_buffer->value(i+1, j) == 1)
                    {
                        new_vel_y += mac_grid.velYBackBufferHalfIndexed(i+1, j);
                        n_valid_neighbors_y++;
                    }

	            	// Average the value for the current cell
	            	if (n_valid_neighbors_y > 0)
	            	{
	            		new_vel_y /= n_valid_neighbors_y;
	            		mac_grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);
						(*_mem_pool->valid_mask_y_back_buffer)(i, j) = 1;
	            	}
            	}
            }
        }
        // Swap valid mask
        _mem_pool->swapValidMaskBuffer();
    }
    mac_grid.swapBuffers();
}


void FluidSolver::extendVelocityAvarageing(MacGrid& mac_grid, int n_iterations)
{
	// First decide which velcities are valid and write all velocities to the
	// Back buffers. Only one valid mask buffer is needed here since it is
	// done for each cell and not x and y velocities individually
	for (int j = 0; j < mac_grid.sizeY(); ++j)
    {
        for (int i = 0; i < mac_grid.sizeX(); ++i)
        {
        	if(mac_grid.cellType(i, j) == LIQUID)
        	{
                (*_mem_pool->valid_mask_x_front_buffer)(i,j) = 1;
                (*_mem_pool->valid_mask_x_back_buffer)(i,j) = 1;
        	}
        	else
        	{
        		(*_mem_pool->valid_mask_x_front_buffer)(i,j) = 0;
				(*_mem_pool->valid_mask_x_back_buffer)(i,j) = 0;
        	}
            mac_grid.setVelXBackBufferHalfIndexed(
            	i, j, mac_grid.velXHalfIndexed(i, j));
            mac_grid.setVelYBackBufferHalfIndexed(
            	i, j, mac_grid.velYHalfIndexed(i, j));
        }
    }

    for (int iter = 0; iter < n_iterations; ++iter)
    {
        for (int j = 0; j < mac_grid.sizeY(); ++j)
        {
            for (int i = 0; i < mac_grid.sizeX(); ++i)
            {
            	if (_mem_pool->valid_mask_x_front_buffer->value(i, j) == 0 &&
            		mac_grid.cellType(i, j) != SOLID)
            	{
	            	MyFloat new_vel_x = 0;
	            	MyFloat new_vel_y = 0;
	            	int n_valid_neighbors = 0;

	            	// Get values of all valid neighbors
                    if(_mem_pool->valid_mask_x_front_buffer->value(i-1, j) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBuffer(i-1, j);
                        new_vel_y += mac_grid.velYBackBuffer(i-1, j);
                        n_valid_neighbors++;
                    }
                    if(_mem_pool->valid_mask_x_front_buffer->value(i, j-1) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBuffer(i, j-1);
                        new_vel_y += mac_grid.velYBackBuffer(i, j-1);
                        n_valid_neighbors++;
                    }
                    if(_mem_pool->valid_mask_x_front_buffer->value(i, j+1) == 1)
                    {
                        new_vel_x += mac_grid.velXBackBuffer(i, j+1);
                        new_vel_y += mac_grid.velYBackBuffer(i, j+1);
                        n_valid_neighbors++;
                    }
                    if(_mem_pool->valid_mask_x_front_buffer->value(i+1, j) == 1)
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
						(*_mem_pool->valid_mask_x_back_buffer)(i, j) = 1;
	            	}
            	}
            }
        }
        // Swap valid mask
        _mem_pool->swapValidMaskBuffer();
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
			_mem_pool->vel_x_sum(i, j) = 0;
			_mem_pool->vel_y_sum(i, j) = 0;
			_mem_pool->weight_vel_x_sum(i, j) = 0;
			_mem_pool->weight_vel_y_sum(i, j) = 0;
		}
	}
	for (auto p = particle_set.begin(); p != particle_set.end(); p++)
	{
		_mem_pool->vel_x_sum.addToValueInterpolated(
			p->posX(), p->posY() - 0.5 * mac_grid.deltaY(), p->velX());
		_mem_pool->vel_y_sum.addToValueInterpolated(
			p->posX() - 0.5 * mac_grid.deltaX(), p->posY(), p->velY());
		_mem_pool->weight_vel_x_sum.addToValueInterpolated(
			p->posX(), p->posY() - 0.5 * mac_grid.deltaY(), 1.0);
		_mem_pool->weight_vel_y_sum.addToValueInterpolated(
			p->posX() - 0.5 * mac_grid.deltaX(), p->posY(), 1.0);
	}
	for (int j = 0; j < mac_grid.sizeY(); ++j)
	{
		for (int i = 0; i < mac_grid.sizeX(); ++i)
		{
			if (_mem_pool->weight_vel_x_sum(i, j) > 0.000001)
			{
				MyFloat new_vel_x =
					_mem_pool->vel_x_sum(i, j) /
					_mem_pool->weight_vel_x_sum(i, j);
				mac_grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);
			}
			if (_mem_pool->weight_vel_y_sum(i, j) > 0.000001)
			{
				MyFloat new_vel_y =
					_mem_pool->vel_y_sum(i, j) /
					_mem_pool->weight_vel_y_sum(i, j);
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
        MyFloat flip_vel_x = p->velX() +
        	mac_grid.velXDiffInterpolated(p->posX(), p->posY());
        MyFloat flip_vel_y = p->velY() +
        	mac_grid.velYDiffInterpolated(p->posX(), p->posY());

        p->setVelocity(
        	pic_vel_x * pic_ratio + flip_vel_x * (1 - pic_ratio),
        	pic_vel_y * pic_ratio + flip_vel_y * (1 - pic_ratio));
	}
}
