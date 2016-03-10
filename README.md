# Eulerian_Fluid_Simulation

## Implemented
* Staggered MAC Grid
* Bilinear interpolation for velocity advection
* Runge Kutta solver
* Pressure solve (using the linear algebra library "Eigen")
* Enforce Dirichlet boundary condition

Example simulation:

![](images/fluid.gif "Fluid Simulation")

## Not Yet Implemented
* Adaptive time steps
* Other interpolation alternatives
* Level set to separate fluid from air
* Velocity extension for level set
* 3D simulation
* Sophisticated rendering
