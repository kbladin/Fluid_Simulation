# Fluid_Simulation

## Implemented
* Staggered MAC Grid
* Different advection schemes
	* Semi Lagrangian advection with
		* Bilinear interpolation
		* Runge Kutta and Euler solver
	* PIC advection with bilinear interpolation
	* FLIP advection with bilinear interpolation
	* PIC / FLIP combination
* Pressure solve using the conjugate gradient method with the linear algebra library "Eigen"
* Enforce Dirichlet boundary condition
* Marker particles for air/fluid separation
* Velocity extension using brute force search through grid

Examples:

![](images/PIC.gif "PIC simulation") ![](images/FLIP.gif "FLIP simulation") ![](images/PICFLIP98.gif "2% PIC 98% FLIP simulation") ![](images/squirt1.gif "Volume preservation") ![](images/squirt2.gif "Volume preservation and high viscisity") ![](images/squirt4_k0.gif "No volume preservation. Smaller delta time")

Youtube video of real time simulation:

[![Real time fluid simulation](https://img.youtube.com/vi/gkGh4Xw5LDQ/0.jpg)](https://www.youtube.com/watch?v=gkGh4Xw5LDQ)

## Not Yet Implemented
* Adaptive time steps
* Vorticity confinement
* Make sure marker particles are outside of obstacles
* Other interpolation alternatives (quadratic, cubic, Catmull Rom)
* Other air/fluid separation techniques (for example level set)
* 3D simulation
