#ifndef MAC_GRID_H
#define MAC_GRID_H

#include <Grid.h>

#include <vector>
#include <iostream>
#include <random>

#include <assert.h>

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <glm/glm.hpp>

enum CellType
{
	LIQUID, AIR, SOLID
};

template class Grid<double>;
template class Grid<CellType>;

//! Let this have dense arrays for now
class MacGrid
{
public:
	MacGrid(int size_x, int size_y, double length_x, double length_y);
	~MacGrid();

	// Simple semi Lagrangian advection with Euler integration
	void advect(double dt);
	void addExternalForce(double dt, double F_x, double F_y);
	void pressureSolve(double dt);
	void enforceDirichlet();

	void clearCellTypeBuffer();

	// Index transforms
	void worldToCell(double x, double y, int* i, int* j) const;
	void cellToWorld(int i, int j, double* x, double* y) const;
	void linearTo2DCellCenter(int idx, int* i, int* j) const;
	int twoDToLinearCellCenter(int i, int j) const;

	// Getters
	double velX(int i, int j) const;
	double velY(int i, int j) const;
	double velXHalfIndexed(int i, int j) const;
	double velYHalfIndexed(int i, int j) const;
	double color(int i, int j) const;
	double velXInterpolated(double x, double y) const;
	double velYInterpolated(double x, double y) const;
	double colorInterpolated(double x, double y) const;
	double divVelX(int i, int j) const;
	double divVelY(int i, int j) const;
	glm::dmat2 computeVelocityGradientMatrix(int i, int j);
	int sizeX() const;
	int sizeY() const;
	double lengthX() const;
	double lengthY() const;
	double deltaX() const;
	double deltaY() const;
	CellType cellType(int i, int j) const;
	CellType cellTypeXHalfIndexed(int i, int j) const;
	CellType cellTypeYHalfIndexed(int i, int j) const;
	
	// Setters
	void setVelX(int i, int j, double vel_x);
	void setVelY(int i, int j, double vel_y);
	void setColor(int i, int j, double color);
	void setCellType(int i, int j, CellType cell_type);
	void addToVelXInterpolated(double x, double y, double vel_x);
	void addToVelYInterpolated(double x, double y, double vel_y);
	void addToColorInterpolated(double x, double y, double color);

	void _swapBuffers();

private:

	void _getAdvectedPositionRK3(
		double x_pos,
		double y_pos,
		double dt,
		double* x,
		double* y);
	void _getAdvectedPositionForwardEuler(
		double x_pos,
		double y_pos,
		double dt,
		double* x,
		double* y);
	void _advectVelX(double dt);
	void _advectVelY(double dt);
	void _advectColor(double dt);

	// Constants
	const int _SIZE_X;
	const int _SIZE_Y;
	const double _LENGTH_X;
	const double _LENGTH_Y;
	const double _DELTA_X;
	const double _DELTA_Y;

	// Always render to back bufer from front buffer, then swap them
	// Since advection can not be done in place, another set of data is needed
	// (except for when adding forces)
	SizedGrid<double> _vel_x_front_buffer;
	SizedGrid<double> _vel_y_front_buffer;
	SizedGrid<double> _color_front_buffer;
	
	SizedGrid<double> _vel_x_back_buffer;
	SizedGrid<double> _vel_y_back_buffer;
	SizedGrid<double> _color_back_buffer;

	Grid<CellType> _cell_type_buffer;

	// Sparse matrix for CG solve
    Eigen::SparseMatrix<double> A;
};

#endif