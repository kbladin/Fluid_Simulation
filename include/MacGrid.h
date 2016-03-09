#ifndef MAC_GRID_H
#define MAC_GRID_H

#include <vector>
#include <iostream>
#include <random>

#include <assert.h>

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

enum CellType
{
	LIQUID, AIR, SOLID
};

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
	void worldToCell(double x, double y, int* i, int* j);
	void cellToWorld(int i, int j, double* x, double* y);
	void linearIndexTo2DIndexXBorders(int idx, int* i, int* j);
	void linearIndexTo2DIndexYBorders(int idx, int* i, int* j);
	void linearIndexTo2DIndexCellCenter(int idx, int* i, int* j);

	// Getters
	double velX(int x, int y) const;
	double velY(int x, int y) const;
	double velXHalfIndexed(int x, int y) const;
	double velYHalfIndexed(int x, int y) const;
	double color(int x, int y) const;
	double velXInterpolated(double x, double y) const;
	double velYInterpolated(double x, double y) const;
	double colorInterpolated(double x, double y) const;
	double divVelX(int x, int y) const;
	double divVelY(int x, int y) const;
	int sizeX() const;
	int sizeY() const;
	double lengthX() const;
	double lengthY() const;
	double deltaX() const;
	double deltaY() const;
	CellType cellType(int i, int j) const;
	CellType cellTypeXHalfIndexed(int x, int y) const;
	CellType cellTypeYHalfIndexed(int x, int y) const;
	
	// Setters
	void setVelX(int x, int y, double vel_x);
	void setVelY(int x, int y, double vel_y);
	void setColor(int x, int y, double color);
	void setCellType(int x, int y, CellType cell_type);
	void addToVelXInterpolated(double x, double y, double vel_x);
	void addToVelYInterpolated(double x, double y, double vel_y);
	void addToColorInterpolated(double x, double y, double color);

private:
	void _swapBuffers();

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
	std::vector<double> _vel_x_front_buffer;
	std::vector<double> _vel_y_front_buffer;
	std::vector<double> _color_front_buffer;
	
	std::vector<double> _vel_x_back_buffer;
	std::vector<double> _vel_y_back_buffer;
	std::vector<double> _color_back_buffer;

	std::vector<CellType> _cell_type_buffer;

	// Sparse matrix for CG solve
    Eigen::SparseMatrix<double> A;
};

#endif