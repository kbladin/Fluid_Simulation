#ifndef GRID_H
#define GRID_H

#include <cassert>

#include <iostream>
#include <vector>

template <class T>
class Grid
{
public:
	// Constructors / Destructor
	Grid(int size_x, int size_y);
	~Grid();
	//Grid(Grid&& rhv);
	Grid<T>& operator=(Grid<T> to_copy);
	
	// Transform
	void linearTo2D(int idx, int* i, int* j) const;
	int twoDToLinear(int i, int j) const;

	// Get
	T value(int i, int j) const;
	int sizeX() const;
	int sizeY() const;

	// Set
	T& operator()(int i, int j);

protected:
	// Constants
	const int _SIZE_X;
	const int _SIZE_Y;

	// Data
	std::vector<T> data;
};

template <class T>
class SizedGrid : public Grid<T>
{
public:
	// Constructors / Destructor
	SizedGrid(int size_x, int size_y, double delta_x, double delta_y);
	~SizedGrid();
	//SizedGrid(SizedGrid&& rhv);	
	SizedGrid<T>& operator=(SizedGrid<T> to_copy);

	// Get
	T valueInterpolated(double x, double y) const;
	double deltaX() const;
	double deltaY() const;

	// Set
	void addToValueInterpolated(double x, double y, T value);
private:
	const double _DELTA_X;
	const double _DELTA_Y;
};

// The functions are defined here due to the template class.

#define CLAMP(x, low, high) x < low ? low : (x > high ? high : x)

template <class T>
Grid<T>::Grid(int size_x, int size_y) :
	_SIZE_X(size_x),
	_SIZE_Y(size_y)
{
	data.resize(_SIZE_X * _SIZE_Y);
	for (int i = 0; i < data.size(); ++i)
	{
		data[i] = T(); // Call default constructor
	}
}

template <class T>
Grid<T>::~Grid()
{

} 
/*
template <class T>
Grid<T>::Grid(Grid&& rhv) :
	_SIZE_X(rhv._SIZE_X),
	_SIZE_Y(rhv._SIZE_Y)
{
	data = std::move(rhv.data);
}
*/
template <class T>
Grid<T>& Grid<T>::operator=(Grid<T> to_copy)
{
	data = to_copy.data;
	return *this;
}

template <class T>
void Grid<T>::linearTo2D(int idx, int* i, int* j) const
{
	*i = idx % (_SIZE_X);
	*j = idx / (_SIZE_X);
}

template <class T>
int Grid<T>::twoDToLinear(int i, int j) const
{
	return i + j * _SIZE_X;
}

/**
	Does not consider border cases for efficiency. That needs to be handled
	outside this function.
*/
template <class T>
T Grid<T>::value(int i, int j) const
{
	return data[twoDToLinear(i, j)];
}

template <class T>
int Grid<T>::sizeX() const
{
	return _SIZE_X;
}

template <class T>
int Grid<T>::sizeY() const
{
	return _SIZE_Y;
}

/**
	Does not consider border cases for efficiency. That needs to be handled
	outside this function.
*/
template <class T>
T& Grid<T>::operator()(int i, int j)
{
	return data[twoDToLinear(i, j)];
}

template <class T>
SizedGrid<T>::SizedGrid(int size_x, int size_y, double delta_x, double delta_y) :
	Grid<T>(size_x, size_y),
	_DELTA_X(delta_x),
	_DELTA_Y(delta_y)
{
}

template <class T>
SizedGrid<T>::~SizedGrid()
{

} 
/*
template <class T>
SizedGrid<T>::SizedGrid(SizedGrid&& rhv) :
	Grid<T>(std::move(rhv)),
	_DELTA_X(rhv._DELTA_X),
	_DELTA_Y(rhv._DELTA_Y)
{
}
*/
template <class T>
SizedGrid<T>& SizedGrid<T>::operator=(SizedGrid<T> to_copy)
{
	this->data = to_copy.data;
	return *this;
}

template <class T>
T SizedGrid<T>::valueInterpolated(double x, double y) const
{
		// Calculate indices
	int i = x / this->_DELTA_X;
	int j = y / this->_DELTA_Y; // -0.5 Due to the MAC grid structure
	int i_plus1 = i + 1;
	int j_plus1 = j + 1;
	double i_frac = x / this->_DELTA_X - i;
	double j_frac = y / this->_DELTA_Y - j;

	assert(i_frac <= 1);
	assert(j_frac <= 1);

	// Border cases
	i = CLAMP(i, 0, this->_SIZE_X - 1);
	j = CLAMP(j, 0, this->_SIZE_Y - 1);
	i_plus1 = CLAMP(i_plus1, 0, this->_SIZE_X - 1);
	j_plus1 = CLAMP(j_plus1, 0, this->_SIZE_Y - 1);
	
	// First interpolate in x, then in y
	double value_00 = this->value(i,j);
	double value_10 = this->value(i_plus1,j);
	double value_01 = this->value(i,j_plus1);
	double value_11 = this->value(i_plus1,j_plus1);

	// Interpolate x
	double value_0 = (1 - i_frac) * value_00 + i_frac * value_10;
	double value_1 = (1 - i_frac) * value_01 + i_frac * value_11; 

	// Interpolate y
	double value = (1 - j_frac) * value_0 + j_frac * value_1;
	return value;
}

template <class T>
double SizedGrid<T>::deltaX() const
{
	return _DELTA_X;
}

template <class T>
double SizedGrid<T>::deltaY() const
{
	return _DELTA_Y;
}

/**
	Writes to the four closest grid points, writes to back buffer.
*/
template <class T>
void SizedGrid<T>::addToValueInterpolated(double x, double y, T value)
{
	// Calculate indices
	int i = x / _DELTA_X; // -0.5 Due to the MAC grid structure
	int j = y / _DELTA_Y;
	int i_plus1 = i + 1;
	int j_plus1 = j + 1;
	double i_frac = x / _DELTA_X - i;
	double j_frac = y / _DELTA_Y - j;

	assert(i_frac <= 1);
	assert(j_frac <= 1);

	// Border cases
	i = CLAMP(i, 0, this->_SIZE_X - 1);
	j = CLAMP(j, 0, this->_SIZE_Y - 1);
	i_plus1 = CLAMP(i_plus1, 0, this->_SIZE_X - 1);
	j_plus1 = CLAMP(j_plus1, 0, this->_SIZE_Y - 1);
	
	// Spread in y
	double value_0 = (1 - j_frac) * value;
	double value_1 = j_frac * value;

	// Spread in x
	double value_00 = (1 - i_frac) * value_0;
	double value_10 = i_frac * value_0;
	double value_01 = (1 - i_frac) * value_1;
	double value_11 = i_frac * value_1;

	// Write data
	(*this)(i, j) 				+= value_00;
	(*this)(i_plus1, j) 		+= value_10;
	(*this)(i, j_plus1) 		+= value_01;
	(*this)(i_plus1, j_plus1) 	+= value_11;
}

#endif
