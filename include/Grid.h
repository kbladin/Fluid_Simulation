#ifndef GRID_H
#define GRID_H

#include <cassert>

#include <iostream>
#include <vector>

#include "MathDefinitions.h"

template <class T>
class Grid
{
public:
	Grid(int size_x, int size_y);
	~Grid();
	
	// Transform
	inline void linearTo2D(int idx, int* i, int* j) const;
	inline int twoDToLinear(int i, int j) const;

	// Get
	inline T value(int i, int j) const;
	inline int sizeX() const;
	inline int sizeY() const;
    inline bool indexIsValid(int i, int j);

	// Set
	inline T& operator()(int i, int j);

protected:
	// Constants
	int _SIZE_X;
	int _SIZE_Y;

	// Data
	std::vector<T> data;
};

template <class T>
class SizedGrid : public Grid<T>
{
public:
	// Constructors / Destructor
	SizedGrid(int size_x, int size_y, MyFloat delta_x, MyFloat delta_y);
	~SizedGrid();
	
	// Get
	/** Value interpolated with world coordinates as input.
		Reads from the four closest grid points.
	*/
	inline T valueInterpolated(MyFloat x, MyFloat y) const;
	inline MyFloat deltaX() const;
	inline MyFloat deltaY() const;
	void worldToCell(MyFloat x, MyFloat y, int* i, int* j) const;
	void cellToWorld(int i, int j, MyFloat* x, MyFloat* y) const;

	// Set
	/** Value interpolated with world coordinates as input.
		Writes to the four closest grid points.
	*/
	inline void addToValueInterpolated(MyFloat x, MyFloat y, T value);
private:
	MyFloat _DELTA_X;
	MyFloat _DELTA_Y;
};

// The functions are defined here due to the template class.

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

template <class T>
inline void Grid<T>::linearTo2D(int idx, int* i, int* j) const
{
	*i = idx % (_SIZE_X);
	*j = idx / (_SIZE_X);
}

template <class T>
inline int Grid<T>::twoDToLinear(int i, int j) const
{
    assert(indexIsValid(i, j));
	return i + j * _SIZE_X;
}

template <class T>
inline T Grid<T>::value(int i, int j) const
{
	return data[twoDToLinear(i, j)];
}

template <class T>
inline int Grid<T>::sizeX() const
{
	return _SIZE_X;
}

template <class T>
inline int Grid<T>::sizeY() const
{
	return _SIZE_Y;
}

template <class T>
inline bool Grid<T>::indexIsValid(int i, int j)
{
    return i >= 0 && i < _SIZE_X && j >= 0 && j < _SIZE_Y;
}

template <class T>
inline T& Grid<T>::operator()(int i, int j)
{
	return data[twoDToLinear(i, j)];
}

template <class T>
SizedGrid<T>::SizedGrid(int size_x, int size_y, MyFloat delta_x, MyFloat delta_y) :
	Grid<T>(size_x, size_y),
	_DELTA_X(delta_x),
	_DELTA_Y(delta_y)
{
}

template <class T>
SizedGrid<T>::~SizedGrid()
{

}

template <class T>
inline T SizedGrid<T>::valueInterpolated(MyFloat x, MyFloat y) const
{
	// Calculate indices
	int i = x / _DELTA_X;
	int j = y / _DELTA_Y;
	MyFloat i_frac = x / _DELTA_X - i;
	MyFloat j_frac = y / _DELTA_Y - j;

    assert(indexIsValid(i, j));
    
	// First interpolate in x, then in y
	MyFloat value_00 = this->value(i,j);
	MyFloat value_10 = this->value(i + 1,j);
	MyFloat value_01 = this->value(i,j + 1);
	MyFloat value_11 = this->value(i + 1,j + 1);

	// Interpolate x
	MyFloat value_0 = (1 - i_frac) * value_00 + i_frac * value_10;
	MyFloat value_1 = (1 - i_frac) * value_01 + i_frac * value_11; 

	// Interpolate y
	MyFloat value = (1 - j_frac) * value_0 + j_frac * value_1;
	return value;
}

template <class T>
inline MyFloat SizedGrid<T>::deltaX() const
{
	return _DELTA_X;
}

template <class T>
inline MyFloat SizedGrid<T>::deltaY() const
{
	return _DELTA_Y;
}

template <class T>
void SizedGrid<T>::worldToCell(MyFloat x, MyFloat y, int* i, int* j) const
{
	*i = x / _DELTA_X;
	*j = y / _DELTA_Y;
}

template <class T>
void SizedGrid<T>::cellToWorld(int i, int j, MyFloat* x, MyFloat* y) const
{
    assert(indexIsValid(i, j));
    
	*x = i * _DELTA_X;
	*y = j * _DELTA_Y;
}

template <class T>
inline void SizedGrid<T>::addToValueInterpolated(MyFloat x, MyFloat y, T value)
{
	// Calculate indices
	int i = x / _DELTA_X;
	int j = y / _DELTA_Y;
	int i_plus1 = i + 1;
	int j_plus1 = j + 1;
	MyFloat i_frac = x / _DELTA_X - i;
	MyFloat j_frac = y / _DELTA_Y - j;

	// Border cases
	i = CLAMP(i, 0, this->_SIZE_X - 1);
	j = CLAMP(j, 0, this->_SIZE_Y - 1);
	i_plus1 = CLAMP(i_plus1, 0, this->_SIZE_X - 1);
	j_plus1 = CLAMP(j_plus1, 0, this->_SIZE_Y - 1);
	
	// Spread in y
	MyFloat value_0 = (1 - j_frac) * value;
	MyFloat value_1 = j_frac * value;

	// Spread in x
	MyFloat value_00 = (1 - i_frac) * value_0;
	MyFloat value_10 = i_frac * value_0;
	MyFloat value_01 = (1 - i_frac) * value_1;
	MyFloat value_11 = i_frac * value_1;

	// Write data
	(*this)(i, j) 				+= value_00;
	(*this)(i_plus1, j) 		+= value_10;
	(*this)(i, j_plus1) 		+= value_01;
	(*this)(i_plus1, j_plus1) 	+= value_11;
}

#endif
