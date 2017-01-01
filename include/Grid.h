#ifndef GRID_H
#define GRID_H

#include <cassert>

#include <iostream>
#include <vector>

#include "MathDefinitions.h"

class GridInterface
{
public:
	GridInterface(int size_x, int size_y, MyFloat delta_x = 1, MyFloat delta_y = 1) :
		_SIZE_X(size_x),
		_SIZE_Y(size_y),
		_DELTA_X(delta_x),
		_DELTA_Y(delta_y) { };
	~GridInterface() {};

	// Transform
	inline void linearTo2D(int idx, int* i, int* j) const
	{
		*i = idx % (_SIZE_X);
		*j = idx / (_SIZE_X);
	};
	inline int twoDToLinear(int i, int j) const
	{
	    assert(this->indexIsValid(i, j));
		return i + j * _SIZE_X;
	};
	inline void worldToCell(MyFloat x, MyFloat y, int* i, int* j) const
	{
		*i = x / _DELTA_X;
		*j = y / _DELTA_Y;
	};
	inline void cellToWorld(int i, int j, MyFloat* x, MyFloat* y) const
	{
	    assert(this->indexIsValid(i, j));   
		*x = i * _DELTA_X;
		*y = j * _DELTA_Y;
	};

	inline bool indexIsValid(int i, int j) const
	{
	    return i >= 0 && i < _SIZE_X && j >= 0 && j < _SIZE_Y;
	};

	// Get
	inline int sizeX() const { return _SIZE_X; };
	inline int sizeY() const { return _SIZE_Y; };
    inline MyFloat deltaX() const { return _DELTA_X; };
	inline MyFloat deltaY() const { return _DELTA_Y; };
    inline MyFloat lengthX() const { return _SIZE_X * _DELTA_X; };
	inline MyFloat lengthY() const { return _SIZE_Y * _DELTA_Y; };

protected:
	// Constants
	int _SIZE_X;
	int _SIZE_Y;
	MyFloat _DELTA_X;
	MyFloat _DELTA_Y;
};

template <class T>
class Grid : public GridInterface
{
public:
	Grid(int size_x, int size_y, MyFloat delta_x = 1, MyFloat delta_y = 1);
	~Grid();

	// Get
	inline T value(int i, int j) const;
    /** Value interpolated with world coordinates as input.
		Reads from the four closest grid points.
	*/
	inline T valueInterpolated(MyFloat x, MyFloat y) const;
	
	// Set
	inline T& operator()(int i, int j);
	/** Value interpolated with world coordinates as input.
		Writes to the four closest grid points.
	*/
	inline void addToValueInterpolated(MyFloat x, MyFloat y, T value);

protected:
	// Data
	std::vector<T> data;
};

// The functions are defined here due to the template class.

template <class T>
Grid<T>::Grid(int size_x, int size_y, MyFloat delta_x, MyFloat delta_y) :
	GridInterface(size_x, size_y, delta_x, delta_y)
{
	int n_elements = _SIZE_X * _SIZE_Y;
    data.resize(n_elements);
    for (int i = 0; i < n_elements; ++i)
	{
		data[i] = T();
	}
}

template <class T>
Grid<T>::~Grid()
{

}

template <class T>
inline T Grid<T>::value(int i, int j) const
{
	return data[twoDToLinear(i, j)];
}

template <class T>
inline T Grid<T>::valueInterpolated(MyFloat x, MyFloat y) const
{
	// Calculate indices
	int i = x / _DELTA_X;
	int j = y / _DELTA_Y;
	MyFloat i_frac = x / _DELTA_X - i;
	MyFloat j_frac = y / _DELTA_Y - j;
    
    i = CLAMP(i, 0, this->_SIZE_X - 1);
    j = CLAMP(j, 0, this->_SIZE_Y - 1);
    int i_plus1 = CLAMP(i + 1, 0, this->_SIZE_X - 1);
    int j_plus1 = CLAMP(j + 1, 0, this->_SIZE_Y - 1);
    
	// First interpolate in x, then in y
	MyFloat value_00 = this->value(i,j);
	MyFloat value_10 = this->value(i_plus1,j);
	MyFloat value_01 = this->value(i,j_plus1);
	MyFloat value_11 = this->value(i_plus1,j_plus1);

	// Interpolate x
	MyFloat value_0 = (1 - i_frac) * value_00 + i_frac * value_10;
	MyFloat value_1 = (1 - i_frac) * value_01 + i_frac * value_11; 

	// Interpolate y
	MyFloat value = (1 - j_frac) * value_0 + j_frac * value_1;
	return value;
}

template <class T>
inline T& Grid<T>::operator()(int i, int j)
{
	return data[twoDToLinear(i, j)];
}

template <class T>
inline void Grid<T>::addToValueInterpolated(MyFloat x, MyFloat y, T value)
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
