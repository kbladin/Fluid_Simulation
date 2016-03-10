#ifndef GRID_H
#define GRID_H

#include <vector>

template <class T>
class Grid
{
public:
	// Constructors / Destructor
	Grid(int size_x, int size_y);
	~Grid();
	Grid(Grid&& rhv);
	Grid<T>& operator=(Grid<T> to_copy);
	
	// Transform
	void linearTo2D(int idx, int* i, int* j) const;
	int twoDToLinear(int i, int j) const;

	// Get
	T value(int i, int j) const;

	// Set
	T& operator()(int i, int j);

private:
	// Constants
	const int _SIZE_X;
	const int _SIZE_Y;

	// Data
	std::vector<T> data;
};

// The functions are defined here due to the template class.

#define CLAMP(x, low, high) x < low ? low : (x > high ? high : x)

template <class T>
Grid<T>::Grid(int size_x, int size_y) : _SIZE_X(size_x), _SIZE_Y(size_y)
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
Grid<T>::Grid(Grid&& rhv) : _SIZE_X(rhv._SIZE_X), _SIZE_Y(rhv._SIZE_Y)
{
	data = std::move(rhv.data);
}

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

/**
	Does not consider border cases for efficiency. That needs to be handled
	outside this function.
*/
template <class T>
T& Grid<T>::operator()(int i, int j)
{
	return data[twoDToLinear(i, j)];
}


#endif
