#include <LevelSet.h>

LevelSet::LevelSet(int size_x, int size_y, MyFloat length_x, MyFloat length_y) :
	SizedGrid<MyFloat>(size_x, size_y, length_x / size_x, length_y / size_y),
	_LENGTH_X(length_x),
	_LENGTH_Y(length_y)
{
	for (int j = 0; j < _SIZE_Y; ++j)
	{
		for (int i = 0; i < _SIZE_X; ++i)
		{
			(*this)(i, j) = 0;
		}
	}
}

LevelSet::~LevelSet()
{

}

MyFloat LevelSet::distance(int from_i, int from_j, int to_i, int to_j)
{
	MyFloat i_diff = (to_i - from_i);
	MyFloat j_diff = (to_j - from_j);
	return sqrt(i_diff * i_diff + j_diff * j_diff);
}

MyFloat LevelSet::computeUpwindGradientX(int i, int j, MyFloat vel_x)
{
	int i_plus1 = i + 1;
	int i_minus1 = i - 1;
	// Border cases
	i = CLAMP(i, 0, _SIZE_X - 1);
	j = CLAMP(j, 0, _SIZE_Y - 1);
	i_plus1 = CLAMP(i_plus1, 0, _SIZE_X - 1);
	i_minus1 = CLAMP(i_minus1, 0, _SIZE_X - 1);
	
	MyFloat grad_x =
		vel_x < 0 ?
		(*this)(i_plus1, j) - (*this)(i, j) :
		(*this)(i, j) - (*this)(i_minus1, j);
	return grad_x;
}

MyFloat LevelSet::computeUpwindGradientY(int i, int j, MyFloat vel_y)
{
	int j_plus1 = j + 1;
	int j_minus1 = j - 1;
	// Border cases
	i = CLAMP(i, 0, _SIZE_X - 1);
	j = CLAMP(j, 0, _SIZE_Y - 1);
	j_plus1 = CLAMP(j_plus1, 0, _SIZE_Y - 1);
	j_minus1 = CLAMP(j_minus1, 0, _SIZE_Y - 1);
	
	MyFloat grad_y =
		vel_y < 0 ?
		(*this)(i, j_plus1) - (*this)(i, j) :
		(*this)(i, j) - (*this)(i, j_minus1);
	return grad_y;
}

MyFloat LevelSet::lengthX() const
{
	return _LENGTH_X;
}

MyFloat LevelSet::lengthY() const
{
	return _LENGTH_Y;
}