#ifndef MATH_DEFINITIONS
#define MATH_DEFINITIONS

#include <Eigen/SparseCore>

#define CLAMP(x, low, high) (x < low ? low : (x > high ? high : x))

#define USE_DOUBLE_PRECISION
#ifdef USE_DOUBLE_PRECISION
typedef double MyFloat;
typedef Eigen::VectorXd VectorX;
#else
typedef float MyFloat;
typedef Eigen::VectorXf VectorX;
#endif

#endif