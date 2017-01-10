#ifndef MATH_DEFINITIONS
#define MATH_DEFINITIONS

#include <Eigen/SparseCore>
#include <math.h>

//#define USE_DOUBLE_PRECISION
#ifdef USE_DOUBLE_PRECISION
typedef double MyFloat;
typedef Eigen::VectorXd VectorX;
#else
typedef float MyFloat;
typedef Eigen::VectorXf VectorX;
#endif

inline MyFloat CLAMP(MyFloat d, MyFloat min, MyFloat max) {
  const MyFloat t = d < min ? min : d;
  return t > max ? max : t;
}

inline MyFloat smoothstep(MyFloat edge0, MyFloat edge1, MyFloat x)
{
    // Scale, bias and saturate x to 0..1 range
    x = CLAMP((x - edge0)/(edge1 - edge0), 0.0, 1.0); 
    // Evaluate polynomial
    return x*x*(3 - 2*x);
}

inline MyFloat gaussian(MyFloat x, MyFloat sigma, MyFloat mu)
{
	MyFloat a = 1 / (sigma * sqrt(2 * M_PI));
	MyFloat x_minus_b = x - mu;
	return a * exp(-(x_minus_b * x_minus_b) / (2 * sigma*sigma));
}

#endif
