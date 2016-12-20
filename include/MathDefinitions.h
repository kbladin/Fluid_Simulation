#ifndef MATH_DEFINITIONS
#define MATH_DEFINITIONS

#define USE_DOUBLE_PRECISION

#ifdef USE_DOUBLE_PRECISION
typedef double MyFloat;
#else
typedef float MyFloat;
#endif

#endif