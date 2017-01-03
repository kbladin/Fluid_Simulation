#ifndef ODE_SOLVER
#define ODE_SOLVER

#include <functional>
#include "MathDefinitions.h"

template <class input_t, class output_t>
class OdeSolver
{
public:
	OdeSolver(){};
	~OdeSolver(){};
	
	virtual output_t step(
		std::function<output_t(input_t, output_t)> f,
		input_t x,
		output_t y,
		MyFloat h) = 0;
};

template <class input_t, class output_t>
class EulerExplicit
{
public:
	EulerExplicit(){};
	~EulerExplicit(){};
	
	virtual output_t step(std::function<output_t(
		input_t, output_t)> f,
		input_t x,
		output_t y,
		MyFloat h);
};

template <class input_t, class output_t>
output_t EulerExplicit<input_t, output_t>::step(
	std::function<output_t(input_t, output_t)> f,
	input_t x,
	output_t y,	
	MyFloat h)
{
	return f(x + h, y);
}

template <class input_t, class output_t>
class RK3
{
public:
	RK3(){};
	~RK3(){};
	
	virtual output_t step(
		std::function<output_t(input_t, output_t)> f,
		input_t x,
		output_t y,
		MyFloat h);
};

template <class input_t, class output_t>
output_t RK3<input_t, output_t>::step(
	std::function<output_t(input_t, output_t)> f,
	input_t x,
	output_t y,	
	MyFloat h)
{
	output_t k_1 = f(x, y);
	output_t k_2 = f(x + k_1 * 1.0/2 * h, y);
	output_t k_3 = f(x + k_2 * 3.0/4 * h, y);
	return (k_1 * 2 + k_2 * 3 + k_3 * 4) * 1.0/9 * h;
}

template <class input_t, class output_t>
class RK4
{
public:
	RK4(){};
	~RK4(){};
	
	virtual output_t step(
		std::function<output_t(input_t, output_t)> f,
		input_t x,
		output_t y,
		MyFloat h);
};

template <class input_t, class output_t>
output_t RK4<input_t, output_t>::step(
	std::function<output_t(input_t, output_t)> f,
	input_t x,
	output_t y,	
	MyFloat h)
{
	output_t f1 = f(x, y);
	output_t f2 = f(x + h/2.0, y + f1/2.0 * h);
	output_t f3 = f(x + h/2.0, y + f2/2.0 * h);
	output_t f4 = f(x + h, y + f3 * h);

	return (f1 + f2*2 + f3*2 + f4)/6 * h;
}

#endif
