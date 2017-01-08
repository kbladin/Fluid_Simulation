#ifndef ODE_SOLVER
#define ODE_SOLVER

#include <functional>

/**
	The OdeSolver has one function to iteratively solve an ordinary differential
	equation on the form dy/dx = f(x, y) for any function f. The function step
	returns dy given the function f, initial conditions (current state of function)
	parameters x and y as well as a step length h also known as dx or delta x.

	The function can be in any form and its input and output can have arbitrary
	types as long as the input type, output type and step type support the
	+, * and / operation between each other.

	Example implementation using a RK3 solver:

	<code>
	RK3<double, double> solver; // input and output types are double

    std::function<double(double, double)> f =
    	[](double x, double y)
    {
        return 3 * std::exp(-4*x) - 2*y;
    };
    
	double x = 0; // Start value x
    double y = 1; // Start value y
    double dx = 0.1; // Step length
    double end = 4;

    for (x = 0; x < end;)
    {
    	double dy = solver.step(f, x, y, dx);
    	y += dy;
        x += dx;

    	printf("%10f, %10f \n", x, y);
    }
    </code>
*/
template <class input_t, class output_t, class step_t>
class OdeSolver
{
public:
	OdeSolver(){};
	~OdeSolver(){};
	
	/**
		\param f is a function describing the derivative of the function
		f(x, y) = dy/dx.
		\param x is the current state of the input.
		\param y is the current state of the outpur.
		\returns A differential in ouput type so that the equation can be
		iteratively solved.
	*/
	virtual output_t step(
		std::function<output_t(input_t, output_t)> f,
		input_t x,
		output_t y,
		step_t h) = 0;
};

template <class input_t, class output_t, class step_t>
class EulerExplicit
{
public:
	EulerExplicit(){};
	~EulerExplicit(){};
	
	virtual output_t step(std::function<output_t(
		input_t, output_t)> f,
		input_t x,
		output_t y,
		step_t h);
};

template <class input_t, class output_t, class step_t>
output_t EulerExplicit<input_t, output_t, step_t>::step(
	std::function<output_t(input_t, output_t)> f,
	input_t x,
	output_t y,	
	step_t h)
{
	return f(x + h, y) * h;
}

template <class input_t, class output_t, class step_t>
class RK3
{
public:
	RK3(){};
	~RK3(){};
	
	virtual output_t step(
		std::function<output_t(input_t, output_t)> f,
		input_t x,
		output_t y,
		step_t h);
};

template <class input_t, class output_t, class step_t>
output_t RK3<input_t, output_t, step_t>::step(
	std::function<output_t(input_t, output_t)> f,
	input_t x,
	output_t y,	
	step_t h)
{
	output_t k_1 = f(x, y);
	output_t k_2 = f(x + k_1 * h * 1.0/2, y);
	output_t k_3 = f(x + k_2 * h * 3.0/4, y);
	return (k_1 * 2 + k_2 * 3 + k_3 * 4) * h * 1.0/9;
}

template <class input_t, class output_t, class step_t>
class RK4
{
public:
	RK4(){};
	~RK4(){};
	
	virtual output_t step(
		std::function<output_t(input_t, output_t)> f,
		input_t x,
		output_t y,
		step_t h);
};

template <class input_t, class output_t, class step_t>
output_t RK4<input_t, output_t, step_t>::step(
	std::function<output_t(input_t, output_t)> f,
	input_t x,
	output_t y,	
	step_t h)
{
	output_t f1 = f(x, y);
	output_t f2 = f(x + h/2.0, y + h * f1/2.0);
	output_t f3 = f(x + h/2.0, y + h * f2/2.0);
	output_t f4 = f(x + h, y + h * f3);

	return h * (f1 + f2*2 + f3*2 + f4)/6;
}

#endif
