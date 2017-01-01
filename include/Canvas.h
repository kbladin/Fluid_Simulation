#ifndef CANVAS_H
#define CANVAS_H

#include <Grid.h>

#include <math.h>

#include "MathDefinitions.h"

class Color
{
public:
	Color();
	Color(MyFloat red, MyFloat green, MyFloat blue);
	~Color();
	MyFloat r, g, b;

	inline Color& operator+=(Color& c)
	{
		r += c.r;
		g += c.g;
		b += c.b;
		return *this;
	};
	inline Color operator*(MyFloat& s)
	{
		Color c = *this;
		c.r *= s;
		c.g *= s;
		c.b *= s;
		return c;
	};
private:
};

class Canvas
{
public:
	Canvas(int width, int height);
	~Canvas();
	void drawLine(int from_x, int from_y, int to_x, int to_y);
	void drawPoint(int pos_x, int pos_y, int size);
	void fillRectangle(int min_x, int min_y, int max_x, int max_y);
	void setPixel(int i, int j, Color c);
	void addToPixel(int i, int j, Color c);
	void fill(Color c);

	// Getters
	Color lineColor();
	Color fillColor();
	int width();
	int height();
	Color pixel(int i, int j);
	// Setters
	void setLineColor(Color c);
	void setFillColor(Color c);

private:
    const int _WIDTH;
    const int _HEIGHT;
    Color _line_color;
	Color _fill_color;
	Grid<Color> _pixel_data;
};

#endif
