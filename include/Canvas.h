#ifndef CANVAS_H
#define CANVAS_H

#include <math.h>

class Color
{
public:
	Color();
	Color(float red, float green, float blue);
	~Color();
	float r, g, b;
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
	Color _line_color;
	Color _fill_color;
	Color* _pixel_data;
	const int _WIDTH;
	const int _HEIGHT;
};

#endif
