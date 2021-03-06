#include <Canvas.h>

#include "MathDefinitions.h"

Color::Color()
{
	r = 1;
	g = 1;
	b = 1;
}

Color::Color(MyFloat red, MyFloat green, MyFloat blue)
{
	r = red;
	g = green;
	b = blue;
}

Color::~Color()
{

}

Canvas::Canvas(int width, int height) :
	_WIDTH(width),
	_HEIGHT(height),
	_line_color(0,0,0),
	_fill_color(1,1,1),
	_pixel_data(width, height)
{

}

Canvas::~Canvas()
{
	
}

void Canvas::drawLine(int from_x, int from_y, int to_x, int to_y)
{
	MyFloat diff_x = (to_x - from_x);
	MyFloat diff_y = (to_y - from_y);
	
	MyFloat distance = sqrt(diff_x * diff_x + diff_y * diff_y);
	// Step one pixel at the time
	for (int i = 0; i < distance; ++i)
	{
		int x_idx = from_x + (diff_x / distance) * i;
		int y_idx = from_y + (diff_y / distance) * i;

		if (x_idx < 0 || x_idx >= _WIDTH ||
			y_idx < 0 || y_idx >= _HEIGHT)
		{
			break;
		}
		else
		{
			_pixel_data(x_idx, y_idx) = _line_color;
		}
	}
}

void Canvas::drawPoint(int pos_x, int pos_y, int size)
{
	int min_x = pos_x - size / 2;
	int max_x = pos_x + size / 2;
	int min_y = pos_y - size / 2;
	int max_y = pos_y + size / 2;
	
	fillRectangle({min_x, max_x, min_y, max_y});
}

void Canvas::fillRectangle(BBox<int> rect)
{
    rect.x_min = CLAMP(rect.x_min, 0, _WIDTH - 1);
    rect.x_max = CLAMP(rect.x_max, 0, _WIDTH - 1);
    rect.y_min = CLAMP(rect.y_min, 0, _HEIGHT - 1);
    rect.y_max = CLAMP(rect.y_max, 0, _HEIGHT - 1);
    for (int j = rect.y_min; j <= rect.y_max; ++j)
	{
		for (int i = rect.x_min; i <= rect.x_max; ++i)
		{
            if (i == rect.x_min || i == rect.x_max ||
                j == rect.y_min || j == rect.y_max)
            { // On border, use line color
                _pixel_data(i, j) = _line_color;
            }
            else
            { // Inside, use fill color
                _pixel_data(i, j) = _fill_color;
            }
		}
	}
}

void Canvas::setPixel(int i, int j, Color c)
{
	if (i >= 0 && i < _WIDTH &&
		j >= 0 && j < _HEIGHT)
	{
		_pixel_data(i, j) = c;
	}
}

void Canvas::addToPixel(int i, int j, Color c)
{
	if (i >= 0 && i < _WIDTH &&
		j >= 0 && j < _HEIGHT)
	{
		_pixel_data(i, j) += c;
	}
}

void Canvas::fill(Color c)
{
	for (int j = 0; j < _HEIGHT; ++j)
	{
		for (int i = 0; i < _WIDTH; ++i)
		{
			_pixel_data(i, j) = c;
		}
	}
}

// Getters
Color Canvas::lineColor()
{
	return _line_color;
}

Color Canvas::fillColor()
{
	return _fill_color;
}

int Canvas::width()
{
	return _WIDTH;
}

int Canvas::height()
{
	return _HEIGHT;
}

Color Canvas::pixel(int i, int j)
{
	return _pixel_data(i, j);
}

// Setters
void Canvas::setLineColor(Color c)
{
	_line_color = c;
}

void Canvas::setFillColor(Color c)
{
	_fill_color = c;
}
