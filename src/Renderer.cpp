#include <Renderer.h>

#define CLAMP(x, low, high) x < low ? low : (x > high ? high : x)

Renderer::Renderer(MacGrid* grid)
{
	_grid = grid;
	_canvas = new Canvas(200, 200);
}

Renderer::~Renderer()
{
	delete _canvas;
}

void Renderer::renderCanvas()
{
	_canvas->fill(Color(1,1,1));

	// Cell size in pixels
	double cell_size_x = double(_canvas->width()) / _grid->sizeX();
	double cell_size_y = double(_canvas->height()) / _grid->sizeY();
	double line_scale = 20;
	// Render all velocities as lines
	for (int j = 0; j < _grid->sizeY(); ++j)
	{
		for (int i = 0; i < _grid->sizeX(); ++i)
		{
			double vel_x = _grid->velX(i, j);
			double vel_y = _grid->velY(i, j);
			double vel_norm = sqrt(vel_x * vel_x + vel_y * vel_y);
			// Since we draw from the center of each cell, add 0.5
			int from_x = (0.5 + i) * cell_size_x;
			int from_y = (0.5 + j) * cell_size_y;
			int to_x = from_x + vel_x * line_scale;
			int to_y = from_y + vel_y * line_scale;

			_canvas->setLineColor(Color(vel_x / vel_norm, vel_y / vel_norm,0));

			_canvas->drawLine(from_x, from_y, to_x, to_y);
		}
	}
}

void Renderer::writeToPpm( const char* file_path, float min_val, float max_val)
{
	int w = _grid->sizeX();
	int h = _grid->sizeY();

	unsigned char* byte_data = new unsigned char[w * h * 3];
	double value;
	unsigned char byte_value;

	// Copy and cast data
	for (int j = 0; j < h; ++j)
	{
		for (int i = 0; i < w; ++i)
		{
			value = _grid->color(i, j);
			byte_value =
				static_cast<unsigned char>
				CLAMP((value - min_val) / (max_val - min_val), 0, 1) * 255;

			byte_data[(i + j * h) * 3 + 0] = byte_value; // Red
			byte_data[(i + j * h) * 3 + 1] = byte_value; // Green
			byte_data[(i + j * h) * 3 + 2] = byte_value; // Blue
		}
	}

	FILE *fp = fopen(file_path, "wb"); // b - binary mode
	fprintf(fp, "P6\n%d %d\n255\n", w, h);
	fwrite(byte_data, 1, w * h * 3, fp);
	fclose(fp);
	delete[] byte_data;
}

void Renderer::writeCanvasToPpm( const char* file_path)
{
	int w = _canvas->width();
	int h = _canvas->height();

	unsigned char* byte_data = new unsigned char[w * h * 3];
	unsigned char byte_value_r;
	unsigned char byte_value_g;
	unsigned char byte_value_b;

	// Copy and cast data
	for (int j = 0; j < h; ++j)
	{
		for (int i = 0; i < w; ++i)
		{ 
			Color c = _canvas->pixel(i,j);
			byte_value_r = static_cast<unsigned char>(CLAMP(c.r, 0, 1) * 255);
			byte_value_g = static_cast<unsigned char>(CLAMP(c.g, 0, 1) * 255);
			byte_value_b = static_cast<unsigned char>(CLAMP(c.b, 0, 1) * 255);

			byte_data[(i + j * w) * 3 + 0] = byte_value_r; // Red
			byte_data[(i + j * w) * 3 + 1] = byte_value_g; // Green
			byte_data[(i + j * w) * 3 + 2] = byte_value_b; // Blue
		}
	}

	FILE *fp = fopen(file_path, "wb"); // b - binary mode
	fprintf(fp, "P6\n%d %d\n255\n", w, h);
	fwrite(byte_data, 1, w * h * 3, fp);
	fclose(fp);
	delete[] byte_data;
}