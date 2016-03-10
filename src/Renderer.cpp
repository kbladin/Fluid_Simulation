#include <Renderer.h>

#define CLAMP(x, low, high) x < low ? low : (x > high ? high : x)

Renderer::Renderer(double x_min, double y_min, double x_max, double y_max)
{
	_x_min = x_min;
	_y_min = y_min;
	_x_max = x_max;
	_y_max = y_max;

	_canvas = std::unique_ptr<Canvas>(new Canvas(400, 400));
}

Renderer::~Renderer()
{
	
}

void Renderer::clearCanvas()
{
	_canvas->fill(Color(1,1,1));
}

void Renderer::renderGridCellsToCanvas(const MacGrid* grid)
{
	int grid_size_x = grid->sizeX();
	int grid_size_y = grid->sizeY();

	double scale_x = _canvas->width() / (_x_max - _x_min);
	double scale_y = _canvas->height() / (_y_max - _y_min);

	double translate_x = 0.5 * (_x_min * _canvas->width());
	double translate_y = 0.5 * (_y_min * _canvas->height());

	// Cell size in pixels
	double cell_size_x = grid->deltaX() * scale_x;
	double cell_size_y = grid->deltaY() * scale_y;

	_canvas->setLineColor(Color(1,1,1));

	// Render all velocities as lines
	for (int j = 0; j < grid_size_y; ++j)
	{
		for (int i = 0; i < grid_size_x; ++i)
		{
			CellType type = grid->cellType(i, j);
			Color fill_color;
			if (type == LIQUID)
				fill_color = Color(0.7,0.7,1);
			else if(type == AIR)
				fill_color = Color(1,1,1);
			else // SOLID
				fill_color = Color(0.5,0.5,0.5);
			_canvas->setFillColor(fill_color);
			_canvas->fillRectangle(
				- translate_x + i * cell_size_x,
				- translate_y + j * cell_size_y,
				- translate_x + (i + 1) * cell_size_x,
				- translate_y + (j + 1) * cell_size_y);
		}
	}
}

void Renderer::renderColorToCanvas(const MacGrid* grid)
{
	double scale_x = _canvas->width() / (_x_max - _x_min);
	double scale_y = _canvas->height() / (_y_max - _y_min);

	double translate_x = 0.5 * (_x_min * _canvas->width());
	double translate_y = 0.5 * (_y_min * _canvas->height());

	// Cell size in pixels
	double cell_size_x = grid->deltaX() * scale_x;
	double cell_size_y = grid->deltaY() * scale_y;

	int w = _canvas->width();
	int h = _canvas->height();

	// Render all velocities as lines
	for (int j = 0; j < h; ++j)
	{
		for (int i = 0; i < w; ++i)
		{
			double x = (i - translate_x) / scale_x;
			double y = (j - translate_y) / scale_y;
			double color = grid->colorInterpolated(x, y);
			Color c(color, color, color);
			_canvas->setPixel(i, j, c);
		}
	}
}

void Renderer::renderGridVelocitiesToCanvas(const MacGrid* grid)
{
	int grid_size_x = grid->sizeX();
	int grid_size_y = grid->sizeY();

	double scale_x = _canvas->width() / (_x_max - _x_min);
	double scale_y = _canvas->height() / (_y_max - _y_min);

	double translate_x = 0.5 * (_x_min * _canvas->width());
	double translate_y = 0.5 * (_y_min * _canvas->height());

	// Cell size in pixels
	double cell_size_x = grid->deltaX() * scale_x;
	double cell_size_y = grid->deltaY() * scale_y;

	double line_scale = 20;
	// Render all velocities as lines
	for (int j = 0; j < grid_size_y; ++j)
	{
		for (int i = 0; i < grid_size_x; ++i)
		{
			double vel_x = grid->velX(i, j);
			double vel_y = grid->velY(i, j);
			double vel_norm = sqrt(vel_x * vel_x + vel_y * vel_y);
			// Since we draw from the center of each cell, add 0.5
			int from_x = - translate_x + (0.5 + i) * cell_size_x;
			int from_y = - translate_y + (0.5 + j) * cell_size_y;
			int to_x = from_x + vel_x * line_scale;
			int to_y = from_y + vel_y * line_scale;

			_canvas->setLineColor(Color(vel_x / vel_norm, vel_y / vel_norm,0));

			_canvas->drawLine(from_x, from_y, to_x, to_y);
		}
	}
}

void Renderer::renderParticlesToCanvas(const MarkerParticleSet* particle_set)
{
	double scale_x = _canvas->width() / (_x_max - _x_min);
	double scale_y = _canvas->height() / (_y_max - _y_min);

	double translate_x = 0.5 * (_x_min * _canvas->width());
	double translate_y = 0.5 * (_y_min * _canvas->height());

	_canvas->setLineColor(Color(0,0,0));
	_canvas->setFillColor(Color(0,0,0));

	for (MarkerParticle* iter = particle_set->getFirst();
		iter != nullptr;
		iter = iter->next)
	{
		// Position in pixel coordinates
		int pos_x = - translate_x + scale_x * iter->posX();
		int pos_y = - translate_y + scale_y * iter->posY();
		
		_canvas->drawPoint(pos_x, pos_y, 6);
	}
}

void Renderer::writeCanvasToPpm(const char* file_path)
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
