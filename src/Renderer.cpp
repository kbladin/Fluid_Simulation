#include <Renderer.h>

Renderer::Renderer(MyFloat x_min, MyFloat y_min, MyFloat x_max, MyFloat y_max)
	: _x_min(x_min)
	, _y_min(y_min)
	, _x_max(x_max)
	, _y_max(y_max)
{

}

Renderer::~Renderer()
{

}

void Renderer::clearCanvas(Canvas& canvas)
{
	canvas.fill(Color(1,1,1));
}

void Renderer::renderGridCellsToCanvas(const MacGrid& grid, Canvas& canvas)
{
	int grid_size_x = grid.sizeX();
	int grid_size_y = grid.sizeY();

	MyFloat scale_x = canvas.width() / (_x_max - _x_min);
	MyFloat scale_y = canvas.height() / (_y_max - _y_min);

	MyFloat translate_x = 0.5 * (_x_min * canvas.width());
	MyFloat translate_y = 0.5 * (_y_min * canvas.height());

	// Cell size in pixels
	MyFloat cell_size_x = grid.deltaX() * scale_x;
	MyFloat cell_size_y = grid.deltaY() * scale_y;

	canvas.setLineColor(Color(1,1,1));

	for (int j = 0; j < grid_size_y; ++j)
	{
		for (int i = 0; i < grid_size_x; ++i)
		{
			CellType type = grid.cellType(i, j);
			Color fill_color;
			if (type == LIQUID)
				fill_color = Color(0.7,0.7,1);
			else if(type == AIR)
				fill_color = Color(1,1,1);
			else // SOLID
				fill_color = Color(0.5,0.5,0.5);
			canvas.setFillColor(fill_color);
			canvas.fillRectangle(
				- translate_x + i * cell_size_x,
				- translate_y + j * cell_size_y,
				- translate_x + (i + 1) * cell_size_x,
				- translate_y + (j + 1) * cell_size_y);
		}
	}
}

void Renderer::renderLevelSetFunctionValuesToCanvas(const LevelSet& level_set, Canvas& canvas)
{
	int grid_size_x = level_set.sizeX();
	int grid_size_y = level_set.sizeY();

	MyFloat scale_x = canvas.width() / (_x_max - _x_min);
	MyFloat scale_y = canvas.height() / (_y_max - _y_min);

	MyFloat translate_x = 0.5 * (_x_min * canvas.width());
	MyFloat translate_y = 0.5 * (_y_min * canvas.height());

	// Cell size in pixels
	MyFloat cell_size_x = level_set.deltaX() * scale_x;
	MyFloat cell_size_y = level_set.deltaY() * scale_y;

	canvas.setLineColor(Color(1,1,1));

	for (int j = 0; j < grid_size_y; ++j)
	{
		for (int i = 0; i < grid_size_x; ++i)
		{
			MyFloat value = level_set.value(i, j);
			Color fill_color = Color(value, value, value);
			canvas.setFillColor(fill_color);
			canvas.fillRectangle(
				- translate_x + i * cell_size_x,
				- translate_y + j * cell_size_y,
				- translate_x + (i + 1) * cell_size_x,
				- translate_y + (j + 1) * cell_size_y);
		}
	}
}

void Renderer::renderGridVelocitiesToCanvas(const MacGrid& grid, Canvas& canvas)
{
	int grid_size_x = grid.sizeX();
	int grid_size_y = grid.sizeY();

	MyFloat scale_x = canvas.width() / (_x_max - _x_min);
	MyFloat scale_y = canvas.height() / (_y_max - _y_min);

	MyFloat translate_x = 0.5 * (_x_min * canvas.width());
	MyFloat translate_y = 0.5 * (_y_min * canvas.height());

	// Cell size in pixels
	MyFloat cell_size_x = grid.deltaX() * scale_x;
	MyFloat cell_size_y = grid.deltaY() * scale_y;

	MyFloat line_scale = 5;
	// Render all velocities as lines
	for (int j = 0; j < grid_size_y; ++j)
	{
		for (int i = 0; i < grid_size_x; ++i)
		{
            MyFloat vel_x = grid.velX(i, j);
            MyFloat vel_y = grid.velY(i, j);
			// Since we draw from the center of each cell, add 0.5
			int from_x = - translate_x + (0.5 + i) * cell_size_x;
			int from_y = - translate_y + (0.5 + j) * cell_size_y;
			int to_x = from_x + vel_x * line_scale;
			int to_y = from_y + vel_y * line_scale;

			// The hue should change with the speed
			float length_in_pixels = sqrt(
				(from_x - to_x) * (from_x - to_x) +
				(from_y - to_y) * (from_y - to_y));

			float length_scaled = length_in_pixels / (canvas.width() / 40);

			MyFloat blue = CLAMP(1 - length_scaled, 0 , 1);
			MyFloat green = (CLAMP(length_scaled, 0 , 1) - CLAMP(length_scaled - 1, 0 , 1));
			MyFloat red = CLAMP(length_scaled - 1, 0 , 1);
	
			canvas.setLineColor(Color(red, green, blue));

			canvas.drawLine(from_x, from_y, to_x, to_y);
		}
	}
}

void Renderer::renderParticlesToCanvas(const MarkerParticleSet& particle_set, Canvas& canvas)
{
	MyFloat scale_x = canvas.width() / (_x_max - _x_min);
	MyFloat scale_y = canvas.height() / (_y_max - _y_min);

	MyFloat translate_x = 0.5 * (_x_min * canvas.width());
	MyFloat translate_y = 0.5 * (_y_min * canvas.height());

	canvas.setLineColor(Color(0,0.5,0.8));
	canvas.setFillColor(Color(0,0.5,0.8));

	for (auto it = particle_set.begin(); it != particle_set.end(); it++)
	{
		// Position in pixel coordinates
		int pos_x = - translate_x + scale_x * it->posX();
		int pos_y = - translate_y + scale_y * it->posY();
		
		canvas.drawPoint(pos_x, pos_y, 3);
	}
}

void Renderer::writeCanvasToPpm(const char* file_path, Canvas& canvas)
{
	int w = canvas.width();
	int h = canvas.height();

	unsigned char* byte_data = new unsigned char[w * h * 3];
	unsigned char byte_value_r;
	unsigned char byte_value_g;
	unsigned char byte_value_b;

	// Copy and cast data
	for (int j = 0; j < h; ++j)
	{
		for (int i = 0; i < w; ++i)
		{ 
			Color c = canvas.pixel(i,j);
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
