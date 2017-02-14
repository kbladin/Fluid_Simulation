#include "FluidRendererGL.h"

#include <gl/glew.h>
#include <gl/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

RenderableFluidMesh::RenderableFluidMesh(MyFloat length_x, MyFloat length_y) :
	_aabb(glm::vec3(0, 0, 0), glm::vec3(length_x, length_y, 0))
{
	_program = std::make_shared<ShaderProgram>(
		"render_cpu_particles",
		(std::string(PROJECT_SOURCE_DIR) + "/shaders/render_cpu_particles.vert").c_str(),
		nullptr,
		nullptr,
		nullptr,
		(std::string(PROJECT_SOURCE_DIR) + "/shaders/render_cpu_particles.frag").c_str());

  auto positions = new std::vector<glm::vec3>;
  positions->resize(1);
  _mesh = std::make_shared<CPUPointCloud>(positions);
}

RenderableFluidMesh::~RenderableFluidMesh()
{

}

void RenderableFluidMesh::updateState(const FluidDomain& fluid_domain)
{
	std::vector<glm::vec3> points;
	points.reserve(fluid_domain.markerParticleSet().size());
	for (auto it = fluid_domain.markerParticleSet().begin();
		it != fluid_domain.markerParticleSet().end(); it++)
	{
		points.push_back(glm::vec3(it->posX(), it->posY(),0));
	}
	_mesh->update(points);
	_color_blend = pow(fluid_domain.picRatio(), 0.2);
}

void RenderableFluidMesh::render(const UsefulRenderData& render_data)
{
  _program->pushUsage();
	// Input to the shader
	glUniformMatrix4fv(
		glGetUniformLocation(_program->id(), "M"),
		1,
		GL_FALSE,
		&absoluteTransform()[0][0]);
	glUniformMatrix4fv(
		glGetUniformLocation(_program->id(), "V"),
		1,
		GL_FALSE,
		&render_data.camera.viewTransform()[0][0]);
	glUniformMatrix4fv(
		glGetUniformLocation(_program->id(), "P"),
		1,
		GL_FALSE,
		&render_data.camera.projectionTransform()[0][0]);
	glUniform1f(glGetUniformLocation(_program->id(), "color_blend"), _color_blend);
    
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_DEPTH_TEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  
	_mesh->render();
  _program->popUsage();
}

bool RenderableFluidMesh::intersects(glm::vec3 origin, glm::vec3 direction, glm::vec2* st) const
{
	glm::vec3 origin_model_space = glm::vec3(
		glm::inverse(absoluteTransform()) * glm::vec4(origin, 1));
	glm::vec3 direction_model_space = glm::vec3(
		glm::inverse(absoluteTransform()) * glm::vec4(direction, 0));
    auto intersection = _aabb.intersects(origin_model_space, direction_model_space);
	if (intersection.first)
	{
        glm::vec3 intersection_point = origin_model_space + intersection.second * direction_model_space;
        *st = glm::vec2(intersection_point);
		return true;
	}
	else
	{
		return false;
	}
}

FluidRendererGL::FluidRendererGL(int size_x, int size_y, MyFloat length_x, MyFloat length_y, int window_width, int window_height) :
	SimpleGraphicsEngine(),
  _renderer(camera(), window_width, window_height),
  _fluid_domain(size_x, size_y, length_x, length_y, 0.005, 0.02),
  _mem_pool(_fluid_domain),
	_fluid_solver(_mem_pool),
  _fluid_mesh(length_x, length_y)
{
// Initialize
	
	// Setup
    //_fluid_domain.addFluidSource(FluidSource( { 4.0 / GRID_X_SIZE, 0.1, 3.0 / GRID_Y_SIZE, 1 - 4.0 / GRID_Y_SIZE }, DELTA_X, DELTA_Y, 0.0, 0.0, 0.0, 1));
  

	scene.addChild(_fluid_mesh);
	//view_space.addChild(&_fluid_mesh);

    //_fluid_mesh.transform_matrix *= glm::scale(2.0f * glm::vec3(1/length_x, 1/length_y, 1));
    _fluid_mesh.setTransform(glm::translate(glm::vec3(-length_x/2, -length_y/2, 0)));
}

FluidRendererGL::~FluidRendererGL()
{
}

FluidDomain& FluidRendererGL::fluidDomain()
{
	return _fluid_domain;
}

void FluidRendererGL::update(double dt_frame)
{
  SimpleGraphicsEngine::update(dt_frame);
  
  MyFloat seconds_per_frame = 1 / 60.0;
  MyFloat dt;
  // Simulate
  for (MyFloat frame_time = 0; frame_time < seconds_per_frame; frame_time += dt)
  {
    // Calculate dt (for now just set it)
    dt = 0.005;
    dt = CLAMP(dt, 0, seconds_per_frame - frame_time);
    // Update fluid domain (creates new fluid from sources)
    _fluid_domain.update(dt);

    // Solve
    //fluid_solver.stepSemiLagrangian(fluid_domain, dt);
    _fluid_solver.stepPICFLIP(_fluid_domain, dt);
  }
  _fluid_mesh.updateState(_fluid_domain);
  _renderer.render(scene);
}

bool FluidRendererGL::intersectsFluidMesh(glm::vec2 ndc_position, glm::vec2* st) const
{
	auto origin_and_direction = perspective_camera.unproject(ndc_position);
	return _fluid_mesh.intersects(origin_and_direction.first, origin_and_direction.second, st);
}
