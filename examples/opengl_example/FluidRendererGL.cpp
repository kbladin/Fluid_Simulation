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

FluidRendererGL::FluidRendererGL(int size_x, int size_y, MyFloat length_x, MyFloat length_y) :
	SimpleGraphicsEngine(),
	_fluid_mesh(length_x, length_y),
	_controller(perspective_camera)
{
	scene.addChild(_fluid_mesh);
	//view_space.addChild(&_fluid_mesh);

    //_fluid_mesh.transform_matrix *= glm::scale(2.0f * glm::vec3(1/length_x, 1/length_y, 1));
    _fluid_mesh.setTransform(glm::translate(glm::vec3(-length_x/2, -length_y/2, 0)));
}

FluidRendererGL::~FluidRendererGL()
{
}

Controller& FluidRendererGL::controller()
{
	return _controller;
}

void FluidRendererGL::update(double dt)
{
  SimpleGraphicsEngine::update(dt);
}

void FluidRendererGL::renderFluid(const FluidDomain& fluid_domain)
{
	_fluid_mesh.updateState(fluid_domain);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable( GL_BLEND );
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
	_fluid_mesh.render({camera()});
}

bool FluidRendererGL::intersectsFluidMesh(glm::vec2 ndc_position, glm::vec2* st) const
{
	auto origin_and_direction = perspective_camera.unproject(ndc_position);
	return _fluid_mesh.intersects(origin_and_direction.first, origin_and_direction.second, st);
}
