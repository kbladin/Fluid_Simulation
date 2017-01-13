#include "FluidRendererGL.h"

#include <gl/glew.h>
#include <gl/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

FluidMesh::FluidMesh()
{
}

FluidMesh::~FluidMesh()
{
}

void FluidMesh::update(const MarkerParticleSet& particle_set)
{
	std::vector<glm::vec3> points;
	points.resize(particle_set.size());
	for (auto it = particle_set.begin(); it != particle_set.end(); it++)
	{
		points.push_back(glm::vec3(it->posX(), it->posY(),0));
	}
	_mesh.update(points);
}

void FluidMesh::render(glm::mat4 M)
{
    GLuint program_ID = ShaderManager::instance()->getShader("render_cpu_particles");
	
	glUseProgram(program_ID);

	// Input to the shader
	glUniformMatrix4fv(
		glGetUniformLocation(program_ID, "M"),
		1,
		GL_FALSE,
		&transform_matrix[0][0]);
    
	glDisable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glDisable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);

	_mesh.render();
}

FluidRendererGL::FluidRendererGL(int size_x, int size_y, MyFloat length_x, MyFloat length_y) :
	SimpleGraphicsEngine(size_x, size_y)
{
	ShaderManager::instance()->loadShader(
		"render_cpu_particles",
		(std::string(PROJECT_SOURCE_DIR) + "/shaders/render_cpu_particles.vert").c_str(),
		nullptr,
		nullptr,
		nullptr,
		(std::string(PROJECT_SOURCE_DIR) + "/shaders/render_cpu_particles.frag").c_str());

	_fluid_mesh = std::make_shared<FluidMesh>();
	
	//camera->addToShader(ShaderManager::instance()->getShader("render_cpu_particles"));
	viewspace_ortho_camera->addToShader(ShaderManager::instance()->getShader("render_cpu_particles"));

	//scene->addChild(_fluid_mesh);
	view_space->addChild(_fluid_mesh.get());

    _fluid_mesh->transform_matrix *= glm::scale(2.0f * glm::vec3(1/length_x, 1/length_y, 1));
    _fluid_mesh->transform_matrix *= glm::translate(glm::vec3(-length_x/2, -length_y/2, 0));
    
    //camera->transform_matrix = glm::inverse(glm::lookAt(glm::vec3(2.0f,1.0f,2.0f), glm::vec3(1.0f,0.0f,0.0f), glm::vec3(0.0f,1.0f,0.0f)));
}

FluidRendererGL::~FluidRendererGL()
{
}

void FluidRendererGL::renderParticles(const MarkerParticleSet& particle_set)
{
	_fluid_mesh->update(particle_set);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable( GL_BLEND );
	render();
}
