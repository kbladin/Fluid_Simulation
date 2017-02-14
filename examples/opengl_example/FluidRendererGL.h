#ifndef FLUID_RENDERER_GL_H
#define FLUID_RENDERER_GL_H

#include <sge/core/simple_graphics_engine.h>
#include <sge/core/shader_manager.h>
#include <sge/core/mesh.h>
#include <sge/core/bounding_box.h>
#include <sge/core/controller.h>
#include <sge/core/mesh.h>
#include <sge/core/simple_forward_3d_renderer.h>

#include <FluidDomain.h>
#include <FluidSolver.h>

#include <memory>

using namespace sge::core;

class RenderableFluidMesh : public RenderableForward
{
public:
	RenderableFluidMesh(MyFloat length_x, MyFloat length_y);
	~RenderableFluidMesh();
	
	void updateState(const FluidDomain& fluid_domain);
	virtual void render(const UsefulRenderData& render_data) override;
	
	bool intersects(glm::vec3 origin, glm::vec3 direction, glm::vec2* st) const;
private:
  std::shared_ptr<ShaderProgram> _program;
	std::shared_ptr<CPUPointCloud> _mesh;
	BoundingBox _aabb;
	float _color_blend;
};

// This class extends SimpleGraphicsEngine,
// Before initializing this object, an OpenGL context must be created
class FluidRendererGL : public SimpleGraphicsEngine
{
public:
	FluidRendererGL(
		int size_x, int size_y, MyFloat length_x, MyFloat length_y,
    int wondow_width, int window_height);
	~FluidRendererGL();

	FluidDomain& fluidDomain();
  SimpleForward3DRenderer& renderer() { return _renderer; };

  void update(double dt);
	bool intersectsFluidMesh(glm::vec2 ndc_position, glm::vec2* st) const;
private:
  SimpleForward3DRenderer _renderer;

  FluidDomain _fluid_domain;
  FluidSolverMemoryPool _mem_pool;
  FluidSolver _fluid_solver;

	RenderableFluidMesh _fluid_mesh;
};

#endif
