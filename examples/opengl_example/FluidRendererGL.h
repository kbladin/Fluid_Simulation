#ifndef FLUID_RENDERER_GL_H
#define FLUID_RENDERER_GL_H

#include <sge/core/simple_graphics_engine.h>
#include <sge/core/shader_manager.h>
#include <sge/core/mesh.h>
#include <sge/core/bounding_box.h>
#include <sge/core/controller.h>
#include <sge/core/new_mesh.h>

#include <FluidDomain.h>

#include <memory>

using namespace sge::core;

class FluidMesh : public Object3D
{
public:
	FluidMesh(MyFloat length_x, MyFloat length_y);
	~FluidMesh();
	
	void updateState(const FluidDomain& fluid_domain);
	virtual void execute();
	
	bool intersects(glm::vec3 origin, glm::vec3 direction, glm::vec2* st) const;
private:
	std::shared_ptr<NewCPUPointCloud> _mesh;
	BoundingBox _aabb;
	float _color_blend;
};

// This class extends SimpleGraphicsEngine,
// Before initializing this object, an OpenGL context must be created
class FluidRendererGL : public SimpleGraphicsEngine
{
public:
	FluidRendererGL(
		int size_x, int size_y, MyFloat length_x, MyFloat length_y);
	~FluidRendererGL();

	Controller& controller();

	void renderFluid(const FluidDomain& fluid_domain);
	bool intersectsFluidMesh(glm::vec2 ndc_position, glm::vec2* st) const;
private:
	FluidMesh _fluid_mesh;
	SphericalController _controller;
};

#endif
