#ifndef FLUID_RENDERER_GL_H
#define FLUID_RENDERER_GL_H

#include <SGE/SimpleGraphicsEngine.h>
#include <SGE/ShaderManager.h>
#include <SGE/Mesh.h>
#include <SGE/BoundingBox.h>
#include <SGE/Controller.h>

#include <FluidDomain.h>

#include <memory>

class FluidMesh : public Object3D
{
public:
	FluidMesh(MyFloat length_x, MyFloat length_y);
	~FluidMesh();
	
	void updateState(const FluidDomain& fluid_domain);
	virtual void execute();
	
	bool intersects(glm::vec3 origin, glm::vec3 direction, glm::vec2* st) const;
private:
	CPUPointCloudMesh _mesh;
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

	Controller* controller();

	void renderFluid(const FluidDomain& fluid_domain);
	bool intersectsFluidMesh(glm::vec2 ndc_position, glm::vec2* st) const;
private:
	FluidMesh _fluid_mesh;
	SphericalController _controller;
};

#endif
