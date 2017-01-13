#ifndef FLUID_RENDERER_GL_H
#define FLUID_RENDERER_GL_H

#include <SGE/SimpleGraphicsEngine.h>
#include <SGE/ShaderManager.h>
#include <SGE/Mesh.h>

#include <MarkerParticleSet.h>

#include <memory>

class FluidMesh : public Object3D
{
public:
	FluidMesh();
	~FluidMesh();
	
	void update(const MarkerParticleSet& particle_set);
	virtual void render(glm::mat4 M);
	
private:
	CPUPointCloudMesh _mesh;
};

// This class extends SimpleGraphicsEngine,
// Before initializing this object, an OpenGL context must be created
class FluidRendererGL : public SimpleGraphicsEngine
{
public:
	FluidRendererGL(int size_x, int size_y, MyFloat length_x, MyFloat length_y);
	~FluidRendererGL();

	void renderParticles(const MarkerParticleSet& particle_set);
private:
	std::shared_ptr<FluidMesh> _fluid_mesh;
};

#endif
