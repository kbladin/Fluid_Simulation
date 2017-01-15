#include "FluidInteractionHandler.h"

#include <iostream>

FluidInteractionHandler::FluidInteractionHandler(
	FluidDomain* fluid_domain, FluidRendererGL* renderer) :
	_fluid_domain(fluid_domain),
	_renderer(renderer)
{
	_mouse_x = 0;
	_mouse_y = 0;
	_emitter_radius = 0.05;
}

FluidInteractionHandler::~FluidInteractionHandler()
{

}

void FluidInteractionHandler::mousePosCallback(double x, double y)
{
	_mouse_x = x;
	_mouse_y = y;
}

void FluidInteractionHandler::mouseButtonCallback(int button, int action, int mods)
{
	glm::vec2 st;
	if (_renderer->intersectsFluidMesh(glm::vec2(_mouse_x, _mouse_y), &st))
	{
		_fluid_domain->addFluidSource(FluidSource(
		{ st.x - _emitter_radius, st.x + _emitter_radius,
			st.y - _emitter_radius, st.y + _emitter_radius },
		_fluid_domain->deltaX(),
		_fluid_domain->deltaY(),
		0.0, 0.0, 0.0, 1));
	}
	/*
	float pos_x = (_mouse_x / 2 + 0.5) * _fluid_domain->lengthX();
	float pos_y = (_mouse_y / 2 + 0.5) * _fluid_domain->lengthY();
	_fluid_domain->addFluidSource(FluidSource(
		{ pos_x - 0.05f, pos_x + 0.05f, pos_y - 0.05f, pos_y + 0.05f },
		_fluid_domain->deltaX(),
		_fluid_domain->deltaY(),
		0.0, 0.0, 0.0, 1));
	*/
}

void FluidInteractionHandler::mouseScrollCallback(double dx, double dy)
{

}

void FluidInteractionHandler::keyCallback(int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS || action == GLFW_REPEAT)
	{
		switch (key)
		{
			case GLFW_KEY_R :
			{
				_fluid_domain->clearFluidSources();
				_fluid_domain->resetParticleSet();
                break;
			}
			case GLFW_KEY_O :
			{
				_fluid_domain->setPicRatio(_fluid_domain->picRatio() * 0.8);
                break;
			}
			case GLFW_KEY_P :
			{
				_fluid_domain->setPicRatio(_fluid_domain->picRatio() * 1.2);
                break;
			}
			case GLFW_KEY_K :
			{
				_emitter_radius -= 0.01;
                break;
			}
			case GLFW_KEY_L :
			{
				_emitter_radius += 0.01;
                break;
			}
			default: break;
		}
	}
	_emitter_radius = glm::clamp(_emitter_radius, 0.01f, 0.2f);
}

void FluidInteractionHandler::windowSizeCallback(int width, int height)
{
	_renderer->setWindowResolution(width, height);
}
