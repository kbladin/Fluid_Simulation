#version 410 core

// In data
// Index to pick the current particle
layout(location = 0) in vec3 position;

// Uniform data
// Transform matrices
uniform mat4 M;
uniform mat4 V;
uniform mat4 P;

void main(){
	// Set camera position
	vec4 vertex_position_viewspace = V * M * vec4(position ,1);

	// Position and size of point
	gl_Position = P * vertex_position_viewspace;
	gl_PointSize = 20;
}