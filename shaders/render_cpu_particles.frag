#version 410 core

// Out data
out vec4 color;

void main(){
  vec2 coord = (gl_PointCoord - vec2(0.5))*2;  // From [0,1] to [-0.5,0.5]
  float r = length(coord);
  if(r > 1) // Outside of circle radius?
    discard;
  color = vec4(0.05,0.2,0.5, (1 - r));
}