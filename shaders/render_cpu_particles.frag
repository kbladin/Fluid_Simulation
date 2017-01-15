#version 410 core

// Out data
out vec4 color;

uniform float color_blend;

void main(){
  vec2 coord = (gl_PointCoord - vec2(0.5))*2;  // From [0,1] to [-0.5,0.5]
  float r = length(coord);
  if(r > 1) // Outside of circle radius?
    discard;
  vec4 color1 = vec4(0.2,0.3,0.5, (1 - r));
  vec4 color2 = vec4(0.3, 0.5,0.2, (1 - r));
  color = color1 * color_blend + color2 * (1 - color_blend);
}