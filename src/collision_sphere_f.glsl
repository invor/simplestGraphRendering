#version 130

uniform float time;

in float collision_time;

out vec4 frag_colour;

void main()
{
    float v = time/collision_time;
    
    float r = min(v * 2.0,1.0);
    float g = (v < 0.5) ? 1.0 : 1.0 - (v-0.5)*2.0;
    
    frag_colour = vec4(r,g,0.0,1.0);
}