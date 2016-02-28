#version 130

uniform float transparency;

in vec3 colour;

out vec4 fragColour;

void main()
{
    fragColour = vec4(colour,transparency);
}