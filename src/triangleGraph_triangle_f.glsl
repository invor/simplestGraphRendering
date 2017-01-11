#version 130

uniform float transparency;

in vec4 colour;

out vec4 fragColour;

void main()
{
    fragColour = vec4(colour[0], colour[1], colour[2], colour[3]*transparency);
}