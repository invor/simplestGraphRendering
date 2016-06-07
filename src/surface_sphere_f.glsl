#version 130

in vec3 position;
in vec3 normal;

out vec4 frag_colour;

void main()
{
    float lambert = max(0.0,dot(normal,normalize(-position)));

    frag_colour = vec4( vec3(lambert) * 0.9, 1.0);
}