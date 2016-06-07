#version 140

uniform mat4 view_matrix;
uniform mat4 projection_matrix;

in vec3 v_position;

out vec3 position;
out vec3 normal;

void main()
{
    position = (view_matrix * vec4(v_position,1.0)).xyz;
    normal = inverse(transpose(mat3(view_matrix))) * normalize(v_position);

    gl_Position = projection_matrix * view_matrix * vec4(normalize(v_position),1.0);
}