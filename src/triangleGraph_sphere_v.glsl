#version 130

#define PI 3.141592653589793238462643383279502884197169399375105820

uniform mat4 model_matrix;
uniform mat4 view_matrix;
uniform mat4 projection_matrix;

in vec3 v_position;

void main()
{					
	gl_Position = projection_matrix * view_matrix * model_matrix * vec4(v_position,1.0);
}