#version 330

in vec3 v_position;

uniform mat4 view_matrix;
uniform mat4 projection_matrix;

void main()
{
	gl_Position = projection_matrix * view_matrix * vec4(v_position,1.0);
	
	//gl_Position = vec4(v_position,1.0);
}