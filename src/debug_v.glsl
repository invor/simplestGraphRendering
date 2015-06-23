#version 130

in vec3 v_position;

uniform mat4 view_matrix;
uniform mat4 projection_matrix;

out float view_depth;

void main()
{
	view_depth = (view_matrix * vec4(v_position,1.0)).z;
	
	gl_Position = projection_matrix * view_matrix * vec4(v_position,1.0);
	
	//gl_Position = vec4(v_position,1.0);
}