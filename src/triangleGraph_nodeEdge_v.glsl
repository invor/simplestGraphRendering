#version 130

#define PI 3.141592653589793238462643383279502884197169399375105820

uniform mat4 view_matrix;
uniform mat4 projection_matrix;

in vec2 v_geoCoords;
in vec4 v_colour;

out vec4 colour;

void main()
{
    colour = v_colour;

    float lat_sin = sin( (PI/180.0) * v_geoCoords.y);
	float lon_sin = sin( (PI/180.0) * v_geoCoords.x);
	
	float lat_cos = cos( (PI/180.0) * v_geoCoords.y);
	float lon_cos = cos( (PI/180.0) * v_geoCoords.x);
	
	//float r = 1.0; //6378137.0;
    float r = 1.001; //6378137.0;
	
	vec3 world_position = vec3( lon_sin * lat_cos * r,
								lat_sin * r,
								lat_cos * lon_cos * r );
								
	gl_Position = projection_matrix * view_matrix * vec4(world_position,1.0);
}