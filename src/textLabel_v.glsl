#version 130

#define PI 3.1415926535897932384

uniform mat4 view_matrix;
uniform mat4 projection_matrix;

uniform vec2 label_geoCoords;
uniform float label_scale;

in vec2 v_position;
in vec2 v_uv;

out vec2 uv;

void main()
{	
	// compute label position on unit sphere
	float lat_sin = sin( (PI/180.0) * label_geoCoords.y);
	float lon_sin = sin( (PI/180.0) * label_geoCoords.x);
	
	float lat_cos = cos( (PI/180.0) * label_geoCoords.y);
	float lon_cos = cos( (PI/180.0) * label_geoCoords.x);
	
	float r = 1.0; //6378137.0;
	
	vec3 world_position = vec3( lon_sin * lat_cos * r,
								lat_sin * r,
								lat_cos * lon_cos * r );
	
	// Original Screen-Space Version:
	/*
	// transform vertex position in DCS to match character position
	vec4 dcs_position = projection_matrix * view_matrix * vec4(world_position,1.0);
	// build base quad for each character in NDCS and add horizontal offset of char position in string
	dcs_position += vec4(v_position*dcs_position.w*label_scale,0.0,0.0);
	// to center lable on geoCoords, shift label to the left first
	dcs_position -= vec4( (0.06*label_scale) * (label_charCount/2.0) * dcs_position.w,0.0,0.0,0.0); 
	dcs_position += vec4( (0.06*label_scale) * float(gl_InstanceID) * dcs_position.w,0.0,0.0,0.0);
	*/
	
	// transform vertex position in DCS to match character position
	float ccs_scale = 0.25;
	
	vec4 ccs_position = view_matrix * vec4(world_position,1.0);
	// build base quad for each character in NDCS and add horizontal offset of char position in string
	ccs_position += vec4(v_position*label_scale*ccs_scale,0.0,0.0);
	
	vec4 dcs_position = projection_matrix * ccs_position;
	
	uv = v_uv;
	
	gl_Position = dcs_position;
}