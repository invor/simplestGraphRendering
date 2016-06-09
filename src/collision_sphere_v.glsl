#version 140

#define PI 3.141592653589793238462643383279502884197169399375105820

uniform mat4 view_matrix;
uniform mat4 projection_matrix;

uniform sampler2D data_tx2D;

uniform float time;
uniform float grow_factor;
uniform uint id_offset;
uniform int mode;

in vec3 v_position;

flat out ivec2 data_idx;
flat out int instance_id;
out float collision_time;
out float distance_to_center;
out float current_radius;

void main()
{
    instance_id = gl_InstanceID;
    int id = gl_InstanceID + int(id_offset);
    int x_idx = id - int(floor(id/8096.0) * 8096);
    int y_idx = int(floor(id/8096.0));
    data_idx = ivec2(x_idx,y_idx);
    vec4 sphere_data = texelFetch(data_tx2D,data_idx,0);
    vec2 v_geoCoords = vec2(sphere_data.y,sphere_data.x);
    float radius = sphere_data.z;
    collision_time = sphere_data.w;
    
    float lat_sin = sin( (PI/180.0) * v_geoCoords.y);
	float lon_sin = sin( (PI/180.0) * v_geoCoords.x);
	
	float lat_cos = cos( (PI/180.0) * v_geoCoords.y);
	float lon_cos = cos( (PI/180.0) * v_geoCoords.x);
	
    float r = 1.001; //6378137.0;
	
	vec3 sphere_center = vec3( lon_sin * lat_cos * r,
								lat_sin * r,
								lat_cos * lon_cos * r );
	
    
    vec3 sphere = v_position * radius * (time/collision_time);
    current_radius = radius * (time/collision_time);

    vec3 position = sphere_center + sphere;

    if( length(position) >= 1.0 && mode == 1)
    {
        position = normalize(position) * (1.001 + (0.005 * radius * (time/collision_time)));
        distance_to_center = length(position-sphere_center);
    }
    else
    {
        distance_to_center = dot(normalize(sphere),normalize(sphere_center));
    }

    gl_Position = projection_matrix * view_matrix * vec4(position,1.0);
}