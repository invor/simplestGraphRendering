#version 130

uniform sampler2D background_tx2D;

in vec2 geoCoords;

out vec4 fragColor;

void main()
{
	vec3 out_color = vec3(1.0);

	out_color = texture(background_tx2D,geoCoords*vec2(10.0,20.0)).rgb;

	fragColor = vec4(out_color,1.0);
}