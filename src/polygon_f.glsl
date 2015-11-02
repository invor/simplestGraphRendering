#version 130

out vec4 fragColor;

void main()
{
	vec3 out_color = vec3(1.0);

	fragColor = vec4(out_color,1.0);
}