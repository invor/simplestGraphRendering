#version 130

in float view_depth;

out vec4 fragColor;

void main()
{
	fragColor = vec4(0.0,0.0, (10+view_depth)/10.0 ,1.0);
}