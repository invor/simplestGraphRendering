#version 330

in float view_depth;

out vec4 gl_FragColor;

void main()
{
	gl_FragColor = vec4(0.0,0.0, (10+view_depth)/10.0 ,1.0);
}