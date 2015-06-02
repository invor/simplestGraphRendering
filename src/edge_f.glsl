#version 330

in float color;

void main()
{
	vec3 out_color = vec3(0.0);

	
	if(color == 0.0)
		out_color = vec3(0.0,0.0,0.0);
	else if(color == 1.0)
		out_color = vec3(0.3,0.55,0.95);
	else if(color == 2.0)
		out_color = vec3(0.95,0.4,0.4);
	else if(color == 3.0)
		out_color = vec3(0.95,0.75,0.45);
	else if(color == 4.0)
		out_color = vec3(0.95,0.9,0.55);
	else if(color == 5.0)
		out_color = vec3(1.0,1.0,1.0);

	gl_FragColor = vec4(out_color,1.0);
}