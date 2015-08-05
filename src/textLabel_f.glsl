#version 130

uniform sampler2D fontAtlas_tx2D;

in vec2 uv;

out vec4 fragColor;

void main()
{

	float character_mask = texture(fontAtlas_tx2D,uv).r;
	
	fragColor = vec4(character_mask);
}