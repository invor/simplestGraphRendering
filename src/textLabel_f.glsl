#version 130

uniform float label_charCount;
uniform sampler2D fontAtlas_tx2D;
uniform sampler2D label_text_tx2D;

in vec2 uv;

out vec4 fragColor;

void main()
{
	vec2 label_text_uv = vec2( 0.5 + floor(uv.x * label_charCount), 0.5)/label_charCount;
	
	vec2 label_uv = texture(label_text_tx2D,label_text_uv).xy + ( vec2(fract(uv.x * label_charCount),uv.y) * vec2(1.0/16.0,1.0/6.0) );

	float character_mask = texture(fontAtlas_tx2D,label_uv).r;
	
	fragColor = vec4(character_mask);
}