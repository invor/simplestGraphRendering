#version 130

uniform float label_charCount;
uniform sampler2D fontAtlas_tx2D;
uniform sampler2D label_text_tx2D;

in vec2 uv;

out vec4 fragColor;

void main()
{
	vec2 label_text_uv = vec2( 0.5 + floor(uv.x * label_charCount), 0.5)/label_charCount;
	
	//vec2 label_uv = texture(label_text_tx2D,label_text_uv).xy + ( vec2( fract(uv.x * label_charCount) ,uv.y) * vec2(1.0/16.0,1.0/6.0) );
	vec2 label_uv = texelFetch(label_text_tx2D,ivec2(floor(uv.x * label_charCount),0),0).xy +
						( vec2( fract(uv.x * label_charCount) , clamp(uv.y,0.05,0.95) ) * vec2(1.0/16.0,1.0/6.0) );

	float character_mask = texture(fontAtlas_tx2D,label_uv).r;
	
	vec3 char_color = vec3(0.0);
	fragColor = vec4(char_color*character_mask,character_mask);
}