#version 130

uniform sampler2D iconAtlas_tx2D;

in vec2 uv;
in vec2 atlasUV;

out vec4 fragColor;

void main()
{
	vec2 icon_uv = atlasUV + uv * vec2(1.0/9.0);

	vec3 icon_rgb = texture(iconAtlas_tx2D,icon_uv).rgb;
	float alpha = (icon_rgb.r < 0.8 && icon_rgb.b < 0.8) ? 1.0 : 0.0;
	
	fragColor = vec4(icon_rgb,alpha);
}