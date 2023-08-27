#shader vertex
#version 460 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 tangent;
layout(location = 3) in vec3 colour;
layout(location = 4) in vec2 texture_coordinate;

out vec2 v_texture_coordinate;

void main()
{
    gl_Position = vec4(position.x, position.y, 0.0, 1.0);
    v_texture_coordinate = texture_coordinate;
}



#shader fragment
#version 460 core

const float Gamma = 2.2;

in vec2 v_texture_coordinate;

out vec4 fragment_colour;

uniform float     u_exposure;
uniform sampler2D u_texture;

vec4 ToneMap(const vec4 colour) { return vec4(vec3(1.0) - exp(-u_exposure * colour.rgb), colour.a); }

vec4 GammaCorrect(const vec4 colour) { return vec4(pow(colour.rgb, vec3(1.0f / Gamma)), colour.a); }

void main()
{
    fragment_colour = texture(u_texture, v_texture_coordinate);
    fragment_colour = ToneMap(fragment_colour);
    fragment_colour = GammaCorrect(fragment_colour);
}