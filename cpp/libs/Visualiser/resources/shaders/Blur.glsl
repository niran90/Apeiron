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

in vec2 v_texture_coordinate;

out vec4 fragment_colour;

uniform sampler2D u_texture;
uniform bool      u_horizontal;
uniform float     u_weights[5] = float[] (0.227027, 0.1945946, 0.1216216, 0.054054, 0.016216);

void main()
{
    fragment_colour = texture(u_texture, v_texture_coordinate);

    vec2 texture_offset = 1.0 / textureSize(u_texture, 0);    // Size of a single texel
    vec3 blurred_colour = fragment_colour.rgb * u_weights[0]; // Current fragment contribution

    if(u_horizontal) // Horizontal blur pass
    {
        for(int i = 1; i < 5; ++i)
        {
            blurred_colour += u_weights[i] * texture(u_texture, v_texture_coordinate + vec2(texture_offset.x * i, 0.0)).rgb;
            blurred_colour += u_weights[i] * texture(u_texture, v_texture_coordinate - vec2(texture_offset.x * i, 0.0)).rgb;
        }
    }
    else // Vertical blur pass
    {
        for(int i = 1; i < 5; ++i)
        {
            blurred_colour += u_weights[i] * texture(u_texture, v_texture_coordinate + vec2(0.0, texture_offset.y * i)).rgb;
            blurred_colour += u_weights[i] * texture(u_texture, v_texture_coordinate - vec2(0.0, texture_offset.y * i)).rgb;
        }
    }

    fragment_colour = vec4(blurred_colour, fragment_colour.a);
}