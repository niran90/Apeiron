#shader vertex
#version 460 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 tangent;
layout(location = 3) in vec4 colour;
layout(location = 4) in vec2 texture_coordinate;

// Output varying data
out Data
{
    vec2 TextureCoordinate;
} v_data_out;

uniform mat4 u_model_matrix;

void main()
{
    gl_Position = u_model_matrix * vec4(position, 1.0);
    v_data_out.TextureCoordinate = texture_coordinate;
}




#shader geometry
#version 460 core

layout (triangles) in;
layout (triangle_strip, max_vertices = 18) out;

// Input varying data
in Data
{
    vec2 TextureCoordinate;
} v_data_in[];

in vec2  i_texture_coordinates[];
out vec4 v_fragment_position;
out vec2 v_texture_coordinate;

uniform mat4 u_plight_space_matrices[6];

void main()
{
    for(int face = 0; face < 6; ++face)
    {
        gl_Layer = face; // Index of the cube map face to render to
        for(int i = 0; i < gl_in.length(); ++i)
        {
            v_fragment_position = gl_in[i].gl_Position;
            gl_Position = u_plight_space_matrices[face] * v_fragment_position;
            v_texture_coordinate = v_data_in[i].TextureCoordinate;

            EmitVertex();
        }
        EndPrimitive();
    }
}




#shader fragment
#version 460 core

in vec4 v_fragment_position;
in vec2 v_texture_coordinate;

uniform float     u_point_light_far_plane;
uniform vec3      u_plight_position;
uniform sampler2D u_diffuse_map;

void main()
{
    // Discard fragment if transparent.
    if(texture(u_diffuse_map, v_texture_coordinate).a < 0.0001) discard;

    gl_FragDepth = length(v_fragment_position.xyz - u_plight_position) / u_point_light_far_plane; // Division maps depth to [0,1]
}