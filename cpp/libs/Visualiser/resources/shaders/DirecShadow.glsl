#shader vertex
#version 460 core

layout(location = 0) in vec3 position;

uniform mat4 u_model_matrix;
uniform mat4 u_dlight_space_matrix;

void main()
{
    gl_Position = u_dlight_space_matrix * u_model_matrix * vec4(position, 1.0);
}



#shader fragment
#version 460 core

void main()
{

}