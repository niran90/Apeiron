#shader vertex
#version 460 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 colour;
layout(location = 3) in vec2 texture_coordinate;

out Data
{
   vec3 Normal;
   vec3 Colour;
   vec3 FragmentPosition;
   vec2 TextureCoordinate;
   mat4 ViewMatrix;
   mat4 ProjectionMatrix;
} v_data_out;

uniform mat4 u_model_matrix;
uniform mat4 u_view_matrix;
uniform mat4 u_projection_matrix;

void main()
{
   const vec4 vertex_position = vec4(position, 1.0);
   gl_Position = u_model_matrix * vertex_position;

   v_data_out.Normal = mat3(transpose(inverse(u_model_matrix))) * normal; // Accounts for model rotation and non-uniform scaling
   v_data_out.Colour = colour;
   v_data_out.FragmentPosition = (u_model_matrix * vertex_position).xyz;
   v_data_out.TextureCoordinate = texture_coordinate;
   v_data_out.ViewMatrix = u_view_matrix;
   v_data_out.ProjectionMatrix = u_projection_matrix;
}




#shader geometry
#version 460 core

layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

in Data
{
   vec3 Normal;
   vec3 Colour;
   vec3 FragmentPosition;
   vec2 TextureCoordinate;
   mat4 ViewMatrix;
   mat4 ProjectionMatrix;
} v_data_in[];

out vec3 v_normal;
out vec3 v_colour;
out vec3 v_fragment_position;
out vec2 v_texture_coordinate;

void main()
{
   const vec3 tangent0 = (gl_in[1].gl_Position - gl_in[0].gl_Position).xyz;
   const vec3 tangent1 = (gl_in[2].gl_Position - gl_in[0].gl_Position).xyz;
   const vec3 face_normal = normalize(cross(tangent0, tangent1));

   for(int i = 0; i < gl_in.length(); i++)
   {
//      gl_Position = v_data_in[i].ProjectionMatrix * v_data_in[i].ViewMatrix * (gl_in[i].gl_Position + 0.5*vec4(face_normal, 0.0));
      gl_Position = v_data_in[i].ProjectionMatrix * v_data_in[i].ViewMatrix * gl_in[i].gl_Position;

      v_normal = face_normal; // Flat shading
//      v_normal = v_data_in[i].Normal; // Smooth shading

      v_colour = v_data_in[i].Colour;
      v_fragment_position = v_data_in[i].FragmentPosition;
      v_texture_coordinate = v_data_in[i].TextureCoordinate;
      EmitVertex();
   }
   EndPrimitive();
}




#shader fragment
#version 460 core
#extension GL_OES_standard_derivatives : enable

in vec3 v_normal;
in vec3 v_colour;
in vec3 v_fragment_position;
in vec2 v_texture_coordinate;

out vec4 colour;

struct Light
{
   vec4 Colour;
   float AmbientIntensity;
   float DiffuseIntensity;
};

struct DirectionalLight
{
   Light Base;
   vec3 Direction;
};

struct PointLight
{
   Light Base;
   vec3 Position;
   vec3 AttenuationCoefficients;
};

struct SpotLight
{
   PointLight Point;
   vec3 Direction;
   float CosConeAngle;
};

struct Material
{
   float SpecularIntensity;
   float Smoothness;
};

uniform vec4 u_colour;
uniform vec3 u_camera_position;
uniform sampler2D u_texture;

uniform Material u_material;

const int Max_Point_Lights = 4;
const int Max_Spot_Lights = 4;
uniform int u_point_light_count;
uniform int u_spot_light_count;
uniform DirectionalLight u_directional_light;
uniform PointLight u_point_lights[Max_Point_Lights];
uniform SpotLight u_spot_lights[Max_Spot_Lights];

vec4 CalculateLightByDirection(Light _light, vec3 _direction)
{
   const vec4 ambient_colour = _light.AmbientIntensity * _light.Colour;

   const float diffuse_factor = max(dot(-v_normal, normalize(_direction)), 0.0f);
   const vec4 diffuse_colour = diffuse_factor * _light.DiffuseIntensity * _light.Colour;

   const vec3 fragment_to_camera = normalize(u_camera_position - v_fragment_position);
   const vec3 reflected_ray = normalize(reflect(_direction, v_normal));
   const float specular_factor = pow(max(dot(fragment_to_camera, reflected_ray), 0.0f), u_material.Smoothness);
   const vec4 specular_colour = specular_factor * u_material.SpecularIntensity * _light.Colour;

   return (ambient_colour + diffuse_colour + specular_colour);
}

vec4 CalculateDirectionalLight()
{
   return CalculateLightByDirection(u_directional_light.Base, u_directional_light.Direction);
}

vec4 CalculatePointLight(PointLight _point_light)
{
   vec3 light_to_fragment = v_fragment_position - _point_light.Position;
   float distance = length(light_to_fragment);
   light_to_fragment = normalize(light_to_fragment);

   const float attenuation = _point_light.AttenuationCoefficients[0] +
   _point_light.AttenuationCoefficients[1] * distance +
   _point_light.AttenuationCoefficients[2] * distance * distance;

   return CalculateLightByDirection(_point_light.Base, light_to_fragment) / attenuation;
}

vec4 CalculatePointLights()
{
   vec4 total_colour = vec4(0.0, 0.0, 0.0, 0.0);
   for(int i = 0; i < u_point_light_count; i++) total_colour += CalculatePointLight(u_point_lights[i]);
   return total_colour;
}

vec4 CalculateSpotLight(SpotLight _spot_light)
{
   vec3 ray_direction = normalize(v_fragment_position - _spot_light.Point.Position);
   float spot_light_factor = dot(ray_direction, _spot_light.Direction);

   if(spot_light_factor > _spot_light.CosConeAngle)
      return (1.0f - (1.0f - spot_light_factor ) / (1.0f - _spot_light.CosConeAngle)) * CalculatePointLight(_spot_light.Point);
   else return vec4(0.0, 0.0, 0.0, 0.0);
}

vec4 CalculateSpotLights()
{
   vec4 total_colour = vec4(0.0, 0.0, 0.0, 0.0);
   for(int i = 0; i < u_spot_light_count; i++) total_colour += CalculateSpotLight(u_spot_lights[i]);
   return total_colour;
}

void main()
{
//   vec4 lighting = CalculatePointLights();
   vec4 lighting = CalculateDirectionalLight() + CalculatePointLights() + CalculateSpotLights();

//   vec4 texture_colour = texture(u_texture, v_texture_coordinate);
//   colour = texture_colour * lighting;

   colour = u_colour * lighting;

//   colour = v_position * lighting;
}