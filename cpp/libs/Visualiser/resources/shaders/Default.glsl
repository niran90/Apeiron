#shader vertex
#version 460 core

// Global constants
const int Max_Point_Lights = 4;
const int Max_Spot_Lights  = 4;

// Input varying data
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 tangent;
layout(location = 3) in vec4 colour;
layout(location = 4) in vec2 texture_coordinate;

// Output varying data
out Data
{
   vec3 Normal;
   vec4 Colour;
   vec2 TextureCoordinate;

   vec3 PointLightPositions[Max_Point_Lights];
   vec3 SpotLightPositions[Max_Spot_Lights];
   vec3 CameraPosition;
   vec3 FragmentPosition;
   vec4 FragmentPositionDlight;

   mat4 ViewProjMatrix;
   mat3 TBNMatrix;
} v_data_out;

// Uniforms
uniform mat4 u_model_matrix;
uniform mat4 u_view_matrix;
uniform mat4 u_projection_matrix;
uniform mat4 u_dlight_space_matrix;
uniform bool u_use_diffuse_map;
uniform bool u_use_normal_map;
uniform bool u_use_displacement_map;
uniform int  u_point_light_count;
uniform int  u_spot_light_count;
uniform vec3 u_camera_position;
uniform vec3 u_point_light_positions[Max_Point_Lights];
uniform vec3 u_spot_light_positions[Max_Spot_Lights];

void main()
{
   gl_Position = u_model_matrix * vec4(position, 1.0);

   v_data_out.Normal = mat3(transpose(inverse(u_model_matrix))) * normal; // Accounts for model rotation and non-uniform scaling
   v_data_out.Colour = colour;
   v_data_out.TextureCoordinate = texture_coordinate;

   const bool transform_to_tangent_space = u_use_diffuse_map && u_use_normal_map;
   v_data_out.TBNMatrix = mat3(1.0f);
   mat3 inverse_TBN_matrix = v_data_out.TBNMatrix;
   if(transform_to_tangent_space)
   {
      const vec3 T = mat3(transpose(inverse(u_model_matrix))) * tangent;
      const vec3 N = v_data_out.Normal;
      const vec3 B = normalize(cross(N, T));

      v_data_out.TBNMatrix = mat3(T, B, N);
      inverse_TBN_matrix   = transpose(v_data_out.TBNMatrix);
   }

   for(int i = 0; i < u_point_light_count; ++i) v_data_out.PointLightPositions[i] = transform_to_tangent_space ? inverse_TBN_matrix * u_point_light_positions[i] : u_point_light_positions[i];
   for(int i = 0; i < u_spot_light_count; ++i)  v_data_out.SpotLightPositions[i]  = transform_to_tangent_space ? inverse_TBN_matrix * u_spot_light_positions[i]  : u_spot_light_positions[i];
   v_data_out.CameraPosition         = transform_to_tangent_space ? inverse_TBN_matrix * u_camera_position : u_camera_position;
   v_data_out.FragmentPosition       = transform_to_tangent_space ? inverse_TBN_matrix * gl_Position.xyz : gl_Position.xyz;
   v_data_out.FragmentPositionDlight = u_dlight_space_matrix * vec4(gl_Position.xyz, 1.0);

   v_data_out.ViewProjMatrix = u_projection_matrix * u_view_matrix;
}




#shader geometry
#version 460 core

layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

// Global constants
const int Max_Point_Lights = 4;
const int Max_Spot_Lights  = 4;

// Input varying data
in Data
{
   vec3 Normal;
   vec4 Colour;
   vec2 TextureCoordinate;

   vec3 PointLightPositions[Max_Point_Lights];
   vec3 SpotLightPositions[Max_Spot_Lights];
   vec3 CameraPosition;
   vec3 FragmentPosition;
   vec4 FragmentPositionDlight;

   mat4 ViewProjMatrix;
   mat3 TBNMatrix;
} v_data_in[];

// Input varying data
out Data
{
   vec3 Normal;
   vec4 Colour;
   vec2 TextureCoordinate;

   vec3 PointLightPositions[Max_Point_Lights];
   vec3 SpotLightPositions[Max_Spot_Lights];
   vec3 CameraPosition;
   vec3 FragmentPosition;
   vec4 FragmentPositionDlight;

   mat3 TBNMatrix;
} v_data_out;

// Uniforms
uniform int u_point_light_count;
uniform int u_spot_light_count;

void main()
{
   const vec3 tangent0    = (gl_in[1].gl_Position - gl_in[0].gl_Position).xyz;
   const vec3 tangent1    = (gl_in[2].gl_Position - gl_in[0].gl_Position).xyz;
   const vec3 face_normal = normalize(cross(tangent0, tangent1));

   for(int i = 0; i < gl_in.length(); i++)
   {
      gl_Position = v_data_in[i].ViewProjMatrix * gl_in[i].gl_Position;
//      gl_Position = v_data_in[i].ViewProjMatrix * (gl_in[i].gl_Position + 0.5*vec4(face_normal, 0.0));

      v_data_out.Normal = face_normal; // Flat shading
//      v_data_out.Normal = v_data_in[i].Normal; // Smooth shading
      v_data_out.Colour = v_data_in[i].Colour;
      v_data_out.TextureCoordinate = v_data_in[i].TextureCoordinate;

      for(int j = 0; j < u_point_light_count; ++j) v_data_out.PointLightPositions[j] = v_data_in[i].PointLightPositions[j];
      for(int j = 0; j < u_spot_light_count; ++j)  v_data_out.SpotLightPositions[j]  = v_data_in[i].SpotLightPositions[j];
      v_data_out.CameraPosition         = v_data_in[i].CameraPosition;
      v_data_out.FragmentPosition       = v_data_in[i].FragmentPosition;
      v_data_out.FragmentPositionDlight = v_data_in[i].FragmentPositionDlight;

      v_data_out.TBNMatrix = v_data_in[i].TBNMatrix;

      EmitVertex();
   }

   EndPrimitive();
}




#shader fragment
#version 460 core

// Global constants
const int   Max_Point_Lights = 4;
const int   Max_Spot_Lights  = 4;
const float Gamma            = 2.2;
const vec3  Sample_Offset_Directions[20] = { vec3( 1, 1,  1), vec3( 1, -1,  1), vec3(-1, -1,  1), vec3(-1, 1,  1),
                                             vec3( 1, 1, -1), vec3( 1, -1, -1), vec3(-1, -1, -1), vec3(-1, 1, -1),
                                             vec3( 1, 1,  0), vec3( 1, -1,  0), vec3(-1, -1,  0), vec3(-1, 1,  0),
                                             vec3( 1, 0,  1), vec3(-1,  0,  1), vec3( 1,  0, -1), vec3(-1, 0, -1),
                                             vec3( 0, 1,  1), vec3( 0, -1,  1), vec3( 0, -1, -1), vec3( 0, 1, -1) };

// Input varying data
in Data
{
   vec3 Normal;
   vec4 Colour;
   vec2 TextureCoordinate;

   vec3 PointLightPositions[Max_Point_Lights];
   vec3 SpotLightPositions[Max_Spot_Lights];
   vec3 CameraPosition;
   vec3 FragmentPosition;
   vec4 FragmentPositionDlight;

   mat3 TBNMatrix;
} v_data_in;

// Output varying data
layout(location = 0) out vec4 fragment_colour;
layout(location = 1) out vec4 bloom_colour;

// Interface Blocks
struct Light            { vec4 Colour; float AmbientIntensity; float DiffuseIntensity; };
struct DirectionalLight { Light Base; vec3 Direction; sampler2D Shadow; };
struct PointLight       { Light Base; vec3 AttenuationCoefficients; samplerCube Shadow; };
struct SpotLight        { PointLight Point; vec3 Direction; float CosConeAngle; };
struct Material         { float SpecularIntensity; float Smoothness; };
struct Texture          { sampler2D SpecularIntensity; float Smoothness; };

// Light/Material Uniforms
uniform int              u_point_light_count;
uniform int              u_spot_light_count;
uniform DirectionalLight u_direct_light;
uniform PointLight       u_point_lights[Max_Point_Lights];
uniform SpotLight        u_spot_lights[Max_Spot_Lights];
uniform Material         u_material;

// Texture Uniforms
uniform bool      u_use_diffuse_map;
uniform bool      u_use_normal_map;
uniform bool      u_use_displacement_map;
uniform bool      u_use_postprocessor;
uniform float     u_point_light_far_plane;
uniform float     u_displacement_map_scale;
uniform sampler2D u_diffuse_map;
uniform sampler2D u_normal_map;
uniform sampler2D u_displacement_map;

bool UseNormalMap() { return u_use_diffuse_map && u_use_normal_map; }
bool UseHeightMap() { return UseNormalMap() && u_use_displacement_map; }
vec3 CalculateFragmentNormal(const vec2 texture_coordinate) { return UseNormalMap() ? normalize(2.0f * texture(u_normal_map, texture_coordinate).rgb - 1.0f) : v_data_in.Normal; }

float CalculateDirectionalLightShadow(const vec2 texture_coordinate)
{
   const vec3 normal = CalculateFragmentNormal(texture_coordinate);

   vec3 projected_coordinates = v_data_in.FragmentPositionDlight.xyz / v_data_in.FragmentPositionDlight.w;
   projected_coordinates = (projected_coordinates + 1.0f) / 2.0f;

   float actual_depth = projected_coordinates.z;
   float shadow = 0.0f;
   if(actual_depth <= 1.0f)
   {
      vec2 texel_size = 1.0f / textureSize(u_direct_light.Shadow, 0);

      const float bias = max(0.05 * (1.0 - dot(normal, normalize(u_direct_light.Direction))), 0.005);
      for(int x = -1; x <= 1; x++)
         for(int y = -1; y <= 1; y++)
         {
            float pcf_depth = bias + texture(u_direct_light.Shadow, projected_coordinates.xy + vec2(x, y) * texel_size).r;
            shadow += actual_depth > pcf_depth ? 1.0f : 0.0f;
         }
      shadow /= 9.0f;
   }

   return shadow;
}

float CalculatePointLightShadow(const PointLight point_light, const vec3 light_position)
{
   const vec3  light_to_fragment = v_data_in.TBNMatrix * (v_data_in.FragmentPosition - light_position);
   const float current_depth     = length(light_to_fragment);

   const float bias          = 0.15;
   const int   n_samples     = 20;
   const float view_distance = length(v_data_in.FragmentPosition - v_data_in.CameraPosition);
   const float disk_radius   = (1.0 + (view_distance / u_point_light_far_plane)) / 25.0;;
   float shadow = 0.0;
   for(int i = 0; i < n_samples; ++i)
   {
      const float closest_depth =
         bias + u_point_light_far_plane * texture(point_light.Shadow, light_to_fragment + disk_radius * Sample_Offset_Directions[i]).r;
      if(current_depth > closest_depth) shadow += 1.0f;
   }
   shadow /= float(n_samples);

   return shadow;
}

vec4 CalculateLightByDirection(Light light, vec3 direction, float shadow_factor, const vec2 texture_coordinate)
{
   const vec3 normal         = CalculateFragmentNormal(texture_coordinate);
   const vec4 ambient_colour = light.AmbientIntensity * light.Colour;

   // Diffuse lighting
   const float diffuse_factor = max(dot(-normal, normalize(direction)), 0.0f);
   const vec4 diffuse_colour  = diffuse_factor * light.DiffuseIntensity * light.Colour;

   // Specular lighting
   const vec3  fragment_to_camera = normalize(v_data_in.CameraPosition - v_data_in.FragmentPosition);
   const vec3  fragment_to_light  = normalize(-direction);
   const vec3  halfway_direction  = normalize(fragment_to_light + fragment_to_camera);
   const float specular_factor    = pow(max(dot(halfway_direction, normal), 0.0f), u_material.Smoothness);
   const vec4  specular_colour    = specular_factor * u_material.SpecularIntensity * light.Colour;

   return ambient_colour + (1.0f - shadow_factor) * (diffuse_colour + specular_colour);
}

vec4 CalculatePointLight(const PointLight point_light, const vec3 light_position, const vec2 texture_coordinate)
{
   vec3 light_to_fragment = v_data_in.FragmentPosition - light_position;
   float distance         = length(light_to_fragment);
   light_to_fragment      = normalize(light_to_fragment);

   const float attenuation = point_light.AttenuationCoefficients[0] +
                             point_light.AttenuationCoefficients[1] * distance +
                             point_light.AttenuationCoefficients[2] * distance * distance;

   return CalculateLightByDirection(point_light.Base, light_to_fragment, CalculatePointLightShadow(point_light, light_position), texture_coordinate) / attenuation;
}

vec4 CalculateSpotLight(const SpotLight spot_light, const vec3 light_position, const vec2 texture_coordinate)
{
   vec3 ray_direction = normalize(v_data_in.FragmentPosition - light_position);
   float spot_light_factor = dot(ray_direction, spot_light.Direction);

   if(spot_light_factor > spot_light.CosConeAngle)
      return (1.0f - (1.0f - spot_light_factor) / (1.0f - spot_light.CosConeAngle)) * CalculatePointLight(spot_light.Point, light_position, texture_coordinate);
   else return vec4(0.0);
}

vec4 CalculateDirectionalLight(const vec2 texture_coordinate)
{
   return CalculateLightByDirection(u_direct_light.Base, u_direct_light.Direction, CalculateDirectionalLightShadow(texture_coordinate), texture_coordinate);
}

vec4 CalculatePointLights(const vec2 texture_coordinate)
{
   vec4 total_colour = vec4(0.0, 0.0, 0.0, 0.0);
   for(int i = 0; i < u_point_light_count; ++i) total_colour += CalculatePointLight(u_point_lights[i], v_data_in.PointLightPositions[i], texture_coordinate);
   return total_colour;
}

vec4 CalculateSpotLights(const vec2 texture_coordinate)
{
   vec4 total_colour = vec4(0.0, 0.0, 0.0, 0.0);
   for(int i = 0; i < u_spot_light_count; i++) total_colour += CalculateSpotLight(u_spot_lights[i], v_data_in.SpotLightPositions[i], texture_coordinate);
   return total_colour;
}

float CalculateHeight(const vec2 texture_coordinate) { return 1.0 - texture(u_displacement_map, texture_coordinate).r; }

vec2 CalculateParallax()
{
   if(!UseHeightMap()) return v_data_in.TextureCoordinate;

   const vec3 fragment_to_camera = normalize(v_data_in.CameraPosition - v_data_in.FragmentPosition);
   const float min_layers        = 64.0f;
   const float max_layers        = 128.0f;
   const float n_layers          = mix(max_layers, min_layers, max(dot(vec3(0.0, 0.0, 1.0), fragment_to_camera), 0.0));
   const float layer_height      = 1.0 / n_layers;
   const vec2  view_dir          = -fragment_to_camera.xy;
   const vec2  offset            = u_displacement_map_scale * view_dir.xy;
   const vec2  delta_offset      = offset / n_layers;

   // Perform steep parallax mapping
   vec2 current_texture_coordinate = v_data_in.TextureCoordinate;
   float current_height            = CalculateHeight(current_texture_coordinate);
   float current_layer_height      = 0.0;
//   const float sign                = 1.0f;
   const float sign                = -1.0f;
   while(current_layer_height < current_height)
   {
      current_texture_coordinate -= sign * delta_offset;
      current_height              = CalculateHeight(current_texture_coordinate);
      current_layer_height       += layer_height;
   }

   // Perform parallax occlusion mapping
   const vec2  prev_tex_coord  = current_texture_coordinate + sign * delta_offset;
   const float after_height    = current_height - current_layer_height;
   const float before_height   = CalculateHeight(prev_tex_coord) - current_layer_height + layer_height;
   const float weight          = after_height / (after_height - before_height);
   const vec2  final_tex_coord = prev_tex_coord * weight + current_texture_coordinate * (1.0 - weight);

   // Discard fragment if texture coordinate is out of bounds.
   if(final_tex_coord.x > 1.0 || final_tex_coord.y > 1.0 || final_tex_coord.x < 0.0 || final_tex_coord.y < 0.0) discard;

   return final_tex_coord;
}

vec4 GammaCorrect(vec4 colour) { return vec4(pow(colour.rgb, vec3(1.0f / Gamma)), colour.a); }

void main()
{
   const vec2 texture_coord   = CalculateParallax();
   const vec4 material_colour = u_use_diffuse_map ? texture(u_diffuse_map, texture_coord) : v_data_in.Colour;

   // Discard fragment if transparent.
   if(material_colour.a < 0.0001) discard;

   const vec4 lighting_colour = CalculateDirectionalLight(texture_coord) +
                                CalculatePointLights(texture_coord) +
                                CalculateSpotLights(texture_coord);

   // Blend lighting and material colours.
   fragment_colour = vec4(lighting_colour.rgb * material_colour.rgb, material_colour.a);

   // Either gamma-correct or post-process.
   if(!u_use_postprocessor) fragment_colour = GammaCorrect(fragment_colour);
   else
   {
      // If the fragment brightness is higher than the threshold, write to bloom colour buffer.
      float brightness = dot(fragment_colour.rgb, vec3(0.2126, 0.7152, 0.0722));
      bloom_colour     = brightness > 1.0 ? fragment_colour : vec4(0.0, 0.0, 0.0, 1.0);
   }

}