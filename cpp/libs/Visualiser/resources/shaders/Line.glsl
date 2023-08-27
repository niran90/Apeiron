#shader vertex
#version 460 core

layout(std430, binding = 0) buffer TVertex { vec4 vertex[]; };

uniform mat4  u_mvp_matrix; // proj * view * model
uniform vec2  u_resolution;
uniform float u_thickness;

void main()
{
   /**
      0 3   5
      o o---o
      |\ \  |
      | \ \ |
      |  \ \|
      o---o o
      1   2 4
   */
   int line_id = gl_VertexID / 6; // Line segment index (each segment has 6 vertices).
   int vert_id = gl_VertexID % 6; // Vertex index in each triangle (indices 0 - 5).

   // The coordinates have to be transformed from model space to window space.
   vec4 vertices[4]; // Four vertices belonging to a string of 3 line segments.
   for(int i = 0; i < 4; ++i)
   {
       vertices[i]      = u_mvp_matrix * vertex[line_id + i];
       vertices[i].xyz /= vertices[i].w; // Perspective divide
       vertices[i].xy   = (vertices[i].xy + 1.0) * 0.5 * u_resolution;
   }

   //   0            1            2            3
   //   o------------o------------o------------o
   //        prev         curr         next
   const vec2 tangent_curr = normalize(vertices[2].xy - vertices[1].xy);
   const vec2 normal_curr  = vec2(-tangent_curr.y, tangent_curr.x);

   vec4 position;
   if(vert_id == 0 || vert_id == 1 || vert_id == 3)
   {
      const vec2  tangent_prev = normalize(vertices[1].xy - vertices[0].xy);
      const vec2  normal_prev  = vec2(-tangent_prev.y, tangent_prev.x);
      const vec2  miter_1      = normalize(normal_curr + normal_prev);
      const float offset       = (vert_id == 1 ? -0.5 : 0.5) * u_thickness / dot(miter_1, normal_curr);

      position     = vertices[1];
      position.xy += offset * miter_1;
   }
   else // vert_id == 2 || vert_id == 4 || vert_id == 5
   {
      const vec2  tangent_next = normalize(vertices[3].xy - vertices[2].xy);
      const vec2  normal_next  = vec2(-tangent_next.y, tangent_next.x);
      const vec2  miter_2      = normalize(normal_curr + normal_next);
      const float offset       = (vert_id == 5 ? 0.5 : -0.5) * u_thickness / dot(miter_2, normal_curr);

      position     = vertices[2];
      position.xy += offset * miter_2;
   }

   // Finally the window coordinates must be transformed back to clip space. Transform from window space to normalized device space.
   // The perspective divide must also be reversed.
   position.xy   = position.xy / u_resolution * 2.0 - 1.0;
   position.xyz *= position.w;

   gl_Position = position;
}



#shader fragment
#version 460 core

layout(location = 0) out vec4 fragment_colour;

void main()
{
   fragment_colour = vec4(1.0);
}