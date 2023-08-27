/***************************************************************************************************************************************************************
* GPL-3.0 License
* Copyright (C) 2022 Niran A. Ilangakoon
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
* published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program.
* If not, see <https://www.gnu.org/licenses/>.
***************************************************************************************************************************************************************/

#include "../include/Model.h"
#include "../include/Action.h"
#include "../include/Animator.h"
#include "../include/Shader.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Model Public Interface
***************************************************************************************************************************************************************/
Model::Model(const Model& other)
{
   ASSERT(!Init_, "Cannot yet copy construct a model if it has already been initialised.")
   *this = other;
}

Model::Model(Model&& other) noexcept
{
   ASSERT(!Init_, "Cannot yet move construct a model if it has already been initialised.")
   *this = std::move(other);
}

Model::~Model() { Delete(); }

/** Set Model Attributes
***************************************************************************************************************************************************************/
Model&
Model::SetName(const std::string& name)
{
   Name_ = name;
   return *this;
}

Model&
Model::SetColour(const SVectorR4& rgba_colour) { return SetColour(Colour{rgba_colour}); }

Model&
Model::SetColour(const Colour& colour)
{
   FOR_EACH(vertex, Mesh_.Vertices_) vertex.Colour = SVectorToGlmVec(colour.Values);
   return *this;
}

Model&
Model::SetMaterial(const std::string& name, const Real specular_intensity, const Real smoothness)
{
   Material_.emplace(name, specular_intensity, smoothness);
   return *this;
}

Model&
Model::SetTexture(const std::string& material, const size_t index, const size_t resolution, const Real dispacement_scale)
{
   SetTexture(material, "", index, resolution, dispacement_scale);
   return *this;
}

Model&
Model::SetTexture(const std::string& material, const std::string& item, const size_t index, const size_t resolution, const Real dispacement_scale)
{
   TextureRequest_.emplace(std::make_pair(TextureName(material, item, index, resolution), dispacement_scale));
   return *this;
}

/** Set Model Actions
***************************************************************************************************************************************************************/
Model&
Model::OffsetPosition(const SVectorR3& displacement)
{
   SPtr<ActionBase> ptr = std::make_shared<Action<ActionType::OffsetPosition>>(Animator_, SVectorToGlmVec(displacement));
   Animator_.Add(ActionType::OffsetPosition, ptr);
   return *this;
}

Model&
Model::OffsetOrientation(const Real angle, const SVectorR3& axis)
{
   SPtr<ActionBase> ptr = std::make_shared<Action<ActionType::OffsetOrientation>>(Animator_, angle, SVectorToGlmVec(axis));
   Animator_.Add(ActionType::OffsetOrientation, ptr);
   return *this;
}

Model&
Model::Scale(const Real factor, const Real start_time, const Real end_time, Reparametriser reparam)
{
   Scale(SVectorR3(factor), start_time, end_time, reparam);
   return *this;
}

Model&
Model::Scale(const SVectorR3& factors, const Real start_time, const Real end_time, Reparametriser reparam)
{
   SPtr<ActionBase> ptr = std::make_shared<Action<ActionType::Scale>>(Animator_, SVectorToGlmVec(factors), start_time, end_time);
   Animator_.Add(ActionType::Scale, ptr);
   return *this;
}

Model&
Model::MoveBy(const SVectorR3& displacement, const Real start_time, const Real end_time, Reparametriser reparam)
{
   SPtr<ActionBase> ptr = std::make_shared<Action<ActionType::MoveBy>>(Animator_, SVectorToGlmVec(displacement), start_time, end_time);
   Animator_.Add(ActionType::MoveBy, ptr);
   return *this;
}

Model&
Model::MoveTo(const SVectorR3& position, const Real start_time, const Real end_time, Reparametriser reparam)
{
   SPtr<ActionBase> ptr = std::make_shared<Action<ActionType::MoveTo>>(Animator_, SVectorToGlmVec(position), start_time, end_time);
   Animator_.Add(ActionType::MoveTo, ptr);
   return *this;
}

Model&
Model::MoveAt(const SVectorR3& velocity, const Real start_time, Reparametriser ramp)
{
   EXIT("TODO")
   return *this;
}

Model&
Model::Trace(std::function<SVectorR3(Real)> path, const Real start_time, const Real end_time, Reparametriser reparam)
{
   SPtr<ActionBase> ptr = std::make_shared<Action<ActionType::Trace>>(Animator_, path, start_time, end_time);
   Animator_.Add(ActionType::Trace, ptr);
   return *this;
}

Model&
Model::RotateBy(const Real angle, const SVectorR3& axis, const Real start_time, const Real end_time, Reparametriser reparam)
{
   SPtr<ActionBase> ptr = std::make_shared<Action<ActionType::RotateBy>>(Animator_, angle, SVectorToGlmVec(axis), start_time, end_time);
   Animator_.Add(ActionType::RotateBy, ptr);
   return *this;
}

Model&
Model::RotateAt(const SVectorR3& angular_velocity, const Real start_time, Reparametriser ramp)
{
   SPtr<ActionBase> ptr = std::make_shared<Action<ActionType::RotateAt>>(Animator_, SVectorToGlmVec(angular_velocity), start_time, ramp);
   Animator_.Add(ActionType::RotateAt, ptr);
   return *this;
}

Model&
Model::RevolveBy(const Real angle, const SVectorR3& axis, const SVectorR3& refe_point, const Real start_time, const Real end_time,
                 Reparametriser reparam)
{
   SPtr<ActionBase> ptr =
      std::make_shared<Action<ActionType::RevolveBy>>(Animator_, angle, SVectorToGlmVec(axis), SVectorToGlmVec(refe_point), start_time, end_time);
   Animator_.Add(ActionType::RevolveBy, ptr);
   return *this;
}

Model&
Model::RevolveAt(const SVectorR3& angular_velocity, const SVectorR3& refe_point, const Real start_time, Reparametriser ramp)
{
   EXIT("TODO")
   return *this;
}

/** Assignment Operators
***************************************************************************************************************************************************************/
Model&
Model::operator=(const Model& model)
{
   ASSERT(!Init_,       "Cannot yet copy assign to a model if it has already been initialised.")
   ASSERT(!model.Init_, "Cannot yet copy assign from a model if it has already been initialised.")

   // NOTE: Should NOT overwrite the original buffer IDs of VAO, VBO, and EBO.
   Mesh_           = model.Mesh_;
   Animator_       = model.Animator_;
   TextureRequest_ = model.TextureRequest_;
   Material_       = model.Material_;
   Centroid_       = model.Centroid_;
   StrokeColour_   = model.StrokeColour_;
   FillColour_     = model.FillColour_;

   return *this;
}

Model&
Model::operator=(Model&& model) noexcept
{
   ASSERT(!Init_,       "Cannot yet move assign to a model if it has already been initialised.")
   ASSERT(!model.Init_, "Cannot yet move assign from a model if it has already been initialised.")

   // NOTE: Should NOT overwrite the original buffer IDs of VAO, VBO, and EBO.
   Mesh_           = std::move(model.Mesh_);
   Animator_       = std::move(model.Animator_);
   TextureRequest_ = std::move(model.TextureRequest_);
   Material_       = std::move(model.Material_);
   Centroid_       = std::move(model.Centroid_);
   StrokeColour_   = std::move(model.StrokeColour_);
   FillColour_     = std::move(model.FillColour_);

   // Reset moved-from model as it is now in an undefined state. Note: to avoid an infinite regress, we need to specifically invoke the copy assigment operator
   // here, NOT the move assignment operator.
   model = Unmove(Model());

   return *this;
}

/***************************************************************************************************************************************************************
* Model Protected Interface
***************************************************************************************************************************************************************/
void
Model::Init()
{
   if(Init_) return;
   ASSERT(Mesh_.Loaded(), "The model cannot be initialised without a mesh.")

   // Compute entry/exit times.
   ComputeLifespan();

   // Compute vertex normals.
   Mesh_.ComputeVertexNormals();

   // Initialise VAO, VBO, and EBO.
   VAO_.Init();
   VBO_.Init(Mesh_.Vertices_);
   EBO_.Init(Mesh_.Indices_);

   // Add vertex buffer to vertex array object.
   VAO_.Bind();
   VAO_.AddBuffer(VBO_, Mesh_.VertexLayout());
   VAO_.Unbind();

   Init_ = true;
}

void
Model::ComputeLifespan() { std::tie(EntryTime_, ExitTime_) = Animator_.ComputeLifespan(); }

void
Model::LoadTextureMap(const std::unordered_map<std::string, Texture&>& texture_map)
{
   Textures_ = texture_map;
   TextureRequest_.reset(); // Can nullify as the request is no longer required.
}

void
Model::Update(const Real global_time)
{
   if(!Init_) return;

   // Direct the animator to update the model.
   Animator_.Update(global_time);

   // Update the vertex buffer if the mesh has been modified.
   VBO_.Update(Mesh_.Vertices_);
}

void
Model::Render(Shader& shader)
{
   if(!Init_) return;

   constexpr int slot_offset(3); // TODO - currently hard-coded.

   if(Material_) shader.UseMaterial(Material_.value());
   if(!Textures_.empty())
   {
      size_t texture_index = 0;
      FOR_EACH(type_string, texture, Textures_)
      {
         // Configure respective texture uniform.
         const auto& uniform_name = TextureUniformString(type_string);
         shader.UseTexture(texture, "u_" + uniform_name, slot_offset + texture_index++);
         shader.SetUniform1i("u_use_" + uniform_name, 1);

         // Set scale if this is a displacement map.
         if(GetTextureType(type_string) == TextureType::Displacement)
         {
            const auto scale = texture.MapScale();
            shader.SetUniform1f("u_" + uniform_name + "_scale", scale);
         }
      }
   }

   shader.UseModel(*this);
   DrawElements();

   // Switch off texture maps and unbind texture.
   if(!Textures_.empty())
      FOR_EACH(type_string, texture, Textures_)
      {
         shader.SetUniform1i("u_use_" + TextureUniformString(type_string), 0);
         texture.Unbind();
      }
}

void
Model::DrawElements() const
{
   if(!Init_) return;

   // Model main body
   VAO_.Bind();
   EBO_.Bind();
   GLCall(glDrawElements(GL_TRIANGLES, EBO_.IndexCount(), GL_UNSIGNED_INT, nullptr))
   EBO_.Unbind();
   VAO_.Unbind();
}

void
Model::Delete()
{
   VBO_.Delete();
   VAO_.Delete();
   EBO_.Delete();
}

}
