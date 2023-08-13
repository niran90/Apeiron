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

#pragma once

#include "../../../include/Global.h"
#include "DataContainer/include/Array.h"
#include "Buffers.h"
#include "Colour.h"
#include "Material.h"
#include "Mesh.h"
#include "Object.h"
#include "Texture.h"

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <unordered_map>

namespace aprn::vis {

class Model : public Object
{
 public:
   Model() = default;

   Model(const Model& other);

   Model(Model&& other) noexcept;

   ~Model();

   /** Set Model Attributes
   ************************************************************************************************************************************************************/
   Model* SetName(const std::string& name) override;

   Model* SetColour(const SVectorR4& rgba_colour) override;

   Model* SetColour(const Colour& colour) override;

   Model* SetMaterial(const std::string& name, Real specular_intensity, Real smoothness) override;

   Model* SetTexture(const std::string& material, size_t index, size_t resolution, Real dispacement_scale) override;

   Model* SetTexture(const std::string& material, const std::string& item, size_t index, size_t resolution, Real dispacement_scale) override;

   /** Set Model Actions
   ************************************************************************************************************************************************************/
   Model* OffsetPosition(const SVectorR3& displacement) override;

   Model* OffsetOrientation(Real angle, const SVectorR3& axis) override;

   Model* Scale(Real factor, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   Model* Scale(const SVectorR3& factors, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   Model* MoveBy(const SVectorR3& displacement, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   Model* MoveTo(const SVectorR3& position, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   Model* MoveAt(const SVectorR3& velocity, Real start_time = Zero, Reparametriser ramp = Identity) override;

   Model* Trace(std::function<SVectorR3(Real)> path, Real start_time, Real end_time = InfFloat<>, Reparametriser reparam = Linear) override;

   Model* RotateBy(Real angle, const SVectorR3& axis, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   Model* RotateAt(const SVectorR3& angular_velocity, Real start_time = Zero, Reparametriser ramp = Identity) override;

   Model* RevolveBy(Real angle, const SVectorR3& axis, const SVectorR3& refe_point, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   Model* RevolveAt(const SVectorR3& angular_velocity, const SVectorR3& refe_point, Real start_time = Zero, Reparametriser ramp = Identity) override;

   /** Assignment Operators
   ************************************************************************************************************************************************************/
   Model& operator=(const Model& model);

   Model& operator=(Model&& model) noexcept;

   /** Other
   ************************************************************************************************************************************************************/
   inline bool Initialised() const override { return Init_; }

   inline const auto& ModelMatrix() const { return Animator_.ModelMatrix(); }

   inline const auto& ModelMesh() const { return Mesh_; }

   inline const auto& TextureRequest() const { return TextureRequest_; }

 protected:
   friend class Animator;
   friend class ObjectFactory;
   friend class PostProcessor;
   friend class Visualiser;

   void Init() override;

   void ComputeLifespan() override;

   void LoadTextureMap(const std::unordered_map<std::string, Texture&>& texture_map) override;

   void Update(Real global_time) override;

   void Render(Shader& shader) override;

   void DrawElements() const;

   void Delete() override;

   template<class T> using UMap = std::unordered_map<std::string, T>;

   Animator                        Animator_;
   Mesh                            Mesh_;
   VertexArray                     VAO_;
   VertexBuffer                    VBO_;
   IndexBuffer                     EBO_;
   Option<Material>                Material_;
   Option<Pair<std::string, Real>> TextureRequest_; // [texture name, displacement map scale]
   UMap<Texture&>                  Textures_;       // Textures (diffuse, height, normal, etc.) used by this model.
   Colour                          FillColour_;
   glm::vec3                       Centroid_;
};

}
