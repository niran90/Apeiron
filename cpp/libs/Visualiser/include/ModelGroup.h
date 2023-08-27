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
#include "ActionBase.h"
#include "Buffers.h"
#include "Colour.h"
#include "Mesh.h"
#include "Material.h"
#include "RenderObject.h"
#include "Model.h"
#include "Texture.h"

#include <optional>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace aprn::vis {

class Shader;

class ModelGroup : public RenderObject
{
 public:
   ModelGroup() = default;

   ModelGroup(const ModelGroup&) = default;

   ModelGroup(ModelGroup&&) noexcept = default;

   ~ModelGroup();

   /** Set Model Group Attributes
   ************************************************************************************************************************************************************/
   ModelGroup& SetName(const std::string& name) override;

   ModelGroup& SetColour(const SVectorR4& rgba_colour) override;

   ModelGroup& SetColour(const Colour& colour) override;

   ModelGroup& SetMaterial(const std::string& name, Real specular_intensity, Real smoothness) override;

   ModelGroup& SetTexture(const std::string& material, size_t index, size_t resolution, Real dispacement_scale) override;

   ModelGroup& SetTexture(const std::string& material, const std::string& item, size_t index, size_t resolution, Real dispacement_scale) override;

   /** Set Model Group Actions
   ************************************************************************************************************************************************************/
   ModelGroup& OffsetPosition(const SVectorR3& displacement) override;

   ModelGroup& OffsetOrientation(Real angle, const SVectorR3& axis) override;

   ModelGroup& Scale(Real factor, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   ModelGroup& Scale(const SVectorR3& factors, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   ModelGroup& MoveBy(const SVectorR3& displacement, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   ModelGroup& MoveTo(const SVectorR3& position, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   ModelGroup& MoveAt(const SVectorR3& velocity, Real start_time = Zero, Reparametriser ramp = Identity) override;

   ModelGroup& Trace(std::function<SVectorR3(Real)> path, Real start_time, Real end_time = InfFloat<>, Reparametriser reparam = Linear) override;

   ModelGroup& RotateBy(Real angle, const SVectorR3& axis, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   ModelGroup& RotateAt(const SVectorR3& angular_velocity, Real start_time = Zero, Reparametriser ramp = Identity) override;

   ModelGroup& RevolveBy(Real angle, const SVectorR3& axis, const SVectorR3& refe_point, Real start_time, Real end_time, Reparametriser reparam = Linear) override;

   ModelGroup& RevolveAt(const SVectorR3& angular_velocity, const SVectorR3& refe_point, Real start_time = Zero, Reparametriser ramp = Identity) override;

   /** Part/sub-model Addition
   ************************************************************************************************************************************************************/
   ModelGroup& Add(const Model& model);

   ModelGroup& Add(const ModelGroup& model_group);

   ModelGroup& Add(Model&& model);

   ModelGroup& Add(ModelGroup&& model_group);

   /** Assignment Operators
   ************************************************************************************************************************************************************/
   ModelGroup& operator=(const ModelGroup&) = default;

   ModelGroup& operator=(ModelGroup&&) noexcept = default;

   /** Other
   ************************************************************************************************************************************************************/
   bool Initialised() const override;

 protected:
   friend class Visualiser;

   void Init() override;

   void ComputeLifespan() override;

   void LoadTextureMap(const std::unordered_map<std::string, Texture&>& texture_map) override;

   void Update(Real global_time) override;

   void Render(Shader& shader) override;

   void Delete() override;

   DArray<SPtr<RenderObject>> SubModels_;
};

}
