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

#include "../include/Global.h"
//#include "Visualiser/include/Visualiser.h"
//#include "Visualiser/include/Scene.h"
//
//#include "FileManager/include/File.h"
//#include "FileManager/include/FileSystem.h"
//
//#include <string>
//#include <string_view>
//
using namespace aprn;

#include <openvdb/openvdb.h>
#include <openvdb/math/FiniteDifference.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/LevelSetAdvect.h>
#include <openvdb/tools/ValueTransformer.h>
#include <execution>

using namespace openvdb;
using namespace openvdb::math;

void makeCylinder(FloatGrid::Ptr grid, float radius, const openvdb::Vec3d& centre, const CoordBBox& indexBB, double h, float background_value)
{
   FloatGrid::Accessor accessor = grid->getAccessor();

   const auto min = indexBB.min();
   const auto max = indexBB.max();

   for (Int32 i = min.x(); i <= max.x(); ++i) {
      for (Int32 j = min.y(); j <= max.y(); ++j) {
         for (Int32 k = min.z(); k <= max.z(); ++k) {
            // transform point (i, j, k) of index space into world space
            Vec3d p(i * h, j * h, k * h);
            // compute level set function value
            double dx = p.x() - centre.x();
            double dy = p.y() - centre.y();
            float distance = sqrt(dx*dx + dy*dy) - radius;

            if(abs(distance) < background_value) accessor.setValue(Coord(i, j, k), distance);
         }
      }
   }

   grid->setTransform(openvdb::math::Transform::createLinearTransform(h));

//   openvdb::v10_0::tools::extractIsosurfaceMask();
}

float SamplePoint(FloatGrid::Ptr grid, openvdb::Vec3f global_point)
{
   FloatGrid::ConstAccessor accessor = grid->getConstAccessor();
   openvdb::tools::GridSampler<FloatGrid::ConstAccessor, openvdb::tools::QuadraticSampler> sampler(accessor, grid->transform());
   return sampler.wsSample(global_point);
}

void createAndSaveCylinder()
{
   openvdb::initialize();

   float background_value = 0.7;
   openvdb::FloatGrid::Ptr grid0 = openvdb::FloatGrid::create(background_value);
   openvdb::FloatGrid::Ptr grid1 = openvdb::FloatGrid::create(background_value);

   // Common attributes.
//   const double h  = 0.2;
//   CoordBBox indexBB(Coord(-40, -40, -40), Coord(40, 40, 40));
   const double h  = 0.1;
   CoordBBox indexBB(Coord(-80, -80, -5), Coord(80, 80, 5));

   // Make cylinder 0.
   const float r0 = 2.5f;
   const Vec3d c0 = {2.75, 0.0, 0.0};
   makeCylinder(grid0, r0, c0, indexBB, h, background_value);
   grid0->setName("LevelSetCylinder0");
   grid0->setGridClass(openvdb::GRID_LEVEL_SET);

   const openvdb::Vec3f s0 = {-5.0, 0.0, 0.0};
   const openvdb::Vec3f s1 = {-2.5, 0.0, 0.0};
   const openvdb::Vec3f s2 = { 0.0, 0.0, 0.0};
   const openvdb::Vec3f s3 = { 2.5, 0.0, 0.0};
   const openvdb::Vec3f s4 = { 5.0, 0.0, 0.0};
   Print("Grid0 Values:", SamplePoint(grid0, s0), SamplePoint(grid0, s1), SamplePoint(grid0, s2), SamplePoint(grid0, s3), SamplePoint(grid0, s4));

   // Make cylinder 1.
   const float r1 = 2.5f;
   const Vec3d c1 = {-2.75, 0.0, 0.0};
   makeCylinder(grid1, r1, c1, indexBB, h, background_value);
   grid1->setName("LevelSetCylinder1");
   grid1->setGridClass(openvdb::GRID_LEVEL_SET);
   Print("Grid1 Values:", SamplePoint(grid1, s0), SamplePoint(grid1, s1), SamplePoint(grid1, s2), SamplePoint(grid1, s3), SamplePoint(grid1, s4));

   openvdb::tools::csgUnion(*grid0, *grid1);
//   openvdb::tools::csgIntersection(*grid0, *grid1);
//   openvdb::tools::csgDifference(*grid0, *grid1);
   Print("Grid Combination Values:", SamplePoint(grid0, s0), SamplePoint(grid0, s1), SamplePoint(grid0, s2), SamplePoint(grid0, s3), SamplePoint(grid0, s4));

//   // Define a functor that offsets the levelset.
//   const double offset = -1.0;
//   auto func = [&offset](const openvdb::FloatGrid::ValueAllIter& iter){ iter.setValue(iter.getValue() - offset); };
//   openvdb::tools::foreach(grid0->beginValueAll(), func);

   // Compute the gradient grid.
   const VectorGrid::Ptr gradient_grid = openvdb::tools::gradient(*grid0);

   // Advect the level-set field.
   openvdb::tools::DiscreteField<VectorGrid, openvdb::tools::QuadraticSampler> velocity_field(*gradient_grid);
   openvdb::tools::LevelSetAdvection levelset_advector(*grid0, velocity_field);
   levelset_advector.setSpatialScheme(openvdb::math::HJWENO5_BIAS);
   levelset_advector.setTemporalScheme(openvdb::math::TVD_RK2);
   levelset_advector.setTrackerSpatialScheme(openvdb::math::HJWENO5_BIAS);
   levelset_advector.setTrackerTemporalScheme(openvdb::math::TVD_RK1);

   const float dT = 0.7;
   const float dt = 0.01;
   uint iter{};
   for(float t = 0; t < dT; t += dt) Print("Iter:", iter++, "\tTime:", t, "\tSteps:", levelset_advector.advect(t, t + dt));

   // Save grid to file
   openvdb::io::File file("mygrids.vdb");
   openvdb::GridPtrVec grids;
   grids.push_back(grid0);
   file.write(grids);
   file.close();
}

int main()
{
   createAndSaveCylinder();

   return 0;
}