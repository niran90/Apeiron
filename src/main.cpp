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

//using namespace aprn;

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
#include <openvdb/tools/FastSweeping.h>
#include <openvdb/tools/Mask.h>

using namespace openvdb;
using namespace openvdb::math;

void makeSphere(FloatGrid::Ptr grid, float radius, const openvdb::Vec3d& centre, const CoordBBox& indexBB, double h, float background_value)
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
            double dz = p.z() - centre.z();
            float distance = sqrt(dx*dx + dy*dy + dz*dz) - radius;

            if(abs(distance) < background_value) accessor.setValue(Coord(i, j, k), distance);
         }
      }
   }

   grid->setTransform(openvdb::math::Transform::createLinearTransform(h));
}

void makeCylinder(FloatGrid::Ptr grid, float half_height, float radius, const openvdb::Vec3d& centre, const CoordBBox& indexBB, double h, float background_value, bool invert = false)
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
            float l = abs(p.z() - centre.z()) - half_height;
            float r = sqrt(dx * dx + dy * dy) - radius;
            float distance{};

            if(r >= 0.0f && l >= 0.0f) distance = sqrt(r * r + l * l);
            else if(r > l) distance = r;
            else distance = l;

            if(abs(distance) < background_value) accessor.setValue(Coord(i, j, k), (invert ? -1.0 : 1.0) * abs(distance));
         }
      }
   }

   grid->setTransform(openvdb::math::Transform::createLinearTransform(h));
}

void createAndSaveSphere()
{
   openvdb::initialize();

   float background_value = 1.2;
   openvdb::FloatGrid::Ptr grid0 = openvdb::FloatGrid::create(background_value);
   openvdb::FloatGrid::Ptr grid1 = openvdb::FloatGrid::create(background_value);

   // Common attributes.
   const double h  = 0.05;
   CoordBBox indexBB(Coord(-160, -160, -100), Coord(160, 160, 100));

   // Make cylinder 0.
   const float r0 = 2.5f;
   const Vec3d c0 = {2.75, 0.0, 0.0};
   makeSphere(grid0, r0, c0, indexBB, h, background_value);
   grid0->setName("LevelSetCylinder0");
   grid0->setGridClass(openvdb::GRID_LEVEL_SET);

   // Make cylinder 1.
   const float r1 = 2.5f;
   const Vec3d c1 = {-2.75, 0.0, 0.0};
   makeSphere(grid1, r1, c1, indexBB, h, background_value);
   grid1->setName("LevelSetCylinder1");
   grid1->setGridClass(openvdb::GRID_LEVEL_SET);

   openvdb::tools::csgUnion(*grid0, *grid1);
   const auto grid_init = grid0->deepCopy();

   // Define a functor that offsets the levelset.
   double offset = 1.0;
   auto func0 = [&offset](const openvdb::FloatGrid::ValueAllIter& iter){ iter.setValue(iter.getValue() - offset); };
   openvdb::tools::foreach(grid0->beginValueAll(), func0);

   grid0 = openvdb::tools::sdfToSdf(*grid0, 0.0, 1);
   offset *= -1.0;
   openvdb::tools::foreach(grid0->beginValueAll(), func0);

   openvdb::tools::csgDifference(*grid0, *grid_init);
   grid0->tree().prune();

   const auto mask = openvdb::tools::sdfInteriorMask(*grid0);

   // Save grid to file
   openvdb::io::File file("mygrids.vdb");
   openvdb::GridPtrVec grids;
   grids.push_back(mask);
   file.write(grids);
   file.close();
}

void createAndSaveCylinder()
{
   openvdb::initialize();

   const float background_value = 1.2;
   openvdb::FloatGrid::Ptr grid0 = openvdb::FloatGrid::create(background_value);
   openvdb::FloatGrid::Ptr grid1 = openvdb::FloatGrid::create(background_value);
   openvdb::FloatGrid::Ptr grid2 = openvdb::FloatGrid::create(background_value);

   // Common attributes.
//   const double h  = 0.1;
//   CoordBBox indexBB(Coord(-80, -80, -80), Coord(80, 80, 80));
   const double h  = 0.05;
   CoordBBox indexBB(Coord(-160, -160, -160), Coord(160, 160, 160));

   // Make cylinder 0.
   const float r0 = 2.5f;
   const float h0 = 2.5f;
   const Vec3d c0 = {2.75, 0.0, 0.0};
   makeCylinder(grid0, h0, r0, c0, indexBB, h, background_value);
   grid0->setName("LevelSetCylinder0");

   // Make cylinder 1.
   const float r1 = 2.5f;
   const float h1 = 2.5f;
   const Vec3d c1 = {-2.75, 0.0, 0.0};
   makeCylinder(grid1, h1, r1, c1, indexBB, h, background_value);
   grid1->setName("LevelSetCylinder1");

   // Make cylinder 2.
   const float r2 = 6.25f;
   const float h2 = 3.5f;
   const Vec3d c2 = {0.0, 0.0, 0.0};
   makeCylinder(grid2, h2, r2, c2, indexBB, h, background_value);
   grid2->setName("LevelSetCylinder2");

   openvdb::tools::csgUnion(*grid0, *grid1);
   openvdb::tools::csgUnion(*grid0, *grid2);
   const auto grid_init = grid0->deepCopy();

   // Define a functor that offsets the levelset.
   double width = 0.5;
//   double alpha = 0.5;
//   double alpha = 0.55;
   double alpha = 0.6;
//   double alpha = 0.75;
//   double alpha = 1.0;
   double offset = alpha * width;
   auto func0 = [&offset](const openvdb::FloatGrid::ValueAllIter& iter){ iter.setValue(iter.getValue() - offset); };
   openvdb::tools::foreach(grid0->beginValueAll(), func0);

   grid0 = openvdb::tools::sdfToSdf(*grid0, 0.0, 1);
   offset *= -1.0;
   openvdb::tools::foreach(grid0->beginValueAll(), func0);

   openvdb::tools::csgDifference(*grid0, *grid_init);

   offset *= 0.5*(sqrt(2.0) - 1.0);
   openvdb::tools::foreach(grid0->beginValueAll(), func0);
   grid0 = openvdb::tools::sdfToSdf(*grid0, 0.0, 1);

   offset *= -1.0;
   openvdb::tools::foreach(grid0->beginValueAll(), func0);
   grid0 = openvdb::tools::sdfToSdf(*grid0, 0.0, 1);

   openvdb::tools::csgDifference(*grid0, *grid_init);
   grid0 = openvdb::tools::sdfToSdf(*grid0, 0.0, 1);
   auto func1 = [&background_value](const openvdb::FloatGrid::ValueAllIter& iter){ if(background_value < abs(iter.getValue())) iter.setValue(background_value); };
   openvdb::tools::foreach(grid0->beginValueAll(), func1);
   grid0->pruneGrid();

   const auto mask = openvdb::tools::sdfInteriorMask(*grid0);

   // Save grid to file
   openvdb::io::File file("mygrids.vdb");
   openvdb::GridPtrVec grids;
   grids.push_back(grid0);
//   grids.push_back(mask);
   file.write(grids);
   file.close();
}

void processSgtBluff()
{
   openvdb::initialize();

   // Read all grids from a file.
   openvdb::io::File file("sgt_bluff.vdb");
//   openvdb::io::File file("mygrids.vdb");
   file.open();
   openvdb::GridPtrVecPtr myGrids = file.getGrids();
   file.close();

   for(auto& iter : *myGrids)
   {
      auto grid = dynamic_cast<openvdb::FloatGrid*>(iter.get());
      auto grid0 = grid->deepCopy();
      auto grid_init = grid0->deepCopy();

      // Define a functor that offsets the levelset.
      double width = 1.0;
      double alpha = 0.55;
      double offset = alpha * width;
      auto func0 = [&offset](const openvdb::FloatGrid::ValueAllIter& iter){ iter.setValue(iter.getValue() - offset); };
      openvdb::tools::foreach(grid0->beginValueAll(), func0);

      grid0 = openvdb::tools::sdfToSdf(*grid0, 0.0, 1);
      offset *= -1.0;
      openvdb::tools::foreach(grid0->beginValueAll(), func0);

      openvdb::tools::csgDifference(*grid0, *grid_init);

      offset *= 0.5*(sqrt(2.0) - 1.0);
      openvdb::tools::foreach(grid0->beginValueAll(), func0);
      grid0 = openvdb::tools::sdfToSdf(*grid0, 0.0, 1);

      offset *= -1.0;
      openvdb::tools::foreach(grid0->beginValueAll(), func0);
      grid0 = openvdb::tools::sdfToSdf(*grid0, 0.0, 1);

      openvdb::tools::csgDifference(*grid0, *grid_init);
      grid0 = openvdb::tools::sdfToSdf(*grid0, 0.0, 1);
      grid0->pruneGrid();

      offset = 0.2;
      openvdb::tools::foreach(grid_init->beginValueAll(), func0);

      // Write to file.
      openvdb::io::File out("mygrids.vdb");
      openvdb::GridPtrVec grids;
      grids.push_back(grid_init);
      grids.push_back(grid0);
      out.write(grids);
      out.close();
   }
}

int main()
{
//   createAndSaveCylinder();
   processSgtBluff();

   return 0;
}