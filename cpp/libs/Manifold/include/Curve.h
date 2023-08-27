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

#include "LinearAlgebra/include/Vector.h"

namespace aprn::mnfld {

/***************************************************************************************************************************************************************
* Curve Class Definition
***************************************************************************************************************************************************************/
template<size_t ambient_dim = 3>
class Curve
{
   using Vector = SVectorR<ambient_dim>;

 public:
   virtual constexpr Vector Point(const Real param) const = 0;

   virtual constexpr Vector Tangent(const Real param) const = 0;

   virtual constexpr Vector Normal(const Real param) const = 0;

   virtual constexpr Real Length() const = 0;

   constexpr Vector Binormal(const Vector& tangent, const Vector& normal) const;

   constexpr void MakeUnitSpeed() noexcept { UnitSpeed_ = true; }

 protected:
   bool UnitSpeed_{false};
};

/***************************************************************************************************************************************************************
* Linear/Piecewise Linear Curves
***************************************************************************************************************************************************************/

/** Line
***************************************************************************************************************************************************************/
template<size_t ambient_dim = 2>
class Line : public Curve<ambient_dim>
{
   using Vector = SVectorR<ambient_dim>;

 public:
   constexpr Line(const Vector& direction, const Vector& point = Vector{});

   constexpr Vector Point(const Real t) const override;

   constexpr Vector Tangent(const Real t) const override;

   constexpr Vector Normal(const Real t) const override;

   constexpr Real Length() const override { return InfFloat<>; }

 protected:
   Vector Direction;
   Vector Start;
   Real   DirectionNorm_;
   Real   Normaliser_;
};

/** Ray
***************************************************************************************************************************************************************/
template<size_t ambient_dim = 2>
class Ray final : public Line<ambient_dim>
{
   using Vector = SVectorR<ambient_dim>;

 public:
   constexpr Ray(const Vector& direction, const Vector& start = Vector{});

   constexpr Vector Point(const Real t) const override;

   constexpr Real Length() const override { return InfFloat<>; }
};

/** Line Segment
***************************************************************************************************************************************************************/
template<size_t ambient_dim = 2>
class LineSegment final : public Line<ambient_dim>
{
   using Vector = SVectorR<ambient_dim>;

 public:
   constexpr LineSegment(const Vector& start, const Vector& end);

   constexpr Vector Point(const Real t) const override;

   constexpr Real Length() const override { return this->DirectionNorm_; }
};

/** Line Segment Chain
***************************************************************************************************************************************************************/
template<size_t ambient_dim = 2>
class LineSegmentChain final : public Curve<ambient_dim>
{
   using Vector  = SVectorR<ambient_dim>;
   using Segment = LineSegment<ambient_dim>;

 public:
   template<class D>
   LineSegmentChain(const Array<Vector, D>& vertices, bool is_closed = false);

   constexpr Vector Point(const Real l) const override;

   constexpr Vector Tangent(const Real t) const override;

   constexpr Vector Normal(const Real t) const override;

   constexpr Real Length() const override { return ChainLength_; }

 private:
   DArray<Segment> Segments_;
   DArray<Real>    CumulativeLengths_;
   Real            ChainLength_;
   bool            Closed_;
};

/***************************************************************************************************************************************************************
* Elliptical/Circular Curves
***************************************************************************************************************************************************************/

/** Circle
***************************************************************************************************************************************************************/
template<size_t ambient_dim = 2>
class Circle : public Curve<ambient_dim>
{
   using Vector = SVectorR<ambient_dim>;

 public:
   Circle(const Real radius, const Vector& centre = Vector{});

   Circle(const Real radius, const Real start_angle, const Vector& centre = Vector{});

   constexpr Vector Point(const Real t) const override;

   constexpr Vector Tangent(const Real t) const override;

   constexpr Vector Normal(const Real t) const override;

   constexpr Real Length() const override { return Length_; }

 protected:
   constexpr Real Angle(const Real t) const;

   Vector Centre_;
   Real   Radius_;
   Real   StartAngle_{Zero};
   Real   Normaliser_;
   Real   Length_;
};

/** Circular Arc
***************************************************************************************************************************************************************/
template<size_t ambient_dim = 2>
class Arc final : public Circle<ambient_dim>
{
   using Vector = SVectorR<ambient_dim>;

 public:
   Arc(const Real radius, const Real angle, const Vector& centre = Vector{});

   Arc(const Real radius, const Real start_angle, const Real end_angle, const Vector& centre = Vector{});

   constexpr Vector Point(const Real t) const override;

   constexpr Vector Tangent(const Real t) const override;

   constexpr Vector Normal(const Real t) const override;

   constexpr void CheckAngle(const Real t) const;

 private:
   Real EndAngle_;
};

/** Ellipse
***************************************************************************************************************************************************************/
template<size_t ambient_dim = 2>
class Ellipse final : public Curve<ambient_dim>
{
   using Vector = SVectorR<ambient_dim>;

 public:
   Ellipse(const Real radius_x, const Real radius_y, const Vector& centre = Vector{});

   constexpr Vector Point(const Real t) const override;

   constexpr Vector Tangent(const Real t) const override;

   constexpr Vector Normal(const Real t) const override;

   constexpr Real Length() const override { return Length_; }

 private:
   Vector Centre_;
   Real   RadiusX_;
   Real   RadiusY_;
   Real   Length_;

};

/***************************************************************************************************************************************************************
* Polynomial Curves
***************************************************************************************************************************************************************/


/***************************************************************************************************************************************************************
* Trigonometric Curves
***************************************************************************************************************************************************************/


/***************************************************************************************************************************************************************
* Bezier Curves
***************************************************************************************************************************************************************/


/***************************************************************************************************************************************************************
* Other Parametric Curves
***************************************************************************************************************************************************************/

///** Curve Chain
//***************************************************************************************************************************************************************/
//template<size_t ambient_dim = 2>
//class CurveChain final : public Curve<ambient_dim>
//{
//   using Vector = SVectorR<ambient_dim>;
//
// public:
//   template<class D>
//   CurveChain(const Array<Vector, D>& vertices, bool is_closed = false);
//
//   constexpr void MakeUnitSpeed() noexcept { UnitSpeed_ = true; }
//
// private:
//   constexpr Vector ComputePoint(const Parameter& t) const override;
//
//   constexpr Vector ComputeTangent(const Parameter& t) const override;
//
//   constexpr Vector ComputeBitangent(const Parameter& t) const override;
//
//   constexpr Vector ComputeNormal(const Parameter& t) const override;
//
//   DArray<Curve<ambient_dim>> Curves;
//   DArray<Real>               CumulativeLengths_;
//   Real                       ChainLength_;
//   bool                       Closed_;
//   bool                       UnitSpeed_{false};
//};

}

#include "Curve.tpp"
