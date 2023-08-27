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
#include "LinearAlgebra/include/VectorOperations.h"

namespace aprn::mnfld {

/***************************************************************************************************************************************************************
* Curve Class Implementation
***************************************************************************************************************************************************************/

template<size_t D>
constexpr SVectorR<D>
Curve<D>::Binormal(const Vector& tangent, const Vector& normal) const { return CrossProduct(tangent, normal); }

/***************************************************************************************************************************************************************
* Linear/Piecewise Linear Curves
***************************************************************************************************************************************************************/

/** Line
***************************************************************************************************************************************************************/
template<size_t D>
constexpr Line<D>::Line(const Vector& direction, const Vector& point)
   : Direction(direction), Start(point), DirectionNorm_(Magnitude(Direction)), Normaliser_(One / DirectionNorm_) {}

template<size_t D>
constexpr SVectorR<D>
Line<D>::Point(const Real t) const { return Start + t * (this->UnitSpeed_ ? Normaliser_ : One) * Direction; }

template<size_t D>
constexpr SVectorR<D>
Line<D>::Tangent([[maybe_unused]] const Real t) const { return (this->UnitSpeed_ ? Normaliser_ : One) * Direction; }

template<size_t D>
constexpr SVectorR<D>
Line<D>::Normal(const Real t) const
{
   EXIT("TODO")
   return Direction;
}

/** Ray
***************************************************************************************************************************************************************/
template<size_t D>
constexpr Ray<D>::Ray(const Vector& direction, const Vector& start)
   : Line<D>::Line(direction, start) {}

template<size_t D>
constexpr SVectorR<D>
Ray<D>::Point(const Real t) const
{
   return Positive(t) ? Line<D>::Point(t) : throw std::domain_error("The parameter must be positive for rays.");
}

/** Segment
***************************************************************************************************************************************************************/
template<size_t D>
constexpr LineSegment<D>::LineSegment(const Vector& start, const Vector& end)
   : Line<D>::Line(end - start, start) {}

template<size_t D>
constexpr SVectorR<D>
LineSegment<D>::Point(const Real t) const
{
   const Real max_bound = this->UnitSpeed_ ? Length() : One;
   return isBounded<true, true>(t, Zero, max_bound) ? Line<D>::Point(t) :
          throw std::domain_error("The parameter must be in the range [0, " + ToString(max_bound) + "] for this segment.");
}

/** SegmentChain
***************************************************************************************************************************************************************/
template<size_t Dim>
template<class D>
LineSegmentChain<Dim>::LineSegmentChain(const Array<Vector, D>& vertex_list, const bool is_closed)
   : Closed_(is_closed)
{
   const auto& vertices = vertex_list.Derived();

   Segments_.reserve(vertices.size());
   CumulativeLengths_.reserve(vertices.size());
   ChainLength_ = Zero;

   FOR(i, vertices.size() - static_cast<size_t>(!Closed_))
   {
      Segments_.emplace_back(vertices[i], vertices[(i + 1) % vertices.size()]);
      Segments_.back().MakeUnitSpeed();

      ChainLength_ += Segments_.back().Length();
      CumulativeLengths_.emplace_back(ChainLength_);
   }
}

template<size_t D>
constexpr SVectorR<D>
LineSegmentChain<D>::Point(const Real t) const
{
   const Real upper_bound  = this->UnitSpeed_ ? ChainLength_ : One;
   const Real param_length = isBounded<true, true, true>(t, Zero, upper_bound) ? t * (this->UnitSpeed_ ? One : ChainLength_) :
                             throw std::domain_error("The parameter must be in the range [0, " + ToString(upper_bound) + "] for this segment.");
   const auto iter  = std::find_if(CumulativeLengths_.begin(), CumulativeLengths_.end(), [param_length](auto l){ return param_length <= l; });
   const auto index = std::distance(CumulativeLengths_.begin(), iter);
   const Real param = param_length - (index != 0 ? CumulativeLengths_[index - 1] : Zero);

   return isBounded<true, true>(param, Zero, Segments_[index].Length()) ? Segments_[index].Point(param) :
          throw std::domain_error("The parameter for segment " + ToString(index) + " in the chain is out of bounds.");
}

template<size_t D>
constexpr SVectorR<D>
LineSegmentChain<D>::Tangent(const Real t) const
{
   EXIT("TODO")
}

template<size_t D>
constexpr SVectorR<D>
LineSegmentChain<D>::Normal(const Real t) const
{
   EXIT("TODO")
}

/***************************************************************************************************************************************************************
* Circular/Elliptical Curves
***************************************************************************************************************************************************************/

/** Circle
***************************************************************************************************************************************************************/
template<size_t D>
Circle<D>::Circle(const Real radius, const Vector& centre)
   : Circle(radius, Zero, centre) {}

template<size_t D>
Circle<D>::Circle(const Real radius, const Real start_angle, const Vector& centre)
   : Centre_(centre), Radius_(radius), StartAngle_(start_angle), Normaliser_(One / Radius_) { ASSERT(Positive(radius), "A circle's radius cannot be negative.") }

template<size_t D>
constexpr SVectorR<D>
Circle<D>::Point(const Real t) const
{
   const Real max_bound = this->UnitSpeed_ ? TwoPi * Radius_ : One;
   ASSERT((isBounded<true, true>(t, Zero, max_bound)), "The parameter exceeds the expected bounds.")

   const auto theta = Angle(t);
   return ToVector<D>(SVectorR3{Radius_ * std::cos(theta), Radius_ * std::sin(theta), Zero}) + Centre_;
}

template<size_t D>
constexpr SVectorR<D>
Circle<D>::Tangent(const Real t) const
{
   const auto theta = Angle(t);
   return ToVector<D>(SVectorR3{-Radius_ * std::sin(theta), Radius_ * std::cos(theta), Zero});
}

template<size_t D>
constexpr SVectorR<D>
Circle<D>::Normal(const Real t) const
{
   const auto theta = Angle(t);
   return ToVector<D>(SVectorR3{-Radius_ * std::cos(theta), -Radius_ * std::sin(theta), Zero});
}

template<size_t D>
constexpr Real
Circle<D>::Angle(const Real t) const { return StartAngle_ + t * (this->UnitSpeed_ ? Normaliser_ : TwoPi); }

/** Circular Arc
***************************************************************************************************************************************************************/
template<size_t D>
Arc<D>::Arc(const Real radius, const Real angle, const Vector& centre)
   : Arc(radius, Zero, angle, centre) {}

template<size_t D>
Arc<D>::Arc(const Real radius, const Real start_angle, const Real end_angle, const Vector& centre)
   : Circle<D>(radius, start_angle, centre), EndAngle_(end_angle)
{
   ASSERT(Positive(radius), "An arc's radius cannot be negative.")
   ASSERT((isBounded<true, true>(start_angle, Zero, TwoPi)), "An arc's start angle must be in the range [0, 2*PI].")
   ASSERT((isBounded<true, true>(end_angle, Zero, TwoPi)), "An arc's end angle must be in the range [0, 2*PI].")
}

template<size_t D>
constexpr SVectorR<D>
Arc<D>::Point(const Real t) const
{
   CheckAngle(t);
   return Circle<D>::Point(t);
}

template<size_t D>
constexpr SVectorR<D>
Arc<D>::Tangent(const Real t) const
{
   CheckAngle(t);
   return Circle<D>::Tangent(t);
}

template<size_t D>
constexpr SVectorR<D>
Arc<D>::Normal(const Real t) const
{
   CheckAngle(t);
   return Circle<D>::Normal(t);
}

template<size_t D>
constexpr void
Arc<D>::CheckAngle(const Real t) const
{
   const auto theta = this->Angle(t);
   const auto min_max = std::minmax(this->StartAngle_, EndAngle_);

   ASSERT((isBounded<true, true>(theta, min_max.first, min_max.second)), "The computed angle lies outside the min/max angle range of the arc.", theta, " ", min_max.first, " ", min_max.second)
}

/** Ellipse
***************************************************************************************************************************************************************/
template<size_t D>
Ellipse<D>::Ellipse(const Real radius_x, const Real radius_y, const Vector& centre)
   : Centre_(centre), RadiusX_(radius_x), RadiusY_(radius_y) { ASSERT(Positive(radius_x) && Positive(radius_y), "An ellipse's radii cannot be negative.") }

template<size_t D>
constexpr SVectorR<D>
Ellipse<D>::Point(const Real t) const
{
   return ToVector<D>(SVectorR3{RadiusX_ * std::cos(t), RadiusY_ * std::sin(t), Zero}) + Centre_;
}

template<size_t D>
constexpr SVectorR<D>
Ellipse<D>::Tangent(const Real t) const
{
   EXIT("TODO")
}

template<size_t D>
constexpr SVectorR<D>
Ellipse<D>::Normal(const Real t) const
{
   EXIT("TODO")
}

}