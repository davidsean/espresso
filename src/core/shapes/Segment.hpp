/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SHAPES_SEGMENT_HPP
#define SHAPES_SEGMENT_HPP

#include "Shape.hpp"
#include "Point.hpp"

namespace Shapes {

   
class Segment : public Shape {
public:
  Segment() : m_pa( Point({0., 0., 0.})), m_pb(Point({0., 0., 0.})) {}

  Segment(Point a, Point b) {
  m_pa=a;
  m_pb=b;
  }
  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;
  void set_pa(const Vector3d &pos) {
    m_pa = pos;
  }
  
  void set_pb(const Vector3d &pos) {
    m_pb = pos ;
  }
    
  Point const &pa() const { return m_pa; }
  Point const &pb() const { return m_pb; }

private:
  /** two points defining the segment*/
  Point m_pa;
  Point m_pb;
};



} /* namespace Shapes */

#endif
