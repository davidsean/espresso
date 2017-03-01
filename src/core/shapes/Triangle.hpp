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

#ifndef SHAPES_TRIANGLE_HPP
#define SHAPES_TRIANGLE_HPP

#include "Shape.hpp"
#include "Vector.hpp"
#include "Segment.hpp"

namespace Shapes {

class Triangle : public Shape {
public:
    
  Triangle() :  m_va({0.0,0.0,0.0}), m_vb({0.0,0.0,0.0}), m_pc({0.0,0.0,0.0})  {}
  
  Triangle(Vector3d a, Vector3d b, Vector3d c) {
    int i;
    // can pre-calculate some dot products here if needed...
    for (i=0; i<3; i++) {
      m_va[i] = a[i] - c[i];
      m_vb[i] = b[i] - c[i];
      m_pc[i] = c[i];
    }
  }

  Triangle(Point a, Point b, Point c) {
    int i;
    // can pre-calculate some dot products here if needed...
    for (i=0; i<3; i++) {
      m_va[i] = a[i] - c[i];
      m_vb[i] = b[i] - c[i];
      m_pc[i] = c[i];
    }
  }
  
  // carful, this changes the direction vector, to set point A, C needs to be defined first
  void set_va(const Vector3d &pos) {
    int i;  
    for (i=0; i<3; i++)
      m_va[i] = pos[i] - m_pc[i];
  }

  // careful, this changes the direction vector, to set point B, C needs to be defined first
  void set_vb(const Vector3d &pos) { 
    int i;  
    for (i=0; i<3; i++)
      m_vb[i] = pos[i] - m_pc[i];
  }


  void set_pc(const Vector3d &pos) { 
    m_pc = pos;
  }
  // a triangle is a corner+two direction vectors,
  // it is stored this way to prevent calculating 
  // the directions every time
  // vector from corner, to pt A
  Vector3d const &va() const { return m_va; }
  // vector from corner, to pt B
  Vector3d const &vb() const { return m_vb; }
  // corner node, pt C
  Vector3d const &pc() const { return m_pc; }


  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;
 
private:

  /** two vectors defining the triangle ( C to A and C to B */
  Vector3d m_va;
  Vector3d m_vb;
  /** position of the corner origin C*/
  Vector3d m_pc;
};


} /* namespace Shapes */

#endif
