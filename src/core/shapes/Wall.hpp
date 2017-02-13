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

#ifndef SHAPES_WALL_HPP
#define SHAPES_WALL_HPP

#include "Shape.hpp"
#include "Vector.hpp"

namespace Shapes {

class Point : public Shape {
public:
  Point() : m_va({0., 0., 0.}) {}

  Point(Vector3d a) {
   m_va=a;
   }
 
   // overload the index operator  
   const double& operator[] (const int index) {
    return m_va[index];
  }

  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;                   
  Vector3d const &va() const { return m_va; }
private:
  /** vector defining the point*/
  Vector3d m_va;
};


class Segment : public Shape {
public:
  Segment() : m_pa( Point({0., 0., 0.})), m_pb(Point({0., 0., 0.})) {}

  Segment(Point a, Point b) {
  m_pa=a;
  m_pb=b;
  }

  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;

  Point const &pa() const { return m_pa; }
  Point const &pb() const { return m_pb; }

private:
  /** two points defining the segment*/
  Point m_pa;
  Point m_pb;
};


class Triangle : public Shape {
public:
//  Triangle() : m_pa( Point({0., 0., 0.})), m_pb(Point({0., 0., 0.})) { }
  double n[3];
  Triangle(Point a, Point b, Point c) {
  m_sa = Segment (a,b);
  m_sb = Segment (b,c);
  m_sc = Segment (c,a);
  
  // use the cross product ab and bc to find the normal vector
  m_n = Vector3d( {(a[1]*b[2]-a[2]*b[1]), (a[2]*b[0]-a[0]*b[2]), (a[0]*b[1]-a[1]*b[0])} );
  //beore we normalize, we can find the distance...
  m_d = (-m_n[0]*a[0]+m_n[1]*a[1]+m_n[2]*a[2]);
  m_n.normalize();
  }

  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;

  Segment const &sa() const { return m_sa; }
  Segment const &sb() const { return m_sb; }
  Segment const &sc() const { return m_sc; }

  //Segment const &pa() const { return m_sa.pa(); }
  //Segment const &pb() const { return m_sb.pa(); }
  //Segment const &pc() const { return m_sc.pa(); }


private:
  /** three segments defining the triangle*/
  Segment m_sa;
  Segment m_sb;
  Segment m_sc;

  // normal vector and distance (like infinite planar wall)
  Vector3d m_n;
  double m_d;

};


//} /* namespace Shapes */
//#endif


//namespace Shapes {
class Wall : public Shape {
public:
  Wall() : m_n({1., 0., 0.}), m_d(0.0) {}

  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;

  void set_normal(const Vector3d &normal) {
    m_n = normal;

    m_n.normalize();
  }

  Vector3d const &n() const { return m_n; }

  double const &d() const { return m_d; }

  double &d() { return m_d; }

private:
  /** normal vector on the plane. */
  Vector3d m_n;
  /** distance of the wall from the origin. */
  double m_d;
};

} /* namespace Shapes */

#endif
