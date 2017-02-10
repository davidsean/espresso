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

#include "Wall.hpp"
#define SQR(A) ((A) * (A))

namespace Shapes {
/*
int Point::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
//find the distance to the pts 
  *dist=0;
  for (i = 0; i < 3; i++)
    vec[i] = ppos[i] - m_va[i];
    *dist+= SQR(vec[i]);
  *dist=sqrt(*dist);
  printf("d_a: %f \n", *dist);
  return 0;
}


int Segment::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  double d_a, d_b;
  double v_a[3];
  double v_b[3];  
  
// a segment is an infinite line, bounded by two points.
//first, find the distance to the pts
  m_va.calculate_dist(ppos, &d_a, v_a);
  m_va.calculate_dist(ppos, &d_b, v_b);
  printf("d_a: %f \t d_b: %f \n", d_a, d_b);

//then, find the distance to line
 // TODO
 // TODO
 
//then, compare and choose which shape is the main one.
 // TODO
 if (d_a<d_b) {
  *dist=d_a;
  for (i = 0; i < 3; i++)
    vec[i]=v_a[i];
 }
  return 0;
}
*/




int Wall::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
/*  double bleh[3];
  double blehd=0;
  
//   Point test = Point({m_n[0], m_n[1], m_n[2]});
//   test.calculate_dist(ppos,&blehd,bleh);
  printf("hacked wall");
  printf("partt is %f %f %f", ppos[0], ppos[1], ppos[2]);  
  printf("pt is %f %f %f", m_n[0], m_n[1], m_n[2]);
  */

  *dist = -m_d;
  for (i = 0; i < 3; i++)
    *dist += ppos[i] * m_n[i];

  for (i = 0; i < 3; i++)
    vec[i] = m_n[i] * *dist;
  return 0;
}


} /* namespace Shapes */


// namespace Shapes {
// 
// int Wall::calculate_dist(const double *ppos, double *dist, double *vec) const {
//   int i;
// 
//   *dist = -m_d;
//   for (i = 0; i < 3; i++)
//     *dist += ppos[i] * m_n[i];
// 
//   for (i = 0; i < 3; i++)
//     vec[i] = m_n[i] * *dist;
//   return 0;
// }
// 
// } /* namespace Shapes */
