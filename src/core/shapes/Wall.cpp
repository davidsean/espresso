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
  along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include "Wall.hpp"
#define SQR(A) ((A) * (A))

namespace Shapes {

int Point::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
//find the distance to the Point 
  *dist=0;
  for (i = 0; i < 3; i++) {
    vec[i] = ppos[i] - m_va[i];
    *dist+= SQR(vec[i]);
  }
  *dist=sqrt(*dist);
  return 0;
}

int Segment::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  double d_a, d_b, d_ba, dot;
  double v_ap[3];
  double v_bp[3];  
  double v_ba[3];  
  double d_line;  
// a segment is an infinite line, bounded by two points.
//first, find the distance to the pts

  // temporay hack to store *ppos
  v_ap[0] = m_pa.va()[0];
  v_ap[1] = m_pa.va()[1];
  v_ap[2] = m_pa.va()[2];
  m_pb.calculate_dist(v_ap, &d_ba, v_ba);

  m_pa.calculate_dist(ppos, &d_a, v_ap);
  m_pb.calculate_dist(ppos, &d_b, v_bp);
 
// we now have three vectors,
// either: 
// 1) ppos is closer to Point va
// 2) ppos is closer to Point vb
// 3) ppos is closer to a point on the lines bettween va and vb
// for find this, we need to check the normalized dot products

double dot_a = 0;
double dot_b = 0;
  for (i = 0; i<3; i++) {
 //   v_ba[i] = m_vb.m_va[i] - m_va.m_va[i];
    dot_a += v_ap[i]*v_ba[i];
    dot_b += v_bp[i]*v_ba[i];
  }

// if both are negative, pt is closer to b
if (dot_a<0 && dot_b<0) {
  //printf("closer to b\n");
  for (i = 0; i<3; i++) {
    vec[i]=v_ap[i];
  }
  *dist=d_a;
  return 0;

// if both are positive, pt is closer to a
} else if ( dot_a>0 && dot_b>0) {
  //printf("closer to a\n");
  for (i = 0; i<3; i++) {
    vec[i]=v_bp[i];
  }
  *dist=d_b;
  return 0;
}
//then, find the distance to line
else {
//printf("closer to line\n");
  d_ba=0;
  dot=0;
 
  for (i = 0; i<3; i++) {
    dot += v_ap[i]*v_ba[i];
    d_ba += SQR(v_ba[i]);
  }
  dot=SQR(dot);
  d_line=sqrt(SQR(d_a)-dot/d_ba);  
 
  *dist=d_b;
  for (i = 0; i < 3; i++)
    vec[i]=v_ap[i]+v_bp[i];
  return 0;
}

}

int Triangle::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  double d_a, d_b, d_c, d_f;
  double v_a[3], v_b[3], v_c[3];
  Vector3d shift;
  Vector3d base_1;
  Vector3d base_2;

// find if it is bound
  // find the maximum point and use it as the origin
  m_sa.m_pa.calculate_dist(ppos, &d_a, v_a);
  m_sb.m_pa.calculate_dist(ppos, &d_b, v_b);
  m_sc.m_pa.calculate_dist(ppos, &d_c, v_c);

  if (d_a>=d_b && d_a>=d_c) {
    shift = m_sa.m_pa;
    base_1 = m_sa;
    base_2 = -m_sc;
  }else if { d_b>=d_a && d_b>=d_c}{
    shift = m_sb.m_pa;
    base_1 = m_sb;
    base_2 = -m_sa;
  }else if {d_c>=d_a && d_c>=d_b}
    shift = m_sc.m_pa;
    base_1 = m_sc;
    base_2 = -m_sb; 
  } else {
    printf("Could not find furthest vertex");      
  }
  if ppos[0]-
  
  
  
  
  //find distances to the segments
  m_sa.calculate_dist(ppos, &d_a,v_a);
  m_sb.calculate_dist(ppos, &d_b,v_b);
  m_sc.calculate_dist(ppos, &d_c,v_c);

  // find the minimum distance to the surface
  d_f = -m_d;
  for (i = 0; i < 3; i++)
    df += ppos[i] * m_n[i];

  
  



  
  return 0;
}


int Wall::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  double d_v[3];
  double d=0;
  printf("hacked wall\n");
  printf("part is %f %f %f \n", ppos[0], ppos[1], ppos[2]);  
  
  // create a Point at 1 0 0  
  Point p1 = Point({5.0, 2.0, 0.0});
  p1.calculate_dist(ppos, &d, d_v);
  //printf("pt1 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
  
  // create a Point at 10 0 0  
  Point p2 = Point({10.0, 1.0, 0.0});
  p2.calculate_dist(ppos, &d, d_v);
  //printf("pt2 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
 
  Segment s1 = Segment(p1,p2);
  s1.calculate_dist(ppos,&d, d_v);
  //printf("s1 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
  return 0;
/*
  *dist = -m_d;
  for (i = 0; i < 3; i++)
    *dist += ppos[i] * m_n[i];

  for (i = 0; i < 3; i++)
    vec[i] = m_n[i] * *dist;
  return 0;
*/
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
