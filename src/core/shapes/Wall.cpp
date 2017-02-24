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
 
// we now have three options,
// either: 
// 1) ppos is closer to Point va
// 2) ppos is closer to Point vb
// 3) ppos is closer to a point on the lines bettween va and vb
// for find this, we need to check the dot products

double dot_a = 0;
double dot_b = 0;
  for (i = 0; i<3; i++) {
    dot_a += v_ap[i]*v_ba[i];
    dot_b += v_bp[i]*v_ba[i];
  }

// if both are negative, pt is closer to b
if (dot_a<0 && dot_b<0) {
  for (i = 0; i<3; i++) {
    vec[i]=v_bp[i];
  }
  *dist=d_b;
  return 0;

// if both are positive, pt is closer to a
} else if ( dot_a>0 && dot_b>0) {
  for (i = 0; i<3; i++) {
    vec[i]=v_ap[i];
  }
  *dist=d_a;
  return 0;
}
// else, the pt is closer to the line,
// find the distance to line
else {
  double dot=0;
  for (i = 0; i<3; i++) {
    // normalize the line vector
    v_ba[i] /= d_ba;
    //find the projection
    dot += v_ap[i]*(-1.*v_ba[i]);
  }
  *dist=0;
  for (i = 0; i<3; i++) {
    // remove the projection and we have distance vector
    vec[i] = v_ap[i] + (dot*v_ba[i])  ;
    *dist+=SQR(vec[i]);
  }
  *dist=sqrt(*dist);
  return 0;
}

}




int Triangle2::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  
  Vector3d t_0=m_va;  
  Vector3d t_1=m_vb;
  Vector3d t_2;
  for (i = 0; i < 3; i++) {
    t_2[i] = ppos[i] - m_pc[i];
  }
  dot00 = t_0.dot(t_0);
  dot01 = t_0.dot(t_1);
  dot02 = t_0.dot(t_2);
  dot11 = t_1.dot(t_1); 
  dot12 = t_1.dot(t_2);

  double norm=1./(dot00*dot11 - SQR(dot01));
  u = (dot11*dot02 - dot01*dot12)*norm;
  v = (dot00*dot12 - dot01*dot02)*norm;
  if ((u>=0) && (v>=0) && (u+v<1)) {
    // The triangualr surface is the target
    // find the norm 
    Vector3d n =t_0.cross(t_1);
    //beore we normalize, we can find the distance...
    *dist = (n.dot(m_pc));
    n.normalize();
    for (i = 0; i < 3; i++)
      *dist += ppos[i] * n[i];
    for (i = 0; i < 3; i++)
      vec[i] = n[i] * (*dist);
    return 0;
  } elseif ((u>=0) && (v>=0)) {
  // the edge  AB is the target
    Point a,b;
    for (i = 0; i < 3; i++)
      a[i] = m_pc[i] +m_va[i];
      b[i] = m_pc[i] +m_vb[i];
    s_ab = Segment(a,b);
    s_ab.calculate_dist(ppos,dist,vec);
    return 0;
  } elseif ((u<=0) && (v>=0)) {
  // the edge  CB is the target
    Point b;
    for (i = 0; i < 3; i++)
      b[i] = m_pc[i] +m_vb[i];
    s_cb = Segment(m_pc,b);
    s_cb.calculate_dist(ppos,dist,vec);
    return 0;
  } elseif ((v<=0) ) {
  // the edge  AC is the target
    Point a;
    for (i = 0; i < 3; i++)
      a[i] = m_pc[i] +m_vb[i];
    s_ac = Segment(m_pc,a);
    s_ac.calculate_dist(ppos,dist,vec);
    return 0;
  } else {
     printf("error, did not find correct target");
     return 1;
  }
}


int Triangle::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  // testing this method: http://blackpawn.com/texts/pointinpoly/
  double t_0[3];
  double t_1[3];
  double t_2[3];
  double dot00=0;
  double dot01=0;
  double dot02=0;
  double dot11=0;
  double dot12=0;
  double u,v;
  double a[3];

// get some temp vectors
  for (i = 0; i < 3; i++) {
    t_0[i] = m_pc.va()[i] - m_pa.va()[i];
    t_1[i] = m_pb.va()[i] - m_pa.va()[i];
    t_2[i] = ppos[i] - m_pa.va()[i];
    dot00 += SQR(t_0[i]);
    dot01 += t_0[i]*t_1[i];
    dot02 += t_0[i]*t_2[i];
    dot11 += SQR(t_1[i]);
    dot12 += t_1[i]*t_2[i];
  }
  double norm=1./(dot00*dot11 - SQR(dot01));
  u = (dot11*dot02 - dot01*dot12)*norm;
  v = (dot00*dot12 - dot01*dot02)*norm;
  if ((u>=0) && (v>=0) && (u+v<1)) {
    // find the minimum distance to the surface
    *dist = -m_d;
    for (i = 0; i < 3; i++)
      *dist += ppos[i] * m_n[i];
    vec[0]=m_n[0];
    vec[1]=m_n[1];
    vec[2]=m_n[2];
    return 0;
  }
 
// I can do better, but let's worry about performance at some later time...
   else {
  // if outside the planar triangle, find the closest segment
  double d_sa, d_sb, d_sc;
  double v_sa[3], v_sb[3], v_sc[3];
  m_sa.calculate_dist(ppos, &d_sa,v_sa);
  m_sb.calculate_dist(ppos, &d_sb,v_sb);
  m_sc.calculate_dist(ppos, &d_sc,v_sc);
  
    if (d_sa<=d_sb && d_sa<=d_sc) {
      *dist=d_sa;
      vec=v_sa;
      return 0;
    } else if (d_sb<=d_sa && d_sb<=d_sc ) {
      *dist=d_sb;
      vec=v_sb;
      *dist=d_sb;
      return 0;
    } else if (d_sc<=d_sa && d_sc<=d_sb) {
      *dist=d_sc;
      vec=v_sc;
     *dist=d_sc;
      return 0;
    } 
  }
  return 1;
}


int Wall::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  double d_v[3];
  double d=0;

  // create a Point at 1 0 0  
  Point p1 = Point({2.0000, 1.000, 0.0});
  //p1.calculate_dist(ppos, &d, d_v);
  //printf("pt1 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
  
  // create a Point at 10 0 0  
  Point p2 = Point({8.00, 2.000, 0.0});
  //p2.calculate_dist(ppos, &d, d_v);
  //printf("pt2 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
 
  //Segment s1 = Segment(p1,p2);
  //s1.calculate_dist(ppos, &d, d_v);
  //printf("s1 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
 
  Point p3 = Point({3.0, 8.0, 0.0});
  //p3.calculate_dist(ppos, &d, d_v);

    
  Triangle2 t1 = Triangle2(p1,p2,p3);
  t1.calculate_dist(ppos, &d, d_v);
  
  for (i = 0; i < 3; i++)
    vec[i]=d_v[i];
  *dist=d;

/*
   *dist = -m_d;
   for (i = 0; i < 3; i++)
     *dist += ppos[i] * m_n[i];
 
   for (i = 0; i < 3; i++)
     vec[i] = m_n[i] * *dist;
     
*/
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
