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
    vec[i]=v_ap[i];
  }
  *dist=d_b;
  return 0;

// if both are positive, pt is closer to a
} else if ( dot_a>0 && dot_b>0) {
  for (i = 0; i<3; i++) {
    vec[i]=v_bp[i];
  }
  *dist=d_a;
  return 0;
}
// else, the pt is closer to the line,
// find the distance to line
else {
  d_ba=0;
  dot=0;
  for (i = 0; i<3; i++) {
    dot += v_ap[i]*v_ba[i];
    d_ba += SQR(v_ba[i]);
  }
  dot=SQR(dot);
  d_line=sqrt(SQR(d_a)-dot/d_ba);   
  *dist=d_line;
  
  // we have the distance, but we will need the vector
  // there is probably a faster way, but this works.
  double norm_a=0;
  double norm_b=0;

  for (i = 0; i < 3; i++) {
    norm_a+=SQR(v_ap[i]);
    norm_b+=SQR(v_bp[i]);
  }
  norm_a=1./sqrt(norm_a); 
  norm_b=1./sqrt(norm_b); 
  for (i = 0; i < 3; i++) {
    v_ap[i]*=norm_a;
    v_bp[i]*=norm_b;
  }
  
  // add up the two vectors to find twice the distance
  for (i = 0; i < 3; i++)
    vec[i]=0.5*(v_ap[i]+v_bp[i])*(d_line);
  printf("| %f, %f, %f| = %f", vec[0],vec[1],vec[2],*dist);
  
  // simple method, make sure results are the same
  m_pa.calculate_dist(ppos, &d_a, v_ap);
  m_pb.calculate_dist(ppos, &d_b, v_bp);
  *dist=0;
  // add up the two vectors to find twice the distance
  for (i = 0; i < 3; i++)
    vec[i]=0.5*(v_ap[i]+v_bp[i]);
    *dist+=SQR(vec[i]);
  *dist=sqrt(*dist);  
  printf("| %f, %f, %f| = %f", vec[0],vec[1],vec[2],*dist);

    
    
  return 0;
}

}


int Triangle::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  double d_a, d_b, d_c, d_f;
  double v_a[3], v_b[3], v_c[3];

 // if ppos[0]-
  //find distances to the segments
  m_sa.calculate_dist(ppos, &d_a,v_a);
  m_sb.calculate_dist(ppos, &d_b,v_b);
  m_sc.calculate_dist(ppos, &d_c,v_c);

  // testing this method: http://blackpawn.com/texts/pointinpoly/
  double t_0[3];
  double t_1[3];
  double t_2[3];
  double sc00=0;
  double sc01=0;
  double sc02=0;
  double sc11=0;
  double sc12=0;
  double u,v;
  double a[3];


// get some temp vectors
  for (i = 0; i < 3; i++) {
  //  t_0[i] = m_sc.pa()[i] - m_ma.pa()[i];
  //  t_1[i] = m_sb.pa()[i] - m_sa.pa()[i];
    t_0[i] = 1;
    t_1[i] = 1;
    sc00 += SQR(t_0[i]);
    sc01 += t_0[i]*t_1[i];
    sc02 += t_0[i]*v_a[i];
    sc11 += SQR(t_1[i]);
    sc12 += t_1[i]*v_a[i];
  }
  double norm=1./(sc00*sc11 - SQR(sc01));
  u = (sc11*sc02 - sc01*sc12)*norm;
  v = (sc00*sc12 - sc01*sc02)*norm;
  //if ( (u>=0) && (v>=0) && (u+v<1)) {
  if (0) {
    // find the minimum distance to the surface
    d_f = -m_d;
    for (i = 0; i < 3; i++)
    d_f += ppos[i] * m_n[i];
    *dist=d_f;
    vec[0]=m_n[0];
    vec[1]=m_n[1];
    vec[2]=m_n[2];
    //hack!
    *dist=0;
    vec[0]=0;
    vec[1]=0;
    vec[2]=0;
    return 0;
  } else if (d_a<=d_b && d_a<=d_c) {
    *dist=d_a;
    vec=v_a;
    return 0;
  } else if (d_b<=d_a && d_b<=d_c) {
    *dist=d_a;
    vec=v_a;
    return 0;
  } else if (d_c<=d_a && d_c<=d_b) {
    *dist=d_c;
    vec=v_c;
    return 0;
  } else {
    printf("Not sure!");
    //hack!
    *dist=0;
    vec[0]=0;
    vec[1]=0;
    vec[2]=0;
    return 0;
  }
 
// I can do better, but let's worry about performance at some later time...
 
  return 1;
}


int Wall::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  double d_v[3];
  double d=0;
  //printf("hacked wall\n");
  //printf("part is %f %f %f \n", ppos[0], ppos[1], ppos[2]);  
  
  // create a Point at 1 0 0  
  Point p1 = Point({1.000010, 0.00001, 0.0});
  //p1.calculate_dist(ppos, &d, d_v);
  //printf("pt1 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
  
  // create a Point at 10 0 0  
  Point p2 = Point({5.00001, 3.00001, 0.0});
  //p2.calculate_dist(ppos, &d, d_v);
  //printf("pt2 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
 
  Point p3 = Point({1.00001, 5.00001, 0.0});
  //Segment s1 = Segment(p1,p2);
  //s1.calculate_dist(ppos,&d, d_v);
  //printf("s1 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
 
  //s1.calculate_dist(ppos,dist, vec);
  Triangle t1 = Triangle(p1,p2,p3);
  t1.calculate_dist(ppos,dist, vec);
  //printf("s1 is %f %f %f \t d:%f\n", vec[0], vec[1], vec[2], *dist);
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
