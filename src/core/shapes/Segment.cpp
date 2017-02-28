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

#include "Segment.hpp"
#define SQR(A) ((A) * (A))

namespace Shapes {

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



} /* namespace Shapes */
