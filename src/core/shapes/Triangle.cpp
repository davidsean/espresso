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

#include "Triangle.hpp"
#define SQR(A) ((A) * (A))

namespace Shapes {
int Triangle::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  
  Vector3d t_0=m_va;  
  Vector3d t_1=m_vb;
  Vector3d t_2;
  for (i = 0; i < 3; i++) {
    t_2[i] = ppos[i] - m_pc[i];
  }
  // I may want to pre-calculate some of these to keep in the Triangle class
  double dot00 = t_0.dot(t_0); // can precalculate this one
  double dot01 = t_0.dot(t_1); // can precalculate this one
  double dot02 = t_0.dot(t_2);
  double dot11 = t_1.dot(t_1); // can precalculate this one
  double dot12 = t_1.dot(t_2);

  double norm=1./(dot00*dot11 - SQR(dot01)); // can precalculate this one
  double u = (dot11*dot02 - dot01*dot12)*norm;
  double v = (dot00*dot12 - dot01*dot02)*norm;
  if ((u>=0) && (v>=0) && (u+v<1)) {
    // The triangualr surface is the target
    // find the norm 
    Vector3d n;
    n.cross(t_1,t_0,n); // can precalculate this one
    //beore we normalize, we can find the distance...
    *dist = (n.dot(m_pc)); 
    n.normalize();
    for (i = 0; i < 3; i++)
      *dist += ppos[i] * n[i];
    for (i = 0; i < 3; i++)
      vec[i] = n[i] * (*dist);
    return 0;
  } else if ((u>=0) && (v>=0)) {
  // the edge  AB is the target
    Vector3d a;
    Vector3d b;
    for (i = 0; i < 3; i++) {
      a[i] = m_pc[i] + m_va[i];
      b[i] = m_pc[i] + m_vb[i];
    }
    Segment s_ab = Segment(a,b);
    s_ab.calculate_dist(ppos,dist,vec);
    return 0;
  } else if ((u<=0) && (v>=0)) {
  // the edge  CB is the target
    Vector3d b;
    for (i = 0; i < 3; i++)
      b[i] = m_pc[i] +m_vb[i];
    Segment s_cb = Segment(m_pc,b);
    s_cb.calculate_dist(ppos,dist,vec);
    return 0;
  } else if ((v<=0) ) {
  // the edge  AC is the target
    Vector3d a;
    for (i = 0; i < 3; i++)
      a[i] = m_pc[i] +m_va[i];
    Segment s_ac = Segment(a, m_pc);
    s_ac.calculate_dist(ppos,dist,vec);
    return 0;
  } else {
     printf("error, did not find correct target");
     return 1;
  }
}


} /* namespace Shapes */
