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



int Square::calculate_dist(const double *ppos, double *dist, double *vec) const {
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
  if ((u>=0) && (v>=0) && (u<=1) && (v<=1)) {
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
  } else if ((u>=1) && (v<=1) && (v>=0)) {
  // the edge  AD is the target
    Vector3d a,d;
    for (i = 0; i < 3; i++) {
      d[i] = m_pc[i] + m_va[i] + m_vb[i];
      a[i] = m_pc[i] + m_va[i];
    }
    Segment s_ad = Segment(a,d);
    s_ad.calculate_dist(ppos,dist,vec);
    //*dist=0.955;
    return 0; 
} else if ((u>=0) && (v>=0) ) {
  // the edge  DB is the target
    Vector3d d,b;
    for (i = 0; i < 3; i++) {
      d[i] = m_pc[i] + m_va[i] + m_vb[i];
      b[i] = m_pc[i] + m_vb[i];
    }
    Segment s_db = Segment(d,b);
    s_db.calculate_dist(ppos,dist,vec);
    //*dist=0.970;
    return 0;
} else if ((u<=0) && (v>=0)) {
  // the edge  CB is the target
    Vector3d d;
    for (i = 0; i < 3; i++) {
      d[i] = m_pc[i] + m_vb[i];
    }
    Segment s_cb = Segment(m_pc,d);
    s_cb.calculate_dist(ppos,dist,vec);
    //*dist=0.955;
    return 0; 
} else if ((v<=0)) {
  // the edge  CA is the target
    Vector3d a;
    for (i = 0; i < 3; i++) {
      a[i] = m_pc[i] + m_va[i];
    }
    Segment s_ca = Segment(m_pc,a);
    s_ca.calculate_dist(ppos,dist,vec);
    //*dist=0.955;
    return 0; 
 } else {
     printf("error, did not find correct target");
     return 1;
  }
}



int Voxel::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int j,i,i_min;
  double dx,dy,dz;
  double dx2,dy2,dz2;
  
  double delta[3];
  double delta2[3];
 
  for (i=0; i<3; i++) {
    delta[i] = ppos[i]-m_va[i];
    delta2[i] = SQR(delta[i]);
  }
  //if any distances are negative, this eliminates 3 of the 6 faces
  for (j=0; j<3; j++) {
    if (delta[j]<=0){
      //printf("first corner!\n");
      //i_min determines which square of the 3 is active
      i_min=0;
      for (i=0; i<3; i++) {
        if (delta2[i]<delta2[i_min])
          i_min=i;
      }
      //create the active Square
      Vector3d a,b;
      for (i=0; i<3; i++) {
          a[i]=m_va[i];
          b[i]=m_va[i];
      }
      // dont change the active component!
      a[(i_min+1)%3]+=m_l;
      b[(i_min+2)%3]+=m_l;
      Square voxel_face = Square(a,b,m_va);
      voxel_face.calculate_dist(ppos,dist,vec);
      return 0;
    }
  }
  
  // two options remain: 1)outside the 3 other faces or 2) Inside the voxel.
  // test outside first, as is it most probable
  for (i=0; i<3; i++) {
    delta[i]-=m_l;
    delta2[i] = SQR(delta[i]);
  }
  // if any are positive, it is outside.
  for (j=0; j<3; j++) {
    if (delta[j]>=0) {
      //printf("second corner!\n");
      //i_min determines which square of the 3 is active
      i_min=0;
      for (i=0; i<3; i++) {
        if (delta2[i]<delta2[i_min])
          i_min=i;
      }
      //create the active Square
      Vector3d a,b,c;
      for (i=0; i<3; i++) {
          a[i]=m_va[i]+m_l;
          b[i]=a[i]; // don't need to add a second time!
          c[i]=a[i]; // ditto
      }
      a[(i_min+1)%3]-=m_l;
      b[(i_min+2)%3]-=m_l;
      Square voxel_face = Square(a,b,c);
      voxel_face.calculate_dist(ppos,dist,vec);
      return 0;
    }
  }
  // if the control flow reached here, the particle is INSIDE the voxel.
  // not sure what do do now except return a negative distance
  // this will at least set the LB boundaries correctly.
  // TODO do something better
  //printf("Inside\n");
  *dist=0.955;
  return 0;
}



int Wall::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i;
  double d_v[3];
  double d=0;

  // create a Point at 1 0 0  
  //Point p1 = Point({2.00, 2.0, 0.0});
  //p1.calculate_dist(ppos, &d, d_v);
  //printf("pt1 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
  
  // create a Point at 10 0 0  
  //Point p2 = Point({8.0, 8.0, 0.0});
  //p2.calculate_dist(ppos, &d, d_v);
  //printf("pt2 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
 
  //Segment s1 = Segment(p1,p2);
  //s1.calculate_dist(ppos, &d, d_v);
  //printf("s1 is %f %f %f \t d:%f\n", d_v[0], d_v[1], d_v[2], d);
 
  //Point p3 = Point({2.0, 8.1, 0.0});
  //p3.calculate_dist(ppos, &d, d_v);
    
  //Square S1 = Square(p1,p2,p3);
  //S1.calculate_dist(ppos, dist, vec);
  
  //Voxel v1= Voxel({4.0,4.0,0.01},3.0);
  //v1.calculate_dist(ppos, dist, vec);
  // Hack the wall framework to add a voxel 2x2 at vector: d_v*d
  Voxel v1= Voxel({m_n[0]*m_d, m_n[1]*m_d, m_n[2]*m_d},2.0);
  v1.calculate_dist(ppos,dist,vec);
  return 0;


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
