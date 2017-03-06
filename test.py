#
# Copyright (C) 2013,2014 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from __future__ import print_function
import espressomd
from espressomd import thermostat
from espressomd import shapes
from espressomd import code_info
from espressomd import interactions
from espressomd import constraints

from espressomd.visualizationOpenGL import *
from threading import Thread

import numpy as np

print("""
=======================================================
=                      test.py                         =
=======================================================

Program Information:""")
print(code_info.features())
system = espressomd.System()


visualizer = openGLLive(system,
        constraint_type_colors = [[0.0,1.0,0.0,0.5]],
        background_color = [1.0,1.0,1.0], 
        )

system.box_l =[50,50,50]
system.time_step = 0.01
system.cell_system.skin = 10
system.thermostat.set_langevin(kT=0.2, gamma=1.0)

system.non_bonded_inter[1, 0].lennard_jones.set_params(
     epsilon=1, sigma=1,
     cutoff=2**(1. / 6), shift="auto")




harmonic = interactions.HarmonicBond(k=10., r_0=1.0)
system.bonded_inter.add(harmonic)

#fene = interactions.FeneBond(k=30., d_r_max=1.5)
#system.bonded_inter.add(fene)

#system.cell_system.set_n_square(use_verlet_lists=False)


#system.part.add(id=0, pos=[box_l/2.0,box_l/2.0,box_l/2.0], fix=[1,1,1])
#system.part.add(id=0, pos=[2,1,1],type=0, ext_force=[0,0,0])


# helper
def draw_triangle(fp, a,  b,  c) :
  fp.write("draw triangle \"{} {} {}\" \"{} {} {}\" \"{} {} {}\"\n".format(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]))

def draw_xsquare(fp, a, b) :
   draw_triangle(fp,a, [a[0], b[1], a[2]], b)
   draw_triangle(fp,a, [a[0], a[1], b[2]], b)

def draw_ysquare(fp, a, b) :
   draw_triangle(fp, a, [a[0], a[1], b[2]], b)
   draw_triangle(fp, a, [b[0], a[1], a[2]], b)

def draw_zsquare(fp, a, b) :
   draw_triangle(fp, a, [a[0], b[1], a[2]], b)
   draw_triangle(fp, a, [b[0], a[1], a[2]], b)


def draw_Voxel(fp, pos, l) :
  a=l*.5
  for corner1 in (np.zeros(3), np.ones(3)):
    corner2=np.ones(3)-corner1
    for index in range(3):
      corner2[index] = 1-corner2[index]
      if index==0:
        draw_xsquare(fp, pos+corner1*a, pos+corner2*a)
      if index==1:
        draw_ysquare(fp, pos+corner1*a, pos+corner2*a)
      if index==2:
        draw_zsquare(fp, pos+corner1*a, pos+corner2*a)
      corner2[index] = 1-corner2[index]


def open_ASCII_STL(fp):
  vertex=[]
  for line in fp:
    parsed=line.split()
    if len(parsed)>1:
      if (parsed[0]=="vertex"):
        vertex.append(float(parsed[1]))
        vertex.append(float(parsed[2]))
        vertex.append(float(parsed[3]))
  return np.array(vertex)




zero=0.01
fp=open("test.con" , "w")
fp.write("draw delete all\n")
fp.write("draw color 3\n")
fp.write("draw material Opaque\n")


#voxels=np.loadtxt("cow.voxel")
#print("Found {} voxels".format(np.shape(voxels)[0]))
#for v in voxels:
#  draw_Voxel(fp,v,1.)


#vert= open_ASCII_STL(open("cow_ASCII.stl","r"))
vert= open_ASCII_STL(open("icosahedron_ASCII.stl","r"))
#vert= open_ASCII_STL(open("ESPResSo_ASCII.stl","r"))
vert = np.reshape(vert,(len(vert)/9,3,3))

print (" Need to print {} triangles".format(len(vert[:,0,0])))
# rescale?
xmin=np.min(vert[:,:,0])
ymin=np.min(vert[:,:,1])
zmin=np.min(vert[:,:,2])

vert[:,:,0]-=xmin
vert[:,:,1]-=ymin
vert[:,:,2]-=zmin


xmax=np.max(vert[:,:,0])
ymax=np.max(vert[:,:,1])
zmax=np.max(vert[:,:,2])

max= np.max([xmax, ymax, zmax])

rescale=(system.box_l[0])*1./(max+2)
vert*=rescale*.75

# center in y
ymean=np.mean(vert[:,:,1])
vert[:,:,1]-=ymean
vert[:,:,1]+=system.box_l[1]*0.5

# center in x
xmean=np.mean(vert[:,:,0])
vert[:,:,0]-=xmean
vert[:,:,0]+=system.box_l[0]*0.5


# shift x by 1
vert[:,:,0]+=1.0

# shift z by 1
vert[:,:,2]+=1.0


#p1=shapes.Point(pos=[8,1.1,1.1])
#p2=shapes.Point(pos=[1.1,1.1,1.1])
#p3=shapes.Point(pos=[5,5,1.1])
#system.constraints.add(particle_type=1, penetrable=1, shape=p1)
#system.constraints.add(particle_type=1, penetrable=1, shape=p2)
#system.constraints.add(particle_type=1, penetrable=1, shape=p3)

#s1=shapes.Segment(a=[5,5,1.1], b=[9,1.1,1.1])
#system.constraints.add(particle_type=1, penetrable=1, shape=s1)

#tri=shapes.Triangle(a=[8,1.1,1.1], b=[1.1,1.1,1.1], c=[5,5,1.1])
#tri=shapes.Triangle(a=[1,1.1,1.1], b=[8.1,1.1,1.1], c=[5,5,1.1])



for v in vert:
  a=v[0]
  b=v[1]
  c=v[2]
#  tri=shapes.Triangle(a=[v[0,0], v[0,1], v[0,2]], b=[v[1,0], v[1,1], v[1,2]], c=[v[2,0], v[2,1], v[2,2]])
  tri=shapes.Triangle(a=[float(a[0]), float(a[1]), float(a[2])], b=[float(b[0]), float(b[1]), float(b[2])], c=[float(c[0]), float(c[1]), float(c[2])])
  system.constraints.add(particle_type=1, penetrable=1, shape=tri)



for c in system.constraints.get():
    if (c.get_params()['shape'].name()=='Shapes::Triangle' ):
      #print("Found a triangle contraint, drawing...")
      a=c.get_params()['shape'].get_params()['a']
      b=c.get_params()['shape'].get_params()['b']
      c=c.get_params()['shape'].get_params()['c']
      draw_triangle(fp, a, b, c)

    else:
      print("Found unknown contraint shape: ")
      print(c.get_params()['shape'].name())



fp.close()

p=0

z=system.box_l[2]-18 #use 8 for ESPResSo
#z=0.5*system.box_l[2] # for cow

for x in np.arange(0,system.box_l[0],0.97):
  for y in np.arange(0,system.box_l[1],0.97):
    system.part.add(id=p, pos=[x,y,z], type=0, ext_force=[0,0,-.1])
    p+=1

num_mono=p-1
p=0
for x in np.arange(0,system.box_l[0],.97):
  for y in np.arange(0,system.box_l[1],.97):
    if (p<num_mono):  
      system.part[p].add_bond((harmonic,p+1))
      p+=1
  if (p<num_mono):  
    system.part[p-1].delete_all_bonds()


#visualizer.update()
#visualizer.start()


sim_ID="test"
vtf_file=open("{}.vtf".format(sim_ID),"w")
vtf_file.write("unitcell {} {} {}\n".format(system.box_l[0], system.box_l[1], system.box_l[2]))
for p in system.part:
  if p.type==0: name="mono"
  else: name="unknown"
  vtf_file.write("atom {} radius 1 name {} type {} \n".format(p.id, name, p.type))
for p in system.part:
  for b in p.bonds:
    vtf_file.write("bond {}:{}\n".format(p.id, b[1]))

vtf_file.write("\ntimestep ordered\n")
for p in system.part:
  vtf_file.write("{} {} {} \n".format(p.pos[0], p.pos[1], p.pos[2]))


print("Simulating...")
for t in range(2000):
  print (t)
  system.integrator.run(10)
  vtf_file.write("\ntimestep ordered\n")
  for p in system.part:
    vtf_file.write("{} {} {} \n".format(p.pos[0], p.pos[1], p.pos[2]))









exit()

