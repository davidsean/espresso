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

system.box_l =[10,10,10]
system.time_step = 0.01
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
    
#  vertex 432.608 69.762 -218.242
#  vertex 432.124 69.754 -217.288
#  vertex 432.616 70.246 -217.288


zero=0.01
fp=open("test.con" , "w")
fp.write("draw delete all\n")
fp.write("draw color 3\n")
fp.write("draw material Opaque\n")


vert= open_ASCII_STL(open("cow_ASCII.stl","r"))
vert= open_ASCII_STL(open("icosahedron_ASCII.stl","r"))
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

# shift x by 1
vert[:,:,0]+=1.0

# shift z by 1
vert[:,:,2]+=1.0


#for v in vert:
#  draw_triangle(fp, v[0], v[1], v[2])


p1=shapes.Point(pos=[8,1.1,1.1])
p2=shapes.Point(pos=[1.1,1.1,1.1])
p3=shapes.Point(pos=[5,5,1.1])
#system.constraints.add(particle_type=1, penetrable=1, shape=p1)
#system.constraints.add(particle_type=1, penetrable=1, shape=p2)
#system.constraints.add(particle_type=1, penetrable=1, shape=p3)

#s1=shapes.Segment(a=[5,5,1.1], b=[9,1.1,1.1])
#system.constraints.add(particle_type=1, penetrable=1, shape=s1)

#tri=shapes.Triangle(a=[8,1.1,1.1], b=[1.1,1.1,1.1], c=[5,5,1.1])
#tri=shapes.Triangle(a=[1,1.1,1.1], b=[8.1,1.1,1.1], c=[5,5,1.1])
a=[5,5,1.1]
b=[8.1,1.1,1.1]
c=[1.,1.1,1.1]

for v in vert:
  a=v[0]
  b=v[1]
  c=v[2]
#  tri=shapes.Triangle(a=[v[0,0], v[0,1], v[0,2]], b=[v[1,0], v[1,1], v[1,2]], c=[v[2,0], v[2,1], v[2,2]])
  tri=shapes.Triangle(a=[float(a[0]), float(a[1]), float(a[2])], b=[float(b[0]), float(b[1]), float(b[2])], c=[float(c[0]), float(c[1]), float(c[2])])
  system.constraints.add(particle_type=1, penetrable=1, shape=tri)
  print("added")

  #draw_triangle(fp, v[0], v[1], v[2])


exit()
for c in system.constraints.get():
    print(c)
    
#draw_triangle(fp, a, b, c)

fp.close()


def main():

    while True: 
        system.integrator.run(1)
        #Update particle information safely here
        visualizer.update()

##Start simulation in seperate thread
#t = Thread(target=main)
#t.daemon = True
#t.start()


p=0
z=8.0
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


for t in range(2000):
  system.integrator.run(10)
  vtf_file.write("\ntimestep ordered\n")
  for p in system.part:
    vtf_file.write("{} {} {} \n".format(p.pos[0], p.pos[1], p.pos[2]))    


exit()
#print('{}\n'.format(energies['non_bonded']))







exit()

z=0.00
fp=open('phiTriangle.dat', 'w')
for x in np.arange(0,system.box_l[0],0.1):
  for y in np.arange(0,system.box_l[1],0.1):
    system.part[0].pos=[x,y,z]
    energies = system.analysis.energy()
    fp.write('{}\n'.format(energies['non_bonded']))
fp.close()

#system.integrator.run(1)
exit()

system.part[0].pos=[3,1,0]
system.integrator.run(1)

system.part[0].pos=[6,1,0]
system.integrator.run(1)

system.part[0].pos=[12,1,0]
system.integrator.run(1)
#system.part[0].pos=[.4,0,0]


   


exit()

