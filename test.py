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

system.box_l =[100,100,100]
system.time_step = 0.01
system.thermostat.set_langevin(kT=0.0, gamma=1.0)

system.non_bonded_inter[0, 0].lennard_jones.set_params(
     epsilon=1, sigma=1,
     cutoff=2**(1. / 6), shift="auto")

fene = interactions.FeneBond(k=30., d_r_max=1.5)
system.bonded_inter.add(fene)

#system.cell_system.set_n_square(use_verlet_lists=False)


#system.part.add(id=0, pos=[box_l/2.0,box_l/2.0,box_l/2.0], fix=[1,1,1])
#system.part.add(id=0, pos=[2,1,1],type=0, ext_force=[0,0,0])



## ut a pillar of voxels

for layer in range(10):
  corner=np.array([0,0,0])+np.array([0,layer,0])
  dist= ( np.sqrt(np.sum((corner**2))) )
  print("floor=shapes.Wall(normal=[{},{},{}], dist={})".format(corner[0], corner[1], corner[2],dist))
  print("system.constraints.add(particle_type=0, penetrable=1, shape=floor)")

for y in np.arange(0.,9.,1.0):
  vec="[%f, %f, %f]" %(0,y,0)
  floor=shapes.Wall(normal=[0,y,0], dist=0.0)
#  system.constraints.add(particle_type=0, penetrable=1, shape=floor)

exit()

floor=shapes.Wall(normal=[0,0,2], dist=0.0)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)
floor=shapes.Wall(normal=[0,1,2], dist=1.0)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)
floor=shapes.Wall(normal=[0,2,2], dist=2.0)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)
floor=shapes.Wall(normal=[0,3,2], dist=3.0)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)
floor=shapes.Wall(normal=[0,4,2], dist=4.0)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)
floor=shapes.Wall(normal=[0,5,2], dist=5.0)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)
floor=shapes.Wall(normal=[0,6,2], dist=6.0)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)
floor=shapes.Wall(normal=[0,7,2], dist=7.0)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)
floor=shapes.Wall(normal=[0,8,2], dist=8.0)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)
floor=shapes.Wall(normal=[0,9,2], dist=9.0)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)

bond_length=0.97
n_poly=20
x=5
y=5.5
z=0.01
for p in range(n_poly):
  system.part.add(id=p, pos=[x,y,z], type=0, ext_force=[-0.1,0,0])
  z+=bond_length
for p in range(n_poly-1):
  system.part[p].add_bond((fene,p+1))
for p in system.part:
  p.pos[2]-=12.0

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

for t in range(1000):
  system.integrator.run(10)
  vtf_file.write("\ntimestep ordered\n")
  for p in system.part:
    vtf_file.write("{} {} {} \n".format(p.pos[0], p.pos[1], p.pos[2]))    


#  Point p1 = Point({1.0, 0.0, 0.0});
#  Point p1 = Point({5.0, 0.0, 0.0});

energies = system.analysis.energy()


exit()
#print('{}\n'.format(energies['non_bonded']))


def main():

    while True: 
        system.integrator.run(1)
        #Update particle information safely here
        visualizer.update()

#Start simulation in seperate thread
#t = Thread(target=main)
#t.daemon = True
#t.start()

#visualizer.update()
#visualizer.start()





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

