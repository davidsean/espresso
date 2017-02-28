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


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


import numpy as np

print("""
=======================================================
=                      test.py                         =
=======================================================

Program Information:""")
print(code_info.features())


system = espressomd.System()


system.box_l =[6,6,6]
system.time_step = 0.01
system.thermostat.set_langevin(kT=0.0, gamma=1.0)

system.non_bonded_inter[0, 0].lennard_jones.set_params(
     epsilon=1, sigma=1,
     cutoff=2**(1. / 6), shift="auto")



#floor=shapes.Wall(normal=[5,4,5], dist=8.124)
#system.constraints.add(particle_type=0, penetrable=1, shape=floor)
floor=shapes.Wall(normal=[2.5,2.5,2.5], dist=4.3301270189221936)
system.constraints.add(particle_type=0, penetrable=1, shape=floor)


fig1 = plt.figure(figsize=plt.figaspect(1.0))
ax1 = fig1.add_subplot(111,projection='3d')


xpos=[]
ypos=[]
zpos=[]

z=0.00
#fp=open('phiTriangle.dat', 'w')
for x in np.arange(0,system.box_l[0],0.2001):
  for y in np.arange(0,system.box_l[1],0.2001):
    for z in np.arange(0,system.box_l[2],0.2001):

      system.part[0].pos=[x,y,z]
      energies = system.analysis.energy()
      if (energies['non_bonded'])>1:
        xpos.append(x)
        ypos.append(y)
        zpos.append(z)

alpha=0.2
ax1.plot3D(xpos, ypos, zpos, ls='', marker='.', color='black',alpha=alpha)
ax1.set_xlim3d(0,system.box_l[0])
ax1.set_ylim3d(0,system.box_l[1])
ax1.set_zlim3d(0,system.box_l[2])
fig1.savefig("test.png")
plt.show()
plt.close(fig1)


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

