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
system.thermostat.set_langevin(kT=0.1, gamma=1.0)

system.non_bonded_inter[0, 0].lennard_jones.set_params(
     epsilon=1, sigma=1,
     cutoff=2**(1. / 6), shift="auto")
#system.cell_system.set_n_square(use_verlet_lists=False)


#system.part.add(id=0, pos=[box_l/2.0,box_l/2.0,box_l/2.0], fix=[1,1,1])
system.part.add(id=0, pos=[2,1,1],type=0, ext_force=[0,0,0])

floor=shapes.Wall(normal=[1,0,0], dist=0.001)
c1=system.constraints.add(particle_type=0, penetrable=1, shape=floor)

#  Point p1 = Point({1.0, 0.0, 0.0});
#  Point p1 = Point({5.0, 0.0, 0.0});

energies = system.analysis.energy()
#print('{}\n'.format(energies['non_bonded']))


def main():

    while True: 
        system.integrator.run(1)
        #Update particle information safely here
        visualizer.update()

#Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

visualizer.update()
visualizer.start()

exit()

z=0.1
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

