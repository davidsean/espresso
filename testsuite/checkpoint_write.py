#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
import numpy as np
import unittest as ut
from tests_common import abspath
import espressomd
from espressomd import checkpointing


@ut.skipIf(not espressomd.has_features(["LENNARD_JONES"]) ,
           "Features not available, skipping test!")

class CheckpointTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    checkpoint = checkpointing.Checkpointing(checkpoint_id="test_checkpoint")


    def test(self):
        self.system.seed=1234
        np.random.seed(seed=self.system.seed)
        self.system.box_l = [100.0, 100.0, 100.0]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1./6.), shift="auto")
        num_part=100
        self.system.part.add(id=np.arange(num_part,dtype=int), pos=np.random.random((num_part,3))*self.system.box_l, type=np.zeros(num_part,dtype=int))
        self.system.thermostat.set_langevin(kT=1.0, gamma=1.0)
        
        i = 0
        lj_cap = 10
        min_dist = self.system.analysis.min_dist()
        while (i < 50 and min_dist < 0.9):
            system.integrator.run(100)
            min_dist = self.system.analysis.min_dist()
            i += 1
            lj_cap = lj_cap + 10
            self.system.force_cap = lj_cap
        self.system.force_cap = 0
        self.system.integrator.run(100)

        self.checkpoint.register("CheckpointTest.system")
        self.checkpoint.register("CheckpointTest.system.box_l")
        self.checkpoint.register("CheckpointTest.system.cell_system")
        self.checkpoint.register("CheckpointTest.system.time_step")
        self.checkpoint.register("CheckpointTest.system.non_bonded_inter")
        self.checkpoint.register("CheckpointTest.system.part")
        self.checkpoint.register("CheckpointTest.system.thermostat")
        self.checkpoint.register("CheckpointTest.system.integrator")
        self.checkpoint.register("CheckpointTest.system.random_number_generator_state")
        self.checkpoint.save()
        #self.system.integrator.run(100)
        old_pos = self.system.part[:].pos.copy()
        old_f = self.system.part[:].f.copy()
        old_v = self.system.part[:].v.copy()
        np.savetxt("test_checkpoint.txt", np.transpose([old_pos[:,0], old_pos[:,1], old_pos[:,2],
                                                        old_f[:,0], old_f[:,1], old_f[:,2],
                                                        old_v[:,0], old_v[:,1], old_v[:,2]]),  fmt='%.6f')
        
if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()








