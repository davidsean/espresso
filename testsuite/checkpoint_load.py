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
from espressomd import thermostat
from espressomd import checkpointing
import shutil
import os


@ut.skipIf(not espressomd.has_features(["LENNARD_JONES"]) ,
           "Features not available, skipping test!")

class CheckpointTest(ut.TestCase):
    checkpoint = checkpointing.Checkpointing(checkpoint_id="test_checkpoint")
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    
    def tearDown(self):
        shutil.rmtree("test_checkpoint")
        os.remove("test_checkpoint.txt")
        
    def test(self):
        self.checkpoint.load()
        # load data 100 steps after checkpoint,
        data = np.loadtxt("test_checkpoint.txt")
        old_pos = data[:,0:3]
        old_f = data[:,3:6]
        old_v = data[:,6:9]
        # re-do the integration and hope we get the same results!
        #self.system.integrator.run(100)
        np.testing.assert_allclose(np.copy(self.system.part[:].pos), old_pos, atol=1e-5)
        np.testing.assert_allclose(np.copy(self.system.part[:].f), old_f, atol=1e-5)
        np.testing.assert_allclose(np.copy(self.system.part[:].v), old_v, atol=1e-5)

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()








