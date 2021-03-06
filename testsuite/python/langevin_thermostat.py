#
# Copyright (C) 2013-2018 The ESPResSo project
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from espressomd.interactions import FeneBond
from time import time
from espressomd.accumulators import Correlator
from espressomd.observables import ParticleVelocities, ParticleBodyAngularVelocities
from tests_common import single_component_maxwell


@ut.skipIf(espressomd.has_features("THERMOSTAT_IGNORE_NON_VIRTUAL"),
           "Skipped because of THERMOSTAT_IGNORE_NON_VIRTUAL")
class LangevinThermostat(ut.TestCase):

    """Tests the velocity distribution created by the Langevin thermostat against
       the single component Maxwell distribution."""

    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.cell_system.set_domain_decomposition(use_verlet_lists=True)
    system.cell_system.skin = 0
    system.seed = range(system.cell_system.get_state()["n_nodes"])
    if espressomd.has_features("PARTIAL_PERIODIC"):
        system.periodicity = 0, 0, 0

    @classmethod
    def setUpClass(cls):
        np.random.seed(42)

    def check_velocity_distribution(self, vel, minmax, n_bins, error_tol, kT):
        """check the recorded particle distributions in vel againsta histogram with n_bins bins. Drop velocities outside minmax. Check individual histogram bins up to an accuracy of error_tol agaisnt the analytical result for kT."""
        for i in range(3):
            hist = np.histogram(
                vel[:, i], range=(-minmax, minmax), bins=n_bins, normed=False)
            data = hist[0] / float(vel.shape[0])
            bins = hist[1]
            for j in range(n_bins):
                found = data[j]
                expected = single_component_maxwell(
                    bins[j], bins[j + 1], kT)
                self.assertLessEqual(abs(found - expected), error_tol)

    def test_aa_verify_single_component_maxwell(self):
        """Verifies the normalization of the analytical expression."""
        self.assertLessEqual(
            abs(single_component_maxwell(-10, 10, 4.) - 1.), 1E-4)

    def test_global_langevin(self):
        """Test for global Langevin parameters."""
        N = 200
        system = self.system
        system.part.clear()
        system.time_step = 0.06

        # Place particles
        system.part.add(pos=np.random.random((N, 3)))

        # Enable rotation if compiled in
        if espressomd.has_features("ROTATION"):
            system.part[:].rotation = 1, 1, 1

        kT = 1.1
        gamma = 3.5
        system.thermostat.set_langevin(kT=kT, gamma=gamma)

        # Warmup
        system.integrator.run(100)

        # Sampling
        loops = 400
        v_stored = np.zeros((N * loops, 3))
        omega_stored = np.zeros((N * loops, 3))
        for i in range(loops):
            system.integrator.run(1)
            v_stored[i * N:(i + 1) * N, :] = system.part[:].v
            if espressomd.has_features("ROTATION"):
                omega_stored[i * N:(i + 1) * N, :] = system.part[:].omega_body

        v_minmax = 5
        bins = 4
        error_tol = 0.016
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, kT)
        if espressomd.has_features("ROTATION"):
            self.check_velocity_distribution(
                omega_stored, v_minmax, bins, error_tol, kT)

    @ut.skipIf(not espressomd.has_features("LANGEVIN_PER_PARTICLE"),
               "Test requires LANGEVIN_PER_PARTICLE")
    def test_langevin_per_particle(self):
        """Test for Langevin particle. Covers all combinations of
           particle specific gamma and temp set or not set.
        """
        N = 400
        system = self.system
        system.part.clear()
        system.time_step = 0.06
        system.part.add(pos=np.random.random((N, 3)))
        if espressomd.has_features("ROTATION"):
            system.part[:].rotation = 1, 1, 1

        kT = 0.9
        gamma = 3.2
        gamma2 = 4.3
        kT2 = 1.5
        system.thermostat.set_langevin(kT=kT, gamma=gamma)
        # Set different kT on 2nd half of particles
        system.part[int(N / 2):].temp = kT2
        # Set different gamma on half of the partiles (overlap over both kTs)
        if espressomd.has_features("PARTICLE_ANISOTROPY"):
            system.part[
                int(N / 4):int(3 * N / 4)].gamma = gamma2, gamma2, gamma2
        else:
            system.part[int(N / 4):int(3 * N / 4)].gamma = gamma2

        system.integrator.run(50)
        loops = 600

        v_kT = np.zeros((int(N / 2) * loops, 3))
        v_kT2 = np.zeros((int(N / 2 * loops), 3))

        if espressomd.has_features("ROTATION"):
            omega_kT = np.zeros((int(N / 2) * loops, 3))
            omega_kT2 = np.zeros((int(N / 2 * loops), 3))

        for i in range(loops):
            system.integrator.run(1)
            v_kT[int(i * N / 2):int((i + 1) * N / 2),
                 :] = system.part[:int(N / 2)].v
            v_kT2[int(i * N / 2):int((i + 1) * N / 2),
                  :] = system.part[int(N / 2):].v

            if espressomd.has_features("ROTATION"):
                omega_kT[int(i * N / 2):int((i + 1) * N / 2),
                         :] = system.part[:int(N / 2)].omega_body
                omega_kT2[int(i * N / 2):int((i + 1) * N / 2),
                          :] = system.part[int(N / 2):].omega_body
        v_minmax = 5
        bins = 4
        error_tol = 0.016
        self.check_velocity_distribution(v_kT, v_minmax, bins, error_tol, kT)
        self.check_velocity_distribution(v_kT2, v_minmax, bins, error_tol, kT2)

        if espressomd.has_features("ROTATION"):
            self.check_velocity_distribution(
                omega_kT, v_minmax, bins, error_tol, kT)
            self.check_velocity_distribution(
                omega_kT2, v_minmax, bins, error_tol, kT2)

    def setup_diff_mass_rinertia(self, p):
        if espressomd.has_features("MASS"):
            p.mass = 0.5
        if espressomd.has_features("ROTATION"):
            p.rotation = 1, 1, 1
            # Make sure rinertia does not change diff coeff
            if espressomd.has_features("ROTATIONAL_INERTIA"):
                p.rinertia = 0.4, 0.4, 0.4

    def test_diffusion(self):
        """This tests rotational and translational diffusion coeff via green-kubo"""
        system = self.system
        system.part.clear()

        kT = 1.37
        dt = 0.1
        system.time_step = dt

        # Translational gamma. We cannot test per-component, if rotation is on,
        # because body and space frames become different.
        gamma = 3.1

        # Rotational gamma
        gamma_rot_i = 4.7
        gamma_rot_a = 4.2, 1, 1.2

        # If we have langevin per particle:
        # per particle kT
        per_part_kT = 1.6
        # Translation
        per_part_gamma = 1.63
        # Rotational
        per_part_gamma_rot_i = 2.6
        per_part_gamma_rot_a = 2.4, 3.8, 1.1

        # Particle with global thermostat params
        p_global = system.part.add(pos=(0, 0, 0))
        # Make sure, mass doesn't change diff coeff
        self.setup_diff_mass_rinertia(p_global)

        # particle specific gamma, kT, and both
        if espressomd.has_features("LANGEVIN_PER_PARTICLE"):
            p_gamma = system.part.add(pos=(0, 0, 0))
            self.setup_diff_mass_rinertia(p_gamma)
            if espressomd.has_features("PARTICLE_ANISOTROPY"):
                p_gamma.gamma = per_part_gamma, per_part_gamma, per_part_gamma
                if espressomd.has_features("ROTATION"):
                    p_gamma.gamma_rot = per_part_gamma_rot_a
            else:
                p_gamma.gamma = per_part_gamma
                if espressomd.has_features("ROTATION"):
                    p_gamma.gamma_rot = per_part_gamma_rot_i

            p_kT = system.part.add(pos=(0, 0, 0))
            self.setup_diff_mass_rinertia(p_kT)
            p_kT.temp = per_part_kT

            p_both = system.part.add(pos=(0, 0, 0))
            self.setup_diff_mass_rinertia(p_both)
            p_both.temp = per_part_kT
            if espressomd.has_features("PARTICLE_ANISOTROPY"):
                p_both.gamma = per_part_gamma, per_part_gamma, per_part_gamma
                if espressomd.has_features("ROTATION"):
                    p_both.gamma_rot = per_part_gamma_rot_a
            else:
                p_both.gamma = per_part_gamma
                if espressomd.has_features("ROTATION"):
                    p_both.gamma_rot = per_part_gamma_rot_i

        # Thermostat setup
        if espressomd.has_features("ROTATION"):
            if espressomd.has_features("PARTICLE_ANISOTROPY"):
                # particle anisotropy and rotation
                system.thermostat.set_langevin(
                    kT=kT, gamma=gamma, gamma_rotation=gamma_rot_a)
            else:
                # Rotation without particle anisotropy
                system.thermostat.set_langevin(
                    kT=kT, gamma=gamma, gamma_rotation=gamma_rot_i)
        else:
            # No rotation
            system.thermostat.set_langevin(kT=kT, gamma=gamma)

        system.cell_system.skin = 0.4
        system.integrator.run(100)

        # Correlators
        vel_obs = {}
        omega_obs = {}
        corr_vel = {}
        corr_omega = {}
        all_particles = [p_global]
        if espressomd.has_features("LANGEVIN_PER_PARTICLE"):
            all_particles.append(p_gamma)
            all_particles.append(p_kT)
            all_particles.append(p_both)

        # linear vel
        vel_obs = ParticleVelocities(ids=system.part[:].id)
        corr_vel = Correlator(obs1=vel_obs, tau_lin=20, tau_max=1.4, delta_N=1,
                              corr_operation="componentwise_product", compress1="discard1")
        system.auto_update_accumulators.add(corr_vel)
        # angular vel
        if espressomd.has_features("ROTATION"):
            omega_obs = ParticleBodyAngularVelocities(ids=system.part[:].id)
            corr_omega = Correlator(
                obs1=omega_obs, tau_lin=20, tau_max=1.5, delta_N=1,
                                    corr_operation="componentwise_product", compress1="discard1")
            system.auto_update_accumulators.add(corr_omega)

        system.integrator.run(150000)

        system.auto_update_accumulators.remove(corr_vel)
        corr_vel.finalize()
        if espressomd.has_features("ROTATION"):
            system.auto_update_accumulators.remove(corr_omega)
            corr_omega.finalize()

        # Verify diffusion
        # Translation
        # Cast gammas to vector, to make checks independent of
        # PARTICLE_ANISOTROPY
        gamma = np.ones(3) * gamma
        per_part_gamma = np.ones(3) * per_part_gamma
        self.verify_diffusion(p_global, corr_vel, kT, gamma)
        if espressomd.has_features("LANGEVIN_PER_PARTICLE"):
            self.verify_diffusion(p_gamma, corr_vel, kT, per_part_gamma)
            self.verify_diffusion(p_kT, corr_vel, per_part_kT, gamma)
            self.verify_diffusion(
                p_both, corr_vel, per_part_kT, per_part_gamma)

        # Rotation
        if espressomd.has_features("ROTATION"):
            # Decide on effective gamma rotation, since for rotation it is
            # direction dependent
            eff_gamma_rot = None
            per_part_eff_gamma_rot = None
            if espressomd.has_features("PARTICLE_ANISOTROPY"):
                eff_gamma_rot = gamma_rot_a
                eff_per_part_gamma_rot = per_part_gamma_rot_a
            else:
                eff_gamma_rot = gamma_rot_i * np.ones(3)
                eff_per_part_gamma_rot = per_part_gamma_rot_i * np.ones(3)

            self.verify_diffusion(p_global, corr_omega, kT, eff_gamma_rot)
            if espressomd.has_features("LANGEVIN_PER_PARTICLE"):
                self.verify_diffusion(
                    p_gamma, corr_omega, kT, eff_per_part_gamma_rot)
                self.verify_diffusion(
                    p_kT, corr_omega, per_part_kT, eff_gamma_rot)
                self.verify_diffusion(p_both, corr_omega,
                                      per_part_kT, eff_per_part_gamma_rot)

    def verify_diffusion(self, p, corr, kT, gamma):
        """Verifify diffusion coeff.

           p: particle, corr: dict containing correltor with particle as key,
           kT=kT, gamma=gamma as 3 component vector.
        """
        c = corr
        # Integral of vacf via Green-Kubo
        # D= int_0^infty <v(t_0)v(t_0+t)> dt     (o 1/3, since we work
        # componentwise)
        i = p.id
        acf = c.result()[:, [0, 2 + 3 * i, 2 + 3 * i + 1, 2 + 3 * i + 2]]
        np.savetxt("acf.dat", acf)

        # Integrate w. trapez rule
        for coord in 1, 2, 3:
            I = np.trapz(acf[:, coord], acf[:, 0])
            ratio = I / (kT / gamma[coord - 1])
            self.assertAlmostEqual(ratio, 1., delta=0.07)

    def test_00__friction_trans(self):
        """Tests the translational friction-only part of the thermostat."""

        system = self.system
        # Translation
        gamma_t_i = 2
        gamma_t_a = 0.5, 2, 1.5
        v0 = 5.

        system.time_step = 0.0005
        system.part.clear()
        system.part.add(pos=(0, 0, 0), v=(v0, v0, v0))
        if espressomd.has_features("MASS"):
            system.part[0].mass = 3
        if espressomd.has_features("PARTICLE_ANISOTROPY"):
            system.thermostat.set_langevin(kT=0, gamma=gamma_t_a)
        else:
            system.thermostat.set_langevin(kT=0, gamma=gamma_t_i)

        system.time = 0
        for i in range(100):
            system.integrator.run(10)
            for j in range(3):
                if espressomd.has_features("PARTICLE_ANISOTROPY"):
                    self.assertAlmostEqual(
                        system.part[0].v[j], v0 * np.exp(-gamma_t_a[j] / system.part[0].mass * system.time), places=2)
                else:
                    self.assertAlmostEqual(
                        system.part[0].v[j], v0 * np.exp(-gamma_t_i / system.part[0].mass * system.time), places=2)

    @ut.skipIf(not espressomd.has_features("ROTATION"), "Skipped for lack of ROTATION")
    def test_00__friction_rot(self):
        """Tests the rotational friction-only part of the thermostat."""

        system = self.system
        # Translation
        gamma_t_i = 2
        gamma_t_a = 0.5, 2, 1.5
        gamma_r_i = 3
        gamma_r_a = 1.5, 0.7, 1.2
        o0 = 5.

        system.time_step = 0.0005
        system.part.clear()
        system.part.add(
            pos=(0, 0, 0), omega_body=(o0, o0, o0), rotation=(1, 1, 1))
        if espressomd.has_features("ROTATIONAL_INERTIA"):
            system.part[0].rinertia = 2, 2, 2
        if espressomd.has_features("PARTICLE_ANISOTROPY"):
            system.thermostat.set_langevin(
                kT=0, gamma=gamma_t_a, gamma_rotation=gamma_r_a)
        else:
            system.thermostat.set_langevin(
                kT=0, gamma=gamma_t_i, gamma_rotation=gamma_r_i)

        system.time = 0
        for i in range(100):
            system.integrator.run(10)
            if espressomd.has_features("ROTATIONAL_INERTIA"):
                rinertia = system.part[0].rinertia
            else:
                rinertia = (1, 1, 1)
            for j in range(3):
                if espressomd.has_features("PARTICLE_ANISOTROPY"):
                    self.assertAlmostEqual(
                        system.part[0].omega_body[j], o0 * np.exp(-gamma_r_a[j] / rinertia[j] * system.time), places=2)
                else:
                    self.assertAlmostEqual(
                        system.part[0].omega_body[j], o0 * np.exp(-gamma_r_i / rinertia[j] * system.time), places=2)

    @ut.skipIf(not espressomd.has_features("VIRTUAL_SITES"), "Skipped for lack of VIRTUAL_SITES")
    def test_virtual(self):
        system = self.system
        system.time_step = 0.01
        system.part.clear()

        virtual = system.part.add(pos=[0, 0, 0], virtual=True, v=[1, 0, 0])
        physical = system.part.add(pos=[0, 0, 0], virtual=False, v=[1, 0, 0])

        system.thermostat.set_langevin(
            kT=0, gamma=1, gamma_rotation=1., act_on_virtual=False)

        system.integrator.run(0)

        np.testing.assert_almost_equal(np.copy(virtual.f), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.f), [-1, 0, 0])

        system.thermostat.set_langevin(
            kT=0, gamma=1, gamma_rotation=1., act_on_virtual=True)
        system.integrator.run(0)

        np.testing.assert_almost_equal(np.copy(virtual.f), [-1, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.f), [-1, 0, 0])


if __name__ == "__main__":
    ut.main()
