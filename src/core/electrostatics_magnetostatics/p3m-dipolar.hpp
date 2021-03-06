/*
  Copyright (C) 2010-2018 The ESPResSo project
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
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _P3M_MAGNETOSTATICS_H
#define _P3M_MAGNETOSTATICS_H
/** \file
 * P3M algorithm for long range magnetic dipole-dipole
 * interaction.
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the dipolar Ewald summation. Details of the used method can be found in
 *  Hockney/Eastwood and Deserno/Holm. The file p3m contains only the
 *  Particle-Mesh part.
 *
 *  Further reading:
 *  <ul>
 *  <li> J.J. Cerda, P3M for dipolar interactions. J. Chem. Phys, 129, xxx
 * ,(2008).
 *  </ul>
 *
 */

#include "config.hpp"

#ifdef DP3M
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "p3m-common.hpp"

#include "utils/math/AS_erfc_part.hpp"

typedef struct {
  p3m_parameter_struct params;

  /** local mesh. */
  p3m_local_mesh local_mesh;
  /** real space mesh (local) for CA/FFT.*/
  double *rs_mesh;
  /** real space mesh (local) for CA/FFT of the dipolar field.*/
  double *rs_mesh_dip[3];
  /** k space mesh (local) for k space calculation and FFT.*/
  double *ks_mesh;

  /** number of dipolar particles (only on master node). */
  int sum_dip_part;
  /** Sum of square of magnetic dipoles (only on master node). */
  double sum_mu2;

  /** interpolation of the charge assignment function. */
  double *int_caf[7];

  /** position shift for calc. of first assignment mesh point. */
  double pos_shift;
  /** help variable for calculation of aliasing sums */
  double *meshift;

  /** Spatial differential operator in k-space. We use an i*k differentiation.
   */
  double *d_op;
  /** Force optimised influence function (k-space) */
  double *g_force;
  /** Energy optimised influence function (k-space) */
  double *g_energy;

  /** number of charged particles on the node. */
  int ca_num;

  /** Charge fractions for mesh assignment. */
  double *ca_frac;
  /** index of first mesh point for charge assignment. */
  int *ca_fmp;
  /** number of permutations in k_space */
  int ks_pnum;

  /** send/recv mesh sizes */
  p3m_send_mesh sm;

  /** Field to store grid points to send. */
  double *send_grid;
  /** Field to store grid points to recv */
  double *recv_grid;

  /* Stores the value of the energy correction due to MS effects */
  double energy_correction;
} dp3m_data_struct;

/** dipolar P3M parameters. */
extern dp3m_data_struct dp3m;

/** \name Exported Functions */
/************************************************************/
/*@{*/

void dp3m_pre_init();

void dp3m_set_tune_params(double r_cut, int mesh, int cao, double alpha,
                          double accuracy, int n_interpol);

int dp3m_set_params(double r_cut, int mesh, int cao, double alpha,
                    double accuracy);

int dp3m_set_ninterpol(int n);

int dp3m_set_mesh_offset(double x, double y, double z);

int dp3m_set_eps(double eps);

/** Initialize all structures, parameters and arrays needed for the
 *  P3M algorithm for dipole-dipole interactions.
 */
void dp3m_init(void);

/** Updates \ref p3m_parameter_struct::alpha and \ref
 * p3m_parameter_struct::r_cut if \ref box_l changed. */
void dp3m_scaleby_box_l();

/// sanity checks
int dp3m_sanity_checks();

/** assign the physical dipoles using the tabulated assignment function.
    If Dstore_ca_frac is true, then the charge fractions are buffered in
   Dcur_ca_fmp and
    Dcur_ca_frac. */
void dp3m_dipole_assign(void);

/** reset dipolar p3m core parameters */
void dp3m_deactivate(void);

int dp3m_adaptive_tune(char **log);

/** compute the k-space part of forces and energies for the magnetic
 * dipole-dipole interaction  */
double dp3m_calc_kspace_forces(int force_flag, int energy_flag);

/** Calculate number of magnetic  particles, the sum of the squared
    charges and the squared sum of the charges. */

void dp3m_count_magnetic_particles();

/** assign a single dipole into the current charge grid. cp_cnt gives the a
   running index,
    which may be smaller than 0, in which case the charge is assumed to be
   virtual and is not
    stored in the Dca_frac arrays. */
void dp3m_assign_dipole(double const real_pos[3], double mu,
                        double const dip[3], int cp_cnt);

/** shrink wrap the dipoles grid */
void dp3m_shrink_wrap_dipole_grid(int n_dipoles);

/** Calculate real space contribution of p3m dipolar pair forces and torques.
    If NPT is compiled in, it returns the energy, which is needed for NPT. */
inline double dp3m_add_pair_force(Particle *p1, Particle *p2, double *d,
                                  double dist2, double dist, double force[3]) {
  if ((p1->p.dipm == 0.) || (p2->p.dipm == 0.))
    return 0.;

  double coeff, exp_adist2;
  const Vector3d dip1 = p1->calc_dip();
  const Vector3d dip2 = p2->calc_dip();
  double B_r, C_r, D_r;
  double alpsq = dp3m.params.alpha * dp3m.params.alpha;
#ifdef ROTATION
  double mixmj[3], mixr[3], mjxr[3];
#endif

  if (dist < dp3m.params.r_cut && dist > 0) {
    double adist = dp3m.params.alpha * dist;
#if USE_ERFC_APPROXIMATION
    double erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
#else
    double erfc_part_ri = erfc(adist) / dist;
#endif

    // Calculate scalar multiplications for vectors mi, mj, rij
    double mimj = dip1 * dip2;

    double mir = dip1 * Vector3d{d[0], d[1], d[2]};
    double mjr = dip2 * Vector3d{d[0], d[1], d[2]};

    coeff = 2.0 * dp3m.params.alpha * wupii;
    double dist2i = 1 / dist2;
    exp_adist2 = exp(-adist * adist);

    if (dp3m.params.accuracy > 5e-06)
      B_r = (erfc_part_ri + coeff) * exp_adist2 * dist2i;
    else
      B_r = (erfc(adist) / dist + coeff * exp_adist2) * dist2i;

    C_r = (3 * B_r + 2 * alpsq * coeff * exp_adist2) * dist2i;
    D_r = (5 * C_r + 4 * coeff * alpsq * alpsq * exp_adist2) * dist2i;

    // Calculate real-space forces
    for (int j = 0; j < 3; j++)
      force[j] += coulomb.Dprefactor *
                  ((mimj * d[j] + dip1[j] * mjr + dip2[j] * mir) * C_r -
                   mir * mjr * D_r * d[j]);

// Calculate vector multiplications for vectors mi, mj, rij
#ifdef ROTATION
    mixmj[0] = dip1[1] * dip2[2] - dip1[2] * dip2[1];
    mixmj[1] = dip1[2] * dip2[0] - dip1[0] * dip2[2];
    mixmj[2] = dip1[0] * dip2[1] - dip1[1] * dip2[0];

    mixr[0] = dip1[1] * d[2] - dip1[2] * d[1];
    mixr[1] = dip1[2] * d[0] - dip1[0] * d[2];
    mixr[2] = dip1[0] * d[1] - dip1[1] * d[0];

    mjxr[0] = dip2[1] * d[2] - dip2[2] * d[1];
    mjxr[1] = dip2[2] * d[0] - dip2[0] * d[2];
    mjxr[2] = dip2[0] * d[1] - dip2[1] * d[0];

    // Calculate real-space torques
    for (int j = 0; j < 3; j++) {
      p1->f.torque[j] +=
          coulomb.Dprefactor * (-mixmj[j] * B_r + mixr[j] * mjr * C_r);
      p2->f.torque[j] +=
          coulomb.Dprefactor * (mixmj[j] * B_r + mjxr[j] * mir * C_r);
    }
#endif
#ifdef NPT
#if USE_ERFC_APPROXIMATION
    double fac1 =
        coulomb.Dprefactor * p1->p.dipm * p2->p.dipm * exp(-adist * adist);
#else
    double fac1 = coulomb.Dprefactor * p1->p.dipm * p2->p.dipm;
#endif
    return fac1 * (mimj * B_r - mir * mjr * C_r);
#endif
  }
  return 0.0;
}

/** Calculate real space contribution of dipolar pair energy. */
inline double dp3m_pair_energy(Particle *p1, Particle *p2, double *d,
                               double dist2, double dist) {
  const Vector3d dip1 = p1->calc_dip();
  const Vector3d dip2 = p2->calc_dip();
  double /* fac1,*/ adist, erfc_part_ri, coeff, exp_adist2, dist2i;
  double mimj, mir, mjr;
  double B_r, C_r;
  double alpsq = dp3m.params.alpha * dp3m.params.alpha;

  if (dist < dp3m.params.r_cut && dist > 0) {
    adist = dp3m.params.alpha * dist;
    /*fac1 = coulomb.Dprefactor;*/

#if USE_ERFC_APPROXIMATION
    erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
    /*  fac1 = coulomb.Dprefactor * p1->p.dipm*p2->p.dipm; IT WAS WRONG */ /* *exp(-adist*adist);
                                                                            */
#else
    erfc_part_ri = erfc(adist) / dist;
/* fac1 = coulomb.Dprefactor * p1->p.dipm*p2->p.dipm;  IT WAS WRONG*/
#endif

    // Calculate scalar multiplications for vectors mi, mj, rij
    mimj = dip1 * dip2;
    mir = dip1 * Vector3d{d[0], d[1], d[2]};
    mjr = dip2 * Vector3d{d[0], d[1], d[2]};

    coeff = 2.0 * dp3m.params.alpha * wupii;
    dist2i = 1 / dist2;
    exp_adist2 = exp(-adist * adist);

    if (dp3m.params.accuracy > 5e-06)
      B_r = (erfc_part_ri + coeff) * exp_adist2 * dist2i;
    else
      B_r = (erfc(adist) / dist + coeff * exp_adist2) * dist2i;

    C_r = (3 * B_r + 2 * alpsq * coeff * exp_adist2) * dist2i;

    /*
      printf("(%4i %4i) pair energy = %f (B_r=%15.12f
      C_r=%15.12f)\n",p1->p.identity,p2->p.identity,fac1*(mimj*B_r-mir*mjr*C_r),B_r,C_r);
    */

    /* old line return fac1 * ( mimj*B_r - mir*mjr * C_r );*/
    return coulomb.Dprefactor * (mimj * B_r - mir * mjr * C_r);
  }
  return 0.0;
}

#endif /* DP3M */
#endif /* _P3M_DIPOLES_H */
