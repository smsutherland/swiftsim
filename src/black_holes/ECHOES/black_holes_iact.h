/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_ECHOES_BH_IACT_H
#define SWIFT_ECHOES_BH_IACT_H

#include "black_holes_part.h"
#include "black_holes_properties.h"
#include "cell.h"
#include "cosmology.h"
#include "entropy_floor.h"
#include "error.h"
#include "gravity_properties.h"
#include "kernel_gravity.h"
#include "kernel_hydro.h"

#include <math.h>

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas, not updated).
 * @param xpj The extended data of the second particle (not updated).
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct bpart *bi, const struct part *pj, const struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  bi->density.wcount += wi;
  bi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

#ifdef DEBUG_INTERACTIONS_BH
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_BH)
    bi->ids_ngbs_density[si->num_ngb_density] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_density;
#endif
}

/**
 * @brief Repositioning interaction between two particles (non-symmetric).
 *
 * Function used to identify the gas particle that this BH may move towards.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas)
 * @param xpj The extended data of the second particle.
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 * @param time Current physical time in the simulation.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_repos(
    const float r2, const float dx[3], const float hi, const float hj,
    struct bpart *bi, const struct part *pj, const struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time) {

  /* sutherland TODO: implement conditions on repositioning.
   * EAGLE has a velocity condition based on the sound speed, but we don't have
   * that sort of physics implemented yet.
   * EAGLE also has the option to negate the BH's own contribution to the
   * potential when determining the deepest gas particle. */

  const float max_dist_repos2 =
      kernel_gravity_softening_plummer_equivalent_inv *
      kernel_gravity_softening_plummer_equivalent_inv *
      bh_props->max_reposition_distance_ratio *
      bh_props->max_reposition_distance_ratio * grav_props->epsilon_baryon_cur *
      grav_props->epsilon_baryon_cur;

  /* Are we too far away? */
  if (r2 >= max_dist_repos2) return;

  if (bh_props->enable_repos_v_threshold) {
    /* Now check the velocity */
    const float delta_v[3] = {bi->v[0] - pj->v[0], bi->v[1] - pj->v[1],
                              bi->v[2] - pj->v[2]};
    const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                     delta_v[2] * delta_v[2];
    const float v2_pec = v2 * cosmo->a2_inv;
    const float v2_max = bh_props->repos_v2_threshold;
    if (v2_pec >= v2_max) return;
  }

  float potential = pj->black_holes_data.potential;

  /* Is the potential lower? */
  if (potential < bi->reposition.min_potential) {

    /* Store this as our new best */
    bi->reposition.min_potential = potential;
    bi->reposition.delta_x[0] = -dx[0];
    bi->reposition.delta_x[1] = -dx[1];
    bi->reposition.delta_x[2] = -dx[2];
  }
}

/**
 * @brief Swallowing interaction between two particles (non-symmetric).
 *
 * Function used to flag the gas particles that will be swallowed
 * by the black hole particle.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas)
 * @param xpj The extended data of the second particle.
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_swallow(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct bpart *bi, struct part *pj, struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time) {}

/**
 * @brief Swallowing interaction between two BH particles (non-symmetric).
 *
 * Function used to identify the BH particle that this BH may move towards.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param bj Second particle (black hole)
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_bh_repos(const float r2, const float dx[3],
                               const float hi, const float hj, struct bpart *bi,
                               const struct bpart *bj,
                               const struct cosmology *cosmo,
                               const struct gravity_props *grav_props,
                               const struct black_holes_props *bh_props,
                               const integertime_t ti_current) {

  /* sutherland TODO: implement conditions on repositioning.
   * EAGLE has a velocity condition based on the sound speed, but we don't have
   * that sort of physics implemented yet.
   * EAGLE also has the option to negate the BH's own contribution to the
   * potential when determining the deepest gas particle. */

  /* (Square of) Max repositioning distance allowed based on the softening */
  const float max_dist_repos2 =
      kernel_gravity_softening_plummer_equivalent_inv *
      kernel_gravity_softening_plummer_equivalent_inv *
      bh_props->max_reposition_distance_ratio *
      bh_props->max_reposition_distance_ratio * grav_props->epsilon_baryon_cur *
      grav_props->epsilon_baryon_cur;

  /* Are we too far away? */
  if (r2 >= max_dist_repos2) return;

  /* Now check the velocity */
  const float delta_v[3] = {bi->v[0] - bj->v[0], bi->v[1] - bj->v[1],
                            bi->v[2] - bj->v[2]};
  const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                   delta_v[2] * delta_v[2];
  const float v2_pec = v2 * cosmo->a2_inv;
  const float v2_max = bh_props->repos_v2_threshold;
  if (v2_pec >= v2_max) return;

  float potential = bj->reposition.potential;

  /* Is the potential lower? */
  if (potential < bi->reposition.min_potential) {

    /* Store this as our new best */
    bi->reposition.min_potential = potential;
    bi->reposition.delta_x[0] = -dx[0];
    bi->reposition.delta_x[1] = -dx[1];
    bi->reposition.delta_x[2] = -dx[2];
  }
}

/**
 * @brief Swallowing interaction between two BH particles (non-symmetric).
 *
 * Function used to flag the BH particles that will be swallowed
 * by the black hole particle.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param bj Second particle (black hole)
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_bh_swallow(const float r2, const float dx[3],
                                 const float hi, const float hj,
                                 const struct bpart *bi, struct bpart *bj,
                                 const struct cosmology *cosmo,
                                 const struct gravity_props *grav_props,
                                 const struct black_holes_props *bh_props,
                                 const integertime_t ti_current) {

  if (!bh_props->allow_intergroup_mergers) {
    /* ssutherland: Do we want to store the group_id in the bpart itself to
     * avoid this double indirection? */
    size_t bi_id = bi->gpart->fof_data.group_id;
    size_t bj_id = bj->gpart->fof_data.group_id;

    /* Require being in a group to swallow other BHs */
    if (bi_id == bh_props->group_id_default) return;

    /* Disallow swallowing when BHs are in different FoF groups.
     * If bj isn't in a fof group at all, it can still be swallowed. */
    if (bi_id != bj_id && bj_id != bh_props->group_id_default) {
      return;
    }
#ifdef SWIFT_DEBUG_CHECKS
    if (bi->fof_properties.is_central && bj->fof_properties.is_central) {
      error("BHs %lld and %lld are both in group %ld, but are both central!",
            bi->id, bj->id, bi_id);
    }
#endif
  }

  /* Disallow smaller BHs from swallowing larger ones.
   * In the case of mass ties, the BH with the larger ID swallows the other. */
  /* ssutherland NOTE: EAGLE has a similar condition on the subgrid mass which
   * we did not initially inherit.*/
  if (bi->fof_properties.max_group_mass < bj->fof_properties.max_group_mass ||
      (bi->fof_properties.max_group_mass == bj->fof_properties.max_group_mass &&
       bi->id < bj->id)) {
    return;
  }

  /* ssutherland NOTE:
   * EAGLE uses the mass of the more massive BH when calculating the escape
   * velocity. Is this necessary? I would think that we only check the
   * swallowing BH's mass.
   *
   * In spite of the non-symmetric nature of this function, it seems like SWIFT
   * considers bi swallowing bj to be roughly equivalent to bj swallowing bi.
   * Therefore, it makes sense to use the maximum mass here (I think). */

  /* Get useful constants */
  const float G_Newton = grav_props->G_Newton;

  /* Find the most massive of the two BHs */
  float M = bi->mass;
  if (bj->mass > M) {
    M = bj->mass;
  }

  /* Find the most massive of the two BHs */
  float M_halo = bi->fof_properties.max_group_mass;
  if (bj->fof_properties.max_group_mass > M_halo) {
    M_halo = bj->fof_properties.max_group_mass;
  }

  /* sutherland NOTE: we have a lot of duplicated code here.
   * Once the model is squared away, we can clean this up, optimize for the main
   * condition we use, and potentially even remove the other conditions
   * altogether. */

  int can_merge = 0;
  if (bh_props->merger_threshold_type == BH_mergers_escape_velocity) {

    /* Compute relative velocity */
    const float delta_v[3] = {
        bi->v[0] - bj->v[0],
        bi->v[1] - bj->v[1],
        bi->v[2] - bj->v[2],
    };
    /* |v|^2 */
    const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                     delta_v[2] * delta_v[2];
    /* Peculiar velocity.
     * Velocity in SWIFT is (v_pec * a) */
    const float v2_pec = v2 * cosmo->a2_inv;

    /* If v^2 is below this threshold, the BHs will merge */
    float v2_threshold;
    v2_threshold = 2.f * G_Newton * M / sqrt(r2);

    /* (Square of) max swallowing distance allowed based on the softening */
    const float max_dist_merge2 =
        kernel_gravity_softening_plummer_equivalent_inv *
        kernel_gravity_softening_plummer_equivalent_inv *
        bh_props->max_merging_distance_ratio *
        bh_props->max_merging_distance_ratio * grav_props->epsilon_baryon_cur *
        grav_props->epsilon_baryon_cur;

    can_merge = (v2_pec < v2_threshold) && (r2 <= max_dist_merge2);
  } else if (bh_props->merger_threshold_type == BH_mergers_2_escape_velocity) {

    /* Compute relative velocity */
    const float delta_v[3] = {
        bi->v[0] - bj->v[0],
        bi->v[1] - bj->v[1],
        bi->v[2] - bj->v[2],
    };
    /* |v|^2 */
    const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                     delta_v[2] * delta_v[2];
    /* Peculiar velocity.
     * Velocity in SWIFT is (v_pec * a) */
    const float v2_pec = v2 * cosmo->a2_inv;

    /* If v^2 is below this threshold, the BHs will merge */
    float v2_threshold;
    v2_threshold = 8.f * G_Newton * M / sqrt(r2);

    /* (Square of) max swallowing distance allowed based on the softening */
    const float max_dist_merge2 =
        kernel_gravity_softening_plummer_equivalent_inv *
        kernel_gravity_softening_plummer_equivalent_inv *
        bh_props->max_merging_distance_ratio *
        bh_props->max_merging_distance_ratio * grav_props->epsilon_baryon_cur *
        grav_props->epsilon_baryon_cur;

    can_merge = (v2_pec < v2_threshold) && (r2 <= max_dist_merge2);
  } else if (bh_props->merger_threshold_type == BH_mergers_kernel) {

    can_merge = 1;
  } else if (bh_props->merger_threshold_type == BH_mergers_virial) {

    /* Compute relative velocity */
    const float delta_v[3] = {
        bi->v[0] - bj->v[0],
        bi->v[1] - bj->v[1],
        bi->v[2] - bj->v[2],
    };
    /* |v|^2 */
    const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                     delta_v[2] * delta_v[2];
    /* Peculiar velocity.
     * Velocity in SWIFT is (v_pec * a) */
    const float v2_pec = v2 * cosmo->a2_inv;

    /* If v^2 is below this threshold, the BHs will merge */
    float v2_threshold;
    v2_threshold = 8.f * G_Newton * M / sqrt(r2);

    float virial_density_phys =
        cosmo->critical_density * cosmo->overdensity_BN98;

    /* sutherland NOTE: Is there a faster way to calculate this?
     * We could potentially make use of the custom icbrtf when its enabled. */
    float virial_radius_phys =
        cbrtf(0.75 * M_halo / (virial_density_phys * M_PI));

    float virial_radius2 =
        virial_radius_phys * virial_radius_phys * cosmo->a2_inv;

    can_merge = (v2_pec <= v2_threshold) && (r2 <= 0.15 * virial_radius2);
  } else {
    /* Cannot happen! */
#ifdef SWIFT_DEBUG_CHECKS
    error("Invalid choice of galaxy merger threshold type");
#endif
    can_merge = 0;
  }

  /* If they are close enough and the peculiar velocity is under the escape
   * velocity threshold. */
  if (can_merge) {
    /* This particle is swallowed by the BH with the largest mass of all the
     * candidates wanting to swallow it (we use IDs to break ties)*/
    if ((bj->merger_data.swallow_mass < bi->mass) ||
        (bj->merger_data.swallow_mass == bi->mass &&
         bj->merger_data.swallow_id < bi->id)) {

      message("BH %lld wants to swallow BH particle %lld", bi->id, bj->id);

      bj->merger_data.swallow_id = bi->id;
      bj->merger_data.swallow_mass = bi->mass;

    } else {

      message(
          "BH %lld wants to swallow gas particle %lld BUT CANNOT (old "
          "swallow id=%lld)",
          bi->id, bj->id, bj->merger_data.swallow_id);
    }
  }
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas)
 * @param xpj The extended data of the second particle.
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_feedback(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct bpart *bi, struct part *pj, struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time) {
  /* sutherland: note to DAA - This is where all our feedback prescriptions will
   * go. */

#ifdef DEBUG_INTERACTIONS_BH
  /* Update ngb counters */
  if (si->num_ngb_force < MAX_NUM_OF_NEIGHBOURS_BH)
    bi->ids_ngbs_force[si->num_ngb_force] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_force;
#endif
}

#endif
