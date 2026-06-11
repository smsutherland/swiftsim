/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_ECHOES_BLACK_HOLES_H
#define SWIFT_ECHOES_BLACK_HOLES_H

#include <float.h>

/* Local includes */
#include "black_holes_properties.h"
#include "black_holes_struct.h"
#include "cooling_properties.h"
#include "dimension.h"
#include "error.h"
#include "fof.h"
#include "gravity.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "random.h"

/**
 * @brief Computes the time-step of a given black hole particle.
 *
 * @param bp Pointer to the s-particle data.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 */
__attribute__((always_inline)) INLINE static float black_holes_compute_timestep(
    const struct bpart *const bp, const struct black_holes_props *props,
    const struct phys_const *constants, const struct cosmology *cosmo) {

  return FLT_MAX;
}

/**
 * @brief Initialises the b-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param bp The particle to act upon
 * @param props The properties of the black holes model.
 */
__attribute__((always_inline)) INLINE static void black_holes_first_init_bpart(
    struct bpart *bp, const struct black_holes_props *props) {

  bp->time_bin = 0;

  bp->fof_properties.group_mass = 0.f;
  bp->fof_properties.max_group_mass = 0.f;
  bp->fof_properties.distance_to_CoM = 0.f;
  bp->fof_properties.is_central = 0;
}

/**
 * @brief Prepares a b-particle for its interactions
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_init_bpart(
    struct bpart *bp) {

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_density[i] = -1;
  bp->num_ngb_density = 0;
#endif

  bp->density.wcount = 0.f;
  bp->density.wcount_dh = 0.f;
  bp->reposition.delta_x[0] = -FLT_MAX;
  bp->reposition.delta_x[1] = -FLT_MAX;
  bp->reposition.delta_x[2] = -FLT_MAX;
  bp->reposition.min_potential = FLT_MAX;
  bp->reposition.potential = FLT_MAX;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * The fields do not get predicted but we move the BH to its new position
 * if a new one was calculated in the repositioning loop.
 *
 * @param bp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void black_holes_predict_extra(
    struct bpart *restrict bp, float dt_drift) {
  /* sutherland: This seems to be where the repositioning happens. */

  /* Quit early if we're not repositioning */
  if (bp->reposition.min_potential == FLT_MAX) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (bp->reposition.delta_x[0] == -FLT_MAX ||
      bp->reposition.delta_x[1] == -FLT_MAX ||
      bp->reposition.delta_x[2] == -FLT_MAX) {
    error("Something went wrong with the new repositioning position");
  }

  const double dx = bp->reposition.delta_x[0];
  const double dy = bp->reposition.delta_x[1];
  const double dz = bp->reposition.delta_x[2];
  const double d = sqrt(dx * dx + dy * dy + dz * dz);
  if (d > 1.01 * kernel_gamma * bp->h)
    error("Repositioning BH beyond the kernel support!");
#endif

  /* Move the black hole */
  bp->x[0] += bp->reposition.delta_x[0];
  bp->x[1] += bp->reposition.delta_x[1];
  bp->x[2] += bp->reposition.delta_x[2];

  /* Move its gravity properties as well */
  bp->gpart->x[0] += bp->reposition.delta_x[0];
  bp->gpart->x[1] += bp->reposition.delta_x[1];
  bp->gpart->x[2] += bp->reposition.delta_x[2];

  /* Store the delta position */
  bp->x_diff[0] -= bp->reposition.delta_x[0];
  bp->x_diff[1] -= bp->reposition.delta_x[1];
  bp->x_diff[2] -= bp->reposition.delta_x[2];

  /* Reset the reposition variables */
  bp->reposition.delta_x[0] = -FLT_MAX;
  bp->reposition.delta_x[1] = -FLT_MAX;
  bp->reposition.delta_x[2] = -FLT_MAX;
  bp->reposition.min_potential = FLT_MAX;

  /* Count the jump */
  bp->number_of_repositions++;
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param bp The particle.
 */
__attribute__((always_inline)) INLINE static void
black_holes_reset_predicted_values(struct bpart *restrict bp) {}

/**
 * @brief Kick the additional variables
 *
 * @param bp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void black_holes_kick_extra(
    struct bpart *bp, float dt) {}

/**
 * @brief Finishes the calculation of density on black holes
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_end_density(
    struct bpart *bp, const struct cosmology *cosmo) {

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  bp->density.wcount *= h_inv_dim;
  bp->density.wcount_dh *= h_inv_dim_plus_one;
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
black_holes_bpart_has_no_neighbours(struct bpart *restrict bp,
                                    const struct cosmology *cosmo) {

  warning(
      "BH particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      bp->id, bp->h, bp->density.wcount);

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  bp->density.wcount = kernel_root * h_inv_dim;
  bp->density.wcount_dh = 0.f;
}

/**
 * @brief Return the current instantaneous accretion rate of the BH.
 *
 * Empty BH model --> return 0.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_accretion_rate(const struct bpart *bp) {
  return 0.;
}

/**
 * @brief Return the total accreted gas mass of this BH.
 *
 * Empty BH model --> return 0.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_accreted_mass(const struct bpart *bp) {
  return 0.;
}

/**
 * @brief Return the subgrid mass of this BH.
 *
 * Empty BH model --> return 0.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_subgrid_mass(const struct bpart *bp) {
  return 0.;
}

/**
 * @brief Return the current bolometric luminosity of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_bolometric_luminosity(const struct bpart *bp,
                                      const struct phys_const *constants) {
  return 0.;
}

/**
 * @brief Return the current kinetic jet power of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double black_holes_get_jet_power(
    const struct bpart *bp, const struct phys_const *constants) {
  return 0.;
}

/**
 * @brief Update the properties of a black hole particles by swallowing
 * a gas particle.
 *
 * @param bp The #bpart to update.
 * @param p The #part that is swallowed.
 * @param xp The #xpart that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_swallow_part(
    struct bpart *bp, const struct part *p, const struct xpart *xp,
    const struct cosmology *cosmo) {

  /* Nothing to do here: No swallowing in the default model */
}

/**
 * @brief Update the properties of a black hole particles by swallowing
 * a BH particle.
 *
 * @param bpi The #bpart to update.
 * @param bpj The #bpart that is swallowed.
 * @param cosmo The current cosmological model.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants in internal units.
 */
__attribute__((always_inline)) INLINE static void black_holes_swallow_bpart(
    struct bpart *bpi, const struct bpart *bpj, const struct cosmology *cosmo,
    const double time, const int with_cosmology,
    const struct black_holes_props *props, const struct phys_const *constants) {

  float bpi_mass = bpi->mass;
  float bpj_mass = bpj->mass;

  /* Update mass */
  bpi->mass += bpj_mass;
  bpi->gpart->mass += bpj_mass;

  /* Conservation of momentum */
  const float BH_mom[3] = {bpi_mass * bpi->v[0] + bpj_mass * bpj->v[0],
                           bpi_mass * bpi->v[1] + bpj_mass * bpj->v[1],
                           bpi_mass * bpi->v[2] + bpj_mass * bpj->v[2]};

  bpi->v[0] = BH_mom[0] / bpi->mass;
  bpi->v[1] = BH_mom[1] / bpi->mass;
  bpi->v[2] = BH_mom[2] / bpi->mass;
  bpi->gpart->v_full[0] = bpi->v[0];
  bpi->gpart->v_full[1] = bpi->v[1];
  bpi->gpart->v_full[2] = bpi->v[2];

  bpi->number_of_mergers++;
  bpi->cumulative_number_of_seeds += bpj->cumulative_number_of_seeds;

  bpi->fof_properties.is_central |= bpj->fof_properties.is_central;
  bpi->fof_properties.max_group_mass = fmaxf(
      bpi->fof_properties.max_group_mass, bpj->fof_properties.max_group_mass);
  bpi->fof_properties.group_mass =
      fmaxf(bpi->fof_properties.group_mass, bpj->fof_properties.group_mass);
}

/**
 * @brief Compute the accretion rate of the black hole and all the quantites
 * required for the feedback loop.
 *
 * Nothing to do here.
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param cosmo The cosmological model.
 * @param cooling Properties of the cooling model.
 * @param floor_props Properties of the entropy fllor.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param dt The time-step size (in physical internal units).
 */
__attribute__((always_inline)) INLINE static void black_holes_prepare_feedback(
    struct bpart *restrict bp, const struct black_holes_props *props,
    const struct phys_const *constants, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling,
    const struct entropy_floor_properties *floor_props, const double time,
    const int with_cosmology, const double dt, const integertime_t ti_begin) {
  /* sutherland: note to DAA - This is a good place to update all of our
   * echoes_properties. If we want to do stuff with virialization, the cosmology
   * struct has a field overdensity_BN98 which may help. */
}

/**
 * @brief Finish the calculation of the new BH position.
 *
 * Here, we check that the BH should indeed be moved in the next drift.
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param cosmo The cosmological model.
 * @param dt The black hole particle's time step.
 * @param ti_begin The time at the start of the temp
 */
__attribute__((always_inline)) INLINE static void black_holes_end_reposition(
    struct bpart *restrict bp, const struct black_holes_props *props,
    const struct phys_const *constants, const struct cosmology *cosmo,
    const double dt, const integertime_t ti_begin) {

  /* sutherland TODO: EAGLE has a maximum mass at which repositions happen */

  /* First check: did we find any eligible neighbour particle to jump to? */
  if (bp->reposition.min_potential == FLT_MAX) return;

  /* Record that we have a (possible) repositioning situation */
  bp->number_of_reposition_attempts++;

  /* Is the potential lower (i.e. the BH is at the bottom already) */
  const float potential = gravity_get_comoving_potential(bp->gpart);
  if (potential < bp->reposition.min_potential) {

    /* No need to reposition */
    bp->reposition.min_potential = FLT_MAX;
    bp->reposition.delta_x[0] = -FLT_MAX;
    bp->reposition.delta_x[1] = -FLT_MAX;
    bp->reposition.delta_x[2] = -FLT_MAX;

  } else {
    /* We _should_ reposition, but not fractionally. Here, we will
     * reposition exactly on top of another gas particle - which
     * could cause issues, so we add on a small fractional offset
     * of magnitude 0.001 h in the reposition delta. */

    /* Generate three random numbers in the interval [-0.5, 0.5]; id,
     * id**2, and id**3 are required to give unique random numbers (as
     * random_unit_interval is completely reproducible). */
    const float offset_dx =
        random_unit_interval(bp->id, ti_begin, random_number_BH_reposition) -
        0.5f;
    const float offset_dy = random_unit_interval(bp->id * bp->id, ti_begin,
                                                 random_number_BH_reposition) -
                            0.5f;
    const float offset_dz =
        random_unit_interval(bp->id * bp->id * bp->id, ti_begin,
                             random_number_BH_reposition) -
        0.5f;

    const float length_inv =
        1.0f / sqrtf(offset_dx * offset_dx + offset_dy * offset_dy +
                     offset_dz * offset_dz);

    const float norm = 0.001f * bp->h * length_inv;

    bp->reposition.delta_x[0] += offset_dx * norm;
    bp->reposition.delta_x[1] += offset_dy * norm;
    bp->reposition.delta_x[2] += offset_dz * norm;
  }
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on black hole, therefore no need to use
 * it.
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_reset_feedback(
    struct bpart *restrict bp) {

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_force[i] = -1;
  bp->num_ngb_force = 0;
#endif
}

/**
 * @brief Store the gravitational potential of a black hole by copying it from
 * its #gpart friend.
 *
 * @param bp The black hole particle.
 * @param gp The black hole's #gpart.
 */
__attribute__((always_inline)) INLINE static void
black_holes_store_potential_in_bpart(struct bpart *bp, const struct gpart *gp) {

#ifdef SWIFT_DEBUG_CHECKS
  if (bp->gpart != gp) error("Copying potential to the wrong black hole!");
#endif

  bp->reposition.potential = gp->potential;
}

/**
 * @brief Store the gravitational potential of a particle by copying it from
 * its #gpart friend.
 *
 * @param p_data The black hole data of a gas particle.
 * @param gp The black hole's #gpart.
 */
__attribute__((always_inline)) INLINE static void
black_holes_store_potential_in_part(struct black_holes_part_data *p_data,
                                    const struct gpart *gp) {
  p_data->potential = gp->potential;
}

/**
 * @brief Initialise a BH particle that has just been seeded.
 *
 * @param bp The #bpart to initialise.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param p The #part that became a black hole.
 * @param xp The #xpart that became a black hole.
 * @param ti_current the current time on the time-line.
 */
INLINE static void black_holes_create_from_gas(
    struct bpart *bp, const struct black_holes_props *props,
    const struct phys_const *constants, const struct cosmology *cosmo,
    const struct part *p, const struct xpart *xp,
    const integertime_t ti_current) {

  /* The BH itself is its only seed. */
  bp->cumulative_number_of_seeds = 1;

  /* It's just a baby! No mergers yet. */
  bp->number_of_mergers = 0;

  /* Likewise it's not been swallowed yet either */
  black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);

  /* First initialisation */
  black_holes_init_bpart(bp);

  bp->fof_properties.group_mass = 0.f;
  bp->fof_properties.max_group_mass = 0.f;
  bp->fof_properties.distance_to_CoM = 0.f;
  bp->fof_properties.is_central = 0;
  bp->fof_properties.virial_radius = 0.;
}

/**
 * @brief Store FoF-related properties in a #bpart.
 *
 * @param props The properties of the BH scheme.
 * @param r2 Comoving square distance between the BH and the FoF center of mass.
 * @param group_mass Mass of the FoF group the BH is in. 0. if the BH is not in
 * a group.
 * @param is_central Is the BH central? BHs outside of FoF groups are never
 * central.
 * @param bp The black hole to update.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
black_holes_update_fof_properties(const struct black_holes_props *const props,
                                  float r2, float group_mass, int is_central,
                                  struct bpart *const bp,
                                  const struct cosmology *const cosmo) {
  if (bp->gpart->fof_data.group_id == props->group_id_default) {
    /* BHs that are not in a group have their group data reset */
    bp->fof_properties.group_mass = 0.f;
    bp->fof_properties.distance_to_CoM = 0.f;
    bp->fof_properties.is_central = 0;
  } else {

    /* sutherland TODO: See if we can handle this such
     * that we don't actually need to store r. If storing r^2 is ok, then we
     * could convert to r only when outputting for a snapshot. */
    bp->fof_properties.distance_to_CoM = sqrtf(r2);

    bp->fof_properties.group_mass = group_mass;

    bp->fof_properties.is_central = is_central;
    /* If we are the central BH, then we update the max group data. */
    if (is_central) {
      bp->fof_properties.max_group_mass =
          fmaxf(bp->fof_properties.max_group_mass, group_mass);

      float virial_density_phys =
          cosmo->critical_density * cosmo->overdensity_BN98;

      /* sutherland NOTE: Is there a faster way to calculate this?
       * We could potentially make use of the custom icbrtf when its enabled. */
      float virial_radius_phys =
          cbrtf(0.75 * bp->fof_properties.max_group_mass /
                (virial_density_phys * M_PI));

      /* sutherland NOTE: How often do we want to update this? The target
       * overdensity changes with redshift, so this can fall out of date.
       * black_holes_init_bpart seems like a good place to re-calculate, but
       * that doesn't have access to the cosmology. Should we just add the
       * cosmology there as a parameter?
       * Alternately, we could avoid storing this and just calculate it on the
       * fly whenever we need it. Doing this, we can even still output it in
       * snapshots. */
      bp->fof_properties.virial_radius = virial_radius_phys * cosmo->a_inv;

      if ((bp->fof_properties.virial_radius *
           props->max_merging_distance_ratio) > bp->h)
        warning(
            "BH %lld has a merger distance greater than it's kernel size. (%g "
            "* %g > %g)",
            bp->id, bp->fof_properties.virial_radius,
            props->max_merging_distance_ratio, bp->h);
    }
  }
}

/**
 * @brief Give this #bpart's priority for being considered the central black
 * hole in a FoF group.
 *
 * @param props The properties of the BH scheme.
 * @param bp The black hole give the priority of.
 */
__attribute__((always_inline)) INLINE static float black_holes_central_priority(
    const struct black_holes_props *const props, const struct bpart *const bp) {
  if (props->central_criterion == BH_central_peak_mass) {
    return bp->fof_properties.max_group_mass;
  } else {
#ifdef SWIFT_DEBUG_CHECKS
    error("Invalid choice of central galaxy criterion");
#else
    /* Got to return something if we don't error */
    return -1;
#endif
  }
}

#endif /* SWIFT_ECHOES_BLACK_HOLES_H */
