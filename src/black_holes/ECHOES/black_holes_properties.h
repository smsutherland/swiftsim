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
#ifndef SWIFT_ECHOES_BLACK_HOLES_PROPERTIES_H
#define SWIFT_ECHOES_BLACK_HOLES_PROPERTIES_H

#include "chemistry.h"
#include "fof.h"
#include "hydro_properties.h"

#include <string.h>

enum BH_merger_threshold {
  /*! BHs will merge if their relative velocity is less than their escape
     velocity. */
  BH_mergers_escape_velocity,

  /*! BHs will merge if their relative velocity is less than twice their
   * escape velocity. */
  BH_mergers_2_escape_velocity,

  /*! BHs will merge if one is with the kernel of the other. */
  BH_mergers_kernel,
};

enum BH_central_criterion {
  /*! A BH is considered central if it has the highest peak group mass of all
   * the groups BHs. */
  BH_central_peak_mass,
};

/**
 * @brief Properties of the black hole scheme.
 *
 * In this default scheme, we only have the properties
 * required by the neighbour search.
 */
struct black_holes_props {

  /*! Default group ID to give to particles not in a group.
   *  We steal this from the FoF parameters because we need to check when
   * disallowing integroup mergers. */
  size_t group_id_default;

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /*! Use nibbling? (Always set to 0 in the default model) */
  int use_nibbling;

  /*! Maximal distance over which BHs merge, in units of softening length */
  float max_merging_distance_ratio;

  /*! Maximal distance to reposition, in units of softening length */
  float max_reposition_distance_ratio;

  /*! Which criterion for black hole mergers are we using? */
  enum BH_merger_threshold merger_threshold_type;

  /*! Which criterion for black holes to be considered central are we using? */
  enum BH_central_criterion central_criterion;

  /*! Should black holes be allowed to merge when in different FoF groups? */
  int allow_intergroup_mergers;

  /*! Should we consider relative velocity when doing BH repositioning? */
  int enable_repos_v_threshold;

  /*! Maximum relative peculiar velocity squared that a particle can have for a
   * BH to reposition to it. */
  float repos_v2_threshold;
};

/**
 * @brief Initialise the black hole properties from the parameter file.
 *
 * We read the defaults from the hydro properties.
 *
 * @param bp The #black_holes_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
static INLINE void black_holes_props_init(struct black_holes_props *bp,
                                          const struct phys_const *phys_const,
                                          const struct unit_system *us,
                                          struct swift_params *params,
                                          const struct hydro_props *hydro_props,
                                          const struct cosmology *cosmo) {

  /* Kernel properties */
  bp->eta_neighbours = parser_get_opt_param_float(
      params, "BlackHoles:resolution_eta", hydro_props->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  bp->h_tolerance = parser_get_opt_param_float(params, "BlackHoles:h_tolerance",
                                               hydro_props->h_tolerance);

  /* Get derived properties */
  bp->target_neighbours = pow_dimension(bp->eta_neighbours) * kernel_norm;
  const float delta_eta = bp->eta_neighbours * (1.f + bp->h_tolerance);
  bp->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(bp->eta_neighbours)) *
      kernel_norm;

  /* Number of iterations to converge h */
  bp->max_smoothing_iterations =
      parser_get_opt_param_int(params, "BlackHoles:max_ghost_iterations",
                               hydro_props->max_smoothing_iterations);

  /* Time integration properties */
  const float max_volume_change =
      parser_get_opt_param_float(params, "BlackHoles:max_volume_change", -1);
  if (max_volume_change == -1)
    bp->log_max_h_change = hydro_props->log_max_h_change;
  else
    bp->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

  /* No nibbling in this default model! */
  bp->use_nibbling = 0;

  char temp[40];
  parser_get_param_string(params, "ECHOES:merger_threshold_type", temp);
  if (!strcmp(temp, "EscapeVelocity")) {
    bp->merger_threshold_type = BH_mergers_escape_velocity;
  } else if (!strcmp(temp, "EscapeVelocity2")) {
    bp->merger_threshold_type = BH_mergers_2_escape_velocity;
  } else if (!strcmp(temp, "Kernel")) {
    bp->merger_threshold_type = BH_mergers_kernel;
  } else {
    error(
        "The galaxy merger model must be one of EscapeVelocity, "
        "EscapeVelocity2, or Kernel, "
        "not %s",
        temp);
  }

  bp->max_merging_distance_ratio =
      parser_get_param_float(params, "ECHOES:merger_max_distance_ratio");

  bp->max_reposition_distance_ratio =
      parser_get_param_float(params, "ECHOES:max_reposition_distance_ratio");

  parser_get_param_string(params, "ECHOES:central_galaxy_criterion", temp);
  if (!strcmp(temp, "PeakMass")) {
    bp->central_criterion = BH_central_peak_mass;
  } else {
    error("The galaxy central criterion must be one of PeakMass, not %s", temp);
  }

  bp->allow_intergroup_mergers =
      parser_get_param_int(params, "ECHOES:allow_intergroup_mergers");

  /* We steal this from the FoF props */
  bp->group_id_default = parser_get_opt_param_int(
      params, "FOF:group_id_default", fof_props_default_group_id);

  bp->enable_repos_v_threshold = parser_get_param_int(
      params, "ECHOES:enable_reposition_velocity_threshold");
  if (bp->enable_repos_v_threshold) {
    float repos_v_threshold =
        parser_get_param_float(params,
                               "ECHOES:reposition_velocity_threshold_km_s") *
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));
    bp->repos_v2_threshold = repos_v_threshold * repos_v_threshold;
  } else {
    bp->repos_v2_threshold = FLT_MAX;
  }
}

/**
 * @brief Write a black_holes_props struct to the given FILE as a stream of
 * bytes.
 *
 * @param props the black hole properties struct
 * @param stream the file stream
 */
INLINE static void black_holes_struct_dump(
    const struct black_holes_props *props, FILE *stream) {
  restart_write_blocks((void *)props, sizeof(struct black_holes_props), 1,
                       stream, "black_holes props", "black holes props");
}

/**
 * @brief Restore a black_holes_props struct from the given FILE as a stream
 * of bytes.
 *
 * @param props the black hole properties struct
 * @param stream the file stream
 */
INLINE static void black_holes_struct_restore(
    const struct black_holes_props *props, FILE *stream) {
  restart_read_blocks((void *)props, sizeof(struct black_holes_props), 1,
                      stream, NULL, "black holes props");
}

#endif /* SWIFT_ECHOES_BLACK_HOLES_PROPERTIES_H */
