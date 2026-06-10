/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_ECHOES_BLACK_HOLE_PART_H
#define SWIFT_ECHOES_BLACK_HOLE_PART_H

#include "chemistry_struct.h"
#include "particle_splitting_struct.h"
#include "timeline.h"

/**
 * @brief Particle fields for the black hole particles.
 *
 * All quantities related to gravity are stored in the associate #gpart.
 */
struct bpart {

  /*! Particle ID. */
  long long id;

  /*! Pointer to corresponding gravity part. */
  struct gpart *gpart;

  /*! Particle position. */
  double x[3];

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Particle velocity. */
  float v[3];

  /*! Black hole mass */
  float mass;

  /* Particle cutoff radius. */
  float h;

  /*! Particle time bin */
  timebin_t time_bin;

  /*! Tree-depth at which size / 2 <= h * gamma < size */
  char depth_h;

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density;

  struct {

    /*! Gravitational potential copied from the #gpart. */
    float potential;

    /*! Value of the minimum potential across all neighbours. */
    float min_potential;

    /*! Delta position to apply after the reposition procedure */
    double delta_x[3];

  } reposition;

  /*! Splitting structure */
  struct particle_splitting_data split_data;

  /*! Chemistry information (e.g. metal content at birth, swallowed metal
   * content, etc.) */
  struct chemistry_bpart_data chemistry_data;

  /*! Black holes merger information (e.g. merging ID) */
  struct black_holes_bpart_data merger_data;

  /*! Tracer structure */
  struct tracers_bpart_data tracers_data;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  /*! Number of interactions in the density SELF and PAIR */
  int num_ngb_density;

  /*! List of interacting particles in the density SELF and PAIR */
  long long ids_ngbs_density[MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES];

  /*! Number of interactions in the force SELF and PAIR */
  int num_ngb_force;

  /*! List of interacting particles in the force SELF and PAIR */
  long long ids_ngbs_force[MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES];
#endif

  /* Number of BH mergers this particular particle has experienced. */
  int number_of_mergers;

  /* Total number of BHs which have merged to form this BH. */
  int cumulative_number_of_seeds;

  /*! Total number of times the black hole has been repositioned (excluding
   * repositionings of merged-in black holes) */
  int number_of_repositions;

  /*! Total number of times a black hole attempted repositioning (including
   * cases where it was aborted because the black hole was already at a
   * lower potential than all eligible neighbours) */
  int number_of_reposition_attempts;

  struct {

    /*! Sum of masses of all particles in the group. */
    float group_mass;

    /*! Maximum group mass ever reached by the group. */
    float max_group_mass;

    /*! Boolean flag for whether the BH particle thinks it's a central BH
     * sutherland TODO: Maybe save some spaceby reducing the bit footprint of
     * the bool? For now all the fields are at least 32 bits, so reducing this
     * to 1 bit wouldn't save any space, but if we have more flags in the
     * future, we could combine them into a single int, or use C bitfields. What
     * version of C introduced bitfields? If it's too late we might not want to
     * include them so as not to break support for older compilers. Unless we
     * don't really care? Maybe it's better to say we're using a more recent gcc
     * version. */
    int is_central;

    /* How far is the BH to the center of mass of its group?
     * sutherland TODO: How should we handle cases where a BH is no longer in a
     * group? Options (that I see): Set this to 0, set it to FLT_MAX, keep it at
     * its previous value */
    float distance_to_CoM;

    /* Virial radius of the FoF group.
     *
     * This uses the dark matter mass of the halo as a standin for the virial
     * mass. This assumption works best when the linking length ratio is 0.2.
     * The radius is calculated assuming a sphere of some multiple of the
     * present critical density. The overdensity is based on the fit from Bryan
     * & Norman 1998. */
    float virial_radius;

  } fof_properties;

  /* sutherland: note to DAA - properties that we want to store for trinity can
   * go in here. */
  struct {
  } echos_properties;
} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_ECHOES_BLACK_HOLE_PART_H */
