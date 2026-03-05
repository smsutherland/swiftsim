/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_ECHOES_BLACK_HOLES_CSDS_H
#define SWIFT_ECHOES_BLACK_HOLES_CSDS_H

/* Other Includes */
#include "csds_io.h"
#include "stars.h"

#ifdef WITH_CSDS

/**
 * @brief Defines the fields to write in the CSDS.
 *
 * @param fields (output) The list of fields to write (already allocated).
 *
 * @return The number of fields.
 */
INLINE static int csds_black_holes_define_fields(struct csds_field *fields) {

  int num_fields = 0;

  /* Positions */
  csds_define_standard_field(fields[num_fields++], "Coordinates", struct bpart,
                             x);

  /* Velocities */
  csds_define_standard_field(fields[num_fields++], "Velocities", struct bpart,
                             v);

  /* Masses */
  csds_define_standard_field(fields[num_fields++], "Masses", struct bpart,
                             mass);

  /* Particle IDs */
  csds_define_standard_field(fields[num_fields++], "ParticleIDs", struct bpart,
                             id);

  return num_fields;
}

#endif  // WITH_CSDS
#endif  // SWIFT_ECHOES_BLACK_HOLES_CSDS_H
