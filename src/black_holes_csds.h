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
#ifndef SWIFT_BLACK_HOLES_CSDS_H
#define SWIFT_BLACK_HOLES_CSDS_H

/* Include config */
#include <config.h>

/* Local includes */
#include "./const.h"
#include "align.h"
#include "csds.h"
#include "part_type.h"
#include "timeline.h"

/* Select the correct BH model */
#if defined(BLACK_HOLES_NONE)
#error TODO
#elif defined(BLACK_HOLES_EAGLE)
#error TODO
#elif defined(BLACK_HOLES_SPIN_JET)
#error TODO
#elif defined(BLACK_HOLES_ECHOES)
#include "./black_holes/ECHOES/black_holes_csds.h"
#else
#error "Invalid choice of black hole model"
#endif

#endif /* SWIFT_BLACK_HOLES_CSDS_H */
