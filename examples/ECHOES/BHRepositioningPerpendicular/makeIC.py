###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Sagan Sutherland (sutherland.sagan@gmail.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

import h5py
import numpy as np
import unyt as u
from swiftsimio import Writer, cosmo_array, cosmo_quantity
from swiftsimio.metadata.writer.unit_systems import cosmo_units

random = np.random.default_rng()
n_p = 1024
l_box = 2
a = 1
box_size = cosmo_array(
    [l_box, l_box, l_box],
    u.Mpc,
    comoving=True,
    scale_factor=a,
    scale_exponent=1,
)

w = Writer(unit_system=cosmo_units, boxsize=box_size, scale_factor=a)

halo_center = box_size * np.array([0.5, 0.5, 0.5])
halo_radius = cosmo_quantity(
    l_box * 0.1,
    u.Mpc,
    comoving=True,
    scale_factor=w.scale_factor,
    scale_exponent=1,
)

w.gas.coordinates = halo_radius * random.normal(size=(n_p, 3)) + halo_center
w.gas.coordinates %= box_size
# w.dark_matter.coordinates = random.normal(size=(n_p, 3)) * halo_radius + halo_center

initial_v = cosmo_quantity(
    10.0,
    u.km / u.s,
    comoving=True,
    scale_factor=w.scale_factor,
    scale_exponent=1,
)
v_random = initial_v / 10

w.gas.velocities = cosmo_array(
    random.uniform(-1, 1, (n_p, 3)) * v_random + initial_v * [0, 0, 1]
)
# w.dark_matter.velocities = cosmo_array(random.uniform(-1, 1, (n_p, 3)) * v_random + initial_v)

w.gas.masses = cosmo_array(
    np.full(n_p, 1e6, dtype=float),
    u.solMass,
    comoving=True,
    scale_factor=w.scale_factor,
    scale_exponent=0,
)

w.gas.internal_energy = cosmo_array(
    np.full(n_p, 1e4 / 1e6, dtype=float),
    u.kb * u.K / u.solMass,
    comoving=True,
    scale_factor=w.scale_factor,
    scale_exponent=-2,
)

w.gas.generate_smoothing_lengths()

w.black_holes.coordinates = cosmo_array([box_size * np.array([0.5, 0.5, 0.5])])
w.black_holes.velocities = cosmo_array([initial_v * np.array([0.0, 1.0, 0.0])])
w.black_holes.masses = cosmo_array(
    np.full(1, 1e6, dtype=float),
    u.solMass,
    comoving=True,
    scale_factor=w.scale_factor,
    scale_exponent=0,
)
w.write("perpendicular.hdf5")

with h5py.File("perpendicular.hdf5", "a") as f:
    f["PartType0/SmoothingLength"] = f["PartType0/SmoothingLengths"]
    f["PartType5/SmoothingLength"] = [np.mean(w.gas.smoothing_lengths)]
