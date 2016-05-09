""" @file magic_values.py

    Created 9 May 2016

    Magic values for the shear_measurement_correction project

    ---------------------------------------------------------------------

    Copyright (C) 2016 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

default_ncal = int(1e4)
default_n = int(1e4)

default_seed = 0

# Sigma for the shear and intrinsic ellipticity ("shape") distributions
default_shear_sigma = 0.03
default_shape_sigma = 0.25

# Parameters to handle truncation of ellipticity measurements
ell_trunc_max = 0.9
ell_trunc_p = 4

