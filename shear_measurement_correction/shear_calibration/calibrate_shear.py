""" @file calibrate_shear.py

    Created 13 Aug 2015

    Functions to calibrate a single shear value and get the error on the
    calibrated value.

    ---------------------------------------------------------------------

    Copyright (C) 2015 Bryan R. Gillis

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

def calibrate_shear(g, m=0., c=0., sigma_m=0., sigma_c=0.,
                    second_order = True):
    
    gp = ( g - c ) * ( 1. - m + m**2 ) # First order correction
    if not second_order:
        return gp
    else:
        gpp = gp * ( 1. - sigma_m**2*(1+m) - m**3 ) # Second-order correction
        return gpp

def get_error_of_calibrated_shear(sigma_g, m=0., c=0., sigma_m=0., sigma_c=0.):
    
    return sigma_g*(1. + m + 1.5*m**2 + 1.5*sigma_m**2)