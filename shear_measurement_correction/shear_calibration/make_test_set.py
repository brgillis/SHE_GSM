""" @file make_test_set.py

    Created 9 May 2016

    Function to make a set of test shear data.

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

import numpy as np

from shear_calibration import magic_values as mv

def add_ell_1d(g1,g2):
    return (g1+g2)/(1+g1*g2)

def contract_rand(x, max=1.,p=4.):
    return (x / (1. + (np.abs(x) / max)**p)**(1./p));

def contracted_Gaus_rand(sigma=1.,max=1.,p=4.,n=1):
    
    x = np.random.randn(n)*sigma
    
    return contract_rand(x, max, p)

def make_test_set( n = mv.default_n,
                   
                   shear_sigma = mv.default_shear_sigma,
                   shape_sigma = mv.default_shape_sigma,
                   
                   m = 0.,
                   c = 0.,
                   
                   ell_trunc_max = mv.ell_trunc_max,
                   ell_trunc_p = mv.ell_trunc_p,
                   
                   proper_shear_addition = False,
                   
                   seed = mv.default_seed ):
    """
        @TODO: Function docstring
    """
    
    if seed != 0:
        np.random.seed(seed)
    
    intrinsic_shapes = contracted_Gaus_rand(shape_sigma, ell_trunc_max, ell_trunc_p, n)
    shears = contracted_Gaus_rand(shear_sigma, ell_trunc_max, ell_trunc_p, n)
    measured_shears = (1+m)*shears + c
    
    if proper_shear_addition:
        shapes = add_ell_1d(intrinsic_shapes, measured_shears)
    else:
        shapes = intrinsic_shapes + measured_shears
    
    return shears, shapes