""" @file mock_calibration.py

    Created 9 May 2016

    @TODO: File docstring

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
from shear_calibration.make_test_set import make_test_set
from shear_calibration.regress_test_set import get_m_and_c
from shear_calibration.calibrate_shear import calibrate_shear

def perform_mock_calibration(n = mv.default_n,
                             
                             m = 0.,
                             c = 0.,
                             
                             shear_sigma = mv.default_shear_sigma,
                             shape_sigma = mv.default_shape_sigma,
                   
                             ell_trunc_max = mv.ell_trunc_max,
                             ell_trunc_p = mv.ell_trunc_p,
                               
                             proper_shear_addition = False,
                             
                             seed = mv.default_seed):
    
    # Start by a standard measurement and calibration of a mock test set with the given parameters
    shears, measured_shapes = make_test_set(n, shear_sigma, shape_sigma, m, c,
                                            ell_trunc_max, ell_trunc_p, proper_shear_addition)
    
    mm, cm, sigma_mm, sigma_cm = get_m_and_c(shears, measured_shapes)
    
    calibrated_shapes = calibrate_shear(measured_shapes, mm, cm, sigma_mm, sigma_cm)
    
    # Get the known bias components after the correction
    mmp, cmp, sigma_mmp, sigma_cmp = get_m_and_c(shears, measured_shapes)
    
    
    # Now perform the correction on noise-free data to see what the actual bias is after correction
    noise_free_shears = np.linspace(-5*shear_sigma,5*shear_sigma,100)
    noise_free_measured_shears = (1+m)*noise_free_shears + c
    noise_free_calibrated_shears = calibrate_shear(noise_free_measured_shears, mm, cm,
                                                   sigma_mm, sigma_cm)
    
    # Get the true bias components after the correction
    mp, cp, _sigma_mp, _sigma_cp = get_m_and_c(noise_free_shears, noise_free_calibrated_shears)
    
    
    return ( mm, cm,
             sigma_mm, sigma_cm,
             mp, cp )
    