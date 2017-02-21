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

from multiprocessing import cpu_count, Pool

import numpy as np

from shear_calibration import magic_values as mv
from shear_calibration.make_test_set import make_test_set
from shear_calibration.regress_test_set import get_m_and_c
from shear_calibration.calibrate_shear import calibrate_shear, calibrate_shear_bayesian,\
    get_bayesian_m_error_of_calibrated_shear,\
    get_bayesian_c_error_of_calibrated_shear,\
    get_bayesian_g_error_of_calibrated_shear

def perform_mock_calibration(n = mv.default_n,
                             
                             m = 0.,
                             c = 0.,
                             
                             shear_sigma = mv.default_shear_sigma,
                             shape_sigma = mv.default_shape_sigma,
                   
                             ell_trunc_max = mv.ell_trunc_max,
                             ell_trunc_p = mv.ell_trunc_p,
                               
                             proper_shear_addition = False,
                             
                             second_order = True,
                             bayesian = False,
                             
                             seed = mv.default_seed):
    
    # Start by a standard measurement and calibration of a mock test set with the given parameters
    shears, measured_shapes = make_test_set(n, shear_sigma, shape_sigma, m, c,
                                            ell_trunc_max, ell_trunc_p, proper_shear_addition,
                                            seed=seed)
    
    results = {}
    
    mm, cm, sigma_mm, sigma_cm = get_m_and_c(shears, measured_shapes)
    
    results['mm'] = mm
    results['cm'] = cm
    results['sigma_mm'] = sigma_mm
    results['sigma_cm'] = sigma_cm
    
    if bayesian:
        calibrated_shapes = calibrate_shear_bayesian(measured_shapes, mm, cm, shape_sigma, sigma_mm, sigma_cm)
        results['pred_sigma_mm'] = get_bayesian_m_error_of_calibrated_shear(mm,sigma_mm)
        results['pred_sigma_cm'] = get_bayesian_c_error_of_calibrated_shear(mm,sigma_mm,sigma_cm)
        results['pred_sigma_gm'] = get_bayesian_g_error_of_calibrated_shear(mm,sigma_mm,shape_sigma)
    else:
        calibrated_shapes = calibrate_shear(measured_shapes, mm, cm, sigma_mm, sigma_cm,
                                            second_order=second_order)
    
    # Get the known bias components after the correction
    results['mmp'], results['cmp'], _, _ = get_m_and_c(shears, calibrated_shapes)
    
    
    # Now perform the correction on noise-free data to see what the actual bias is after correction
    noise_free_shears = np.linspace(-5*shear_sigma,5*shear_sigma,100)
    noise_free_measured_shears = (1+m)*noise_free_shears + c
    noise_free_calibrated_shears = calibrate_shear(noise_free_measured_shears, mm, cm,
                                                   sigma_mm, sigma_cm,
                                                   second_order=second_order)
    
    # Get the true bias components after the correction
    results['mp'], results['cp'], _, _ = get_m_and_c(noise_free_shears, noise_free_calibrated_shears)
    
    # Cleanup
    del (shears, measured_shapes, calibrated_shapes, noise_free_shears, noise_free_measured_shears, 
        noise_free_calibrated_shears)
    
    return results

class perform_mock_calibration_caller(object):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
    def __call__(self, seed):
        return perform_mock_calibration(*self.args,seed=seed,**self.kwargs)
    
def perform_multiple_mock_calibrations(ncal = mv.default_ncal,
                                       n = mv.default_n,
                             
                                       m = 0.1,
                                       c = -0.3,
                                         
                                       shear_sigma = mv.default_shear_sigma,
                                       shape_sigma = mv.default_shape_sigma,
                               
                                       ell_trunc_max = mv.ell_trunc_max,
                                       ell_trunc_p = mv.ell_trunc_p,
                                           
                                       proper_shear_addition = False,
                                       
                                       second_order = True,
                                       bayesian = False,
                                         
                                       seed = mv.default_seed,
                                       nproc=-1):
    if(nproc <= 0):
        nproc += cpu_count()
    
    result_arrays = {}
        
    chunksize = int(1e4)
    
    nchunks = ncal // chunksize    
    if chunksize*nchunks < ncal:
        nchunks += 1

    for i in range(nchunks):
        
        size = min(chunksize,ncal-i*chunksize)

        if seed == 0:
            seeds = np.zeros(size, dtype=int)
        else:
            seeds = np.asarray(seed*ncal+np.linspace(i*chunksize, i*chunksize + size - 1, size),dtype=int)
        
        # If we just have one thread, we'll just use a simple function call to ease debugging
        if nproc == 1:
            all_calibration_results = []
            for local_seed in seeds:
                all_calibration_results.append( perform_mock_calibration(n, m, c,
                                                                         shear_sigma,
                                                                         shape_sigma,
                                                                         ell_trunc_max,
                                                                         ell_trunc_p,
                                                                         proper_shear_addition,
                                                                         second_order=second_order,
                                                                         bayesian=bayesian,
                                                                         seed=local_seed) )
        else:
                
            pool = Pool(processes=nproc,maxtasksperchild=10)
            all_calibration_results = pool.map(perform_mock_calibration_caller(n, m, c,
                                                                             shear_sigma,
                                                                             shape_sigma,
                                                                             ell_trunc_max,
                                                                             ell_trunc_p,
                                                                             proper_shear_addition,
                                                                             second_order=second_order,
                                                                             bayesian=bayesian),
                                               seeds,chunksize=10)
            
            pool.terminate()
            
        for calibration_results in all_calibration_results:
            # Append results to the result arrays
            for key in calibration_results:
                if key not in result_arrays:
                    result_arrays[key] = []
                result_arrays[key].append(calibration_results[key])
                
        del all_calibration_results
            
    # Get the results for this
    results = {}
    
    for key in result_arrays:
        
        result_arrays[key] = np.array(result_arrays[key])
        
        # Add the mean and sigma to the results
        results[key+"_mean"] = np.mean(result_arrays[key])
        results[key+"_sigma"] = np.std(result_arrays[key])
        results[key+"_stderr"] = results[key+"_sigma"] / np.sqrt(np.shape(result_arrays[key])[0]-1)
        
    return results