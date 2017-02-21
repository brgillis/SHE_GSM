#!/usr/bin/env python

""" @file run_calibration_testing.py

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

import sys
import cPickle

from shear_calibration.mock_calibration import perform_multiple_mock_calibrations

bayesian = True

def main(argv):
    """ @TODO main docstring
    """
    
    # test_ms = [-0.2, -0.1, 0., 0.1, 0.2]
    test_ms = [-0.2, 0.2]
    test_cs = [-0.1, 0., 0.1]
    
    ncal = int(1e6)
    n = int(1e4)
    
    for m in test_ms:
        for c in test_cs:
            
            print("Testing for m="+str(m)+", c="+str(c)+"...")
            
            results = perform_multiple_mock_calibrations(ncal=ncal, n=n, m=m, c=c, seed=1, nproc=-1,
                                                         second_order=False,
                                                         bayesian=bayesian)
            
            if bayesian:
                filename = "calibration_results_bayesian_m_" + str(m) + "_c_" + str(c) + ".bin"
            else:
                filename = "calibration_results_m_" + str(m) + "_c_" + str(c) + ".bin"
            
            with open(filename,'wb') as fi:
                cPickle.dump(results,fi)
                
            print("Done!")
    
    pass

if __name__ == "__main__":
    main(sys.argv)
