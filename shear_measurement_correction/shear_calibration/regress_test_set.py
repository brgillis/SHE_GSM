""" @file regress_test_set.py

    Created 9 May 2016

    Functions to perform a linear regression on test data to get m, c,
    and their errors.

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
from scipy.stats import linregress

def get_m_and_c(shears, measured_shapes):
    
    n = np.shape(shears)[0]
    
    regression = linregress(shears, measured_shapes)

    m = regression[0]-1
    c = regression[1]

    m_err = regression[4]
    c_err = regression[4]*np.sqrt((np.sum(np.square(shears))) / n)
    
    return m, c, m_err, c_err