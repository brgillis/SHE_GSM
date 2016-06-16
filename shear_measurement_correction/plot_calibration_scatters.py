#!/usr/bin/env python

""" @file /disk2/brg/git/SHE_sim/shear_measurement_correction/plot_calibration_scatters.py

    Created 13 Jun 2016

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
from os.path import join

import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
    
markersize = 100
fontsize = 24
file_format = "eps"
paper_location = "/home/brg/Dropbox/gillis-comp-shared/Papers/Shear_Bias/"

def main(argv):
    """ @TODO main docstring
    """
    
    calibration_set_size = "1e4"
    
    if len(argv) > 1:
        calibration_set_size = argv[1]
        
    if calibration_set_size == "1e6":
        size_label = r"$N = 10^6$"
    else:
        size_label = r"$N = 10^4$"
    
    results_dir_1 = "calibration_results_" + calibration_set_size + "_1st"
    results_dir_2 = "calibration_results_" + calibration_set_size + "_2nd"
    
    # Load in the calculated results
    test_ms = [-0.1, 0., 0.1]
    test_cs = [-0.1, 0., 0.1]
    
    results = {}
    
    for m in test_ms:
        for c in test_cs:
            
            filename_1 = join(results_dir_1,"calibration_results_m_" + str(m) + "_c_" + str(c) + ".bin")
            
            with open(filename_1,'rb') as fi:
                new_res = cPickle.load(fi)
            new_res["c_mean"] = c
            new_res["m_mean"] = m
            new_res["c_sigma"] = 0
            new_res["m_sigma"] = 0
            
            filename_2 = join(results_dir_2,"calibration_results_m_" + str(m) + "_c_" + str(c) + ".bin")
            
            with open(filename_2,'rb') as fi:
                new_res_2 = cPickle.load(fi)
            new_res["cpp_sigma"] = new_res_2["cp_sigma"]
            new_res["mpp_sigma"] = new_res_2["mp_sigma"]
            
            results[(c,m)] = new_res
                
    # We'll plot both the change in measured values and actual values
        
    # Setup the figure
    
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r"$\sigma[c]$",fontsize=fontsize)
    ax.set_ylabel(r"$\sigma[m]$",fontsize=fontsize)
    
    if calibration_set_size=="1e4":
        lim_factor = 1
    else:
        lim_factor = 0.1
    ax.set_xlim(lim_factor*0.002,lim_factor*0.003)
    ax.set_ylim(lim_factor*0.07,lim_factor*0.095)
    
    xticks = [lim_factor*0.002,lim_factor*0.0025,lim_factor*0.003,]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks)
    
    yticks = [lim_factor*0.07,lim_factor*0.08,lim_factor*0.09,]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    
    ax.text(0.05,0.95,size_label,horizontalalignment='left',
            verticalalignment='top',transform=ax.transAxes,
            fontsize=24)

    # Plot the points
    
    for m in test_ms:
        r_color = (m+0.1)/0.2
        for c in test_cs:
            b_color = (c+0.1)/0.2
            
            # Only label the central point, and give it the color black instead
            if m==0 and c==0:
                label_0 = "No correction"
                label_1 = "1st-order correction"
                label_2 = "2nd-order correction"
                color='k'
            else:
                label_0 = None
                label_1 = None
                label_2 = None
                color = (r_color,1-(r_color+b_color)/2,b_color)
            
            
            res = results[(c,m)]
            
            ax.scatter([res["cm_sigma"]],[res["mm_sigma"]],
                        marker='.',s=markersize,label=label_0,color=color)
            
            ax.scatter([res["cp_sigma"]],[res["mp_sigma"]],
                        marker='o',facecolor='None',edgecolor=color,
                        s=markersize,label=label_1,color=color)
            
            ax.scatter([res["cpp_sigma"]],[res["mpp_sigma"]],
                        marker='x',s=markersize,label=label_2,color=color)
            
    ax.legend(loc='lower right',scatterpoints=1)

    # Save the plot
    figname = ("both_" + calibration_set_size + "_scatters." + file_format)
    if file_format == "eps":
        # Save in the paper location
        figname = join(paper_location,figname)
    pyplot.savefig(figname, format=file_format, bbox_inches="tight", pad_inches=0.05)
    
    fig.show()
        
    pyplot.show()
    
    return

if __name__ == "__main__":
    main(sys.argv)
