#!/usr/bin/env python

""" @file /disk2/brg/git/SHE_sim/shear_measurement_correction/plot_calibration_change.py

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
    
markersize = 8
fontsize = 24
file_format = "eps"
paper_location = "/home/brg/Dropbox/gillis-comp-shared/Papers/Shear_Bias/"
calibration_set_size = "1e6"
calibration_order = "1st"

def main(argv):
    """ @TODO main docstring
    """
    
    if len(argv) > 1:
        calibration_set_size = argv[1]
    if len(argv) > 2:
        calibration_order = argv[2]
    
    results_dir = "calibration_results_" + calibration_set_size + "_" + calibration_order
    
    # Load in the calculated results
    test_ms = [-0.1, 0., 0.1]
    test_cs = [-0.1, 0., 0.1]
    
    results = {}
    
    for m in test_ms:
        for c in test_cs:
            
            filename = join(results_dir,"calibration_results_m_" + str(m) + "_c_" + str(c) + ".bin")
            
            with open(filename,'rb') as fi:
                new_res = cPickle.load(fi)
                new_res["c_mean"] = c
                new_res["m_mean"] = m
                new_res["c_sigma"] = 0
                new_res["m_sigma"] = 0
                results[(c,m)] = new_res
                
    # We'll plot both the change in measured values and actual values
    for m_tag, measured_label in (("m", "Measured"),
                                  ("", "Actual")):
        
        # Setup the figure
        
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(measured_label + " c",fontsize=fontsize)
        ax.set_ylabel(measured_label + " m",fontsize=fontsize)
        
        ax.set_xlim(-0.25,0.25)
        ax.set_ylim(-0.31,0.31)
        
        ax.set_yscale("symlog",linthreshy=0.01)
        ax.set_xscale("symlog",linthreshx=0.0001)
        
        # Draw axes
        ax.plot([-1,1],[0,0],label=None,color="k",linestyle="solid",linewidth=0.5)
        ax.plot([0,0],[-1,1],label=None,color="k",linestyle="solid",linewidth=0.5)
        
        ax.plot([-1,1],[0.1,0.1],label=None,color="k",linestyle="dotted",linewidth=0.5)
        ax.plot([-1,1],[-0.1,-0.1],label=None,color="k",linestyle="dotted",linewidth=0.5)
        ax.plot([0.1,0.1],[-1,1],label=None,color="k",linestyle="dotted",linewidth=0.5)
        ax.plot([-0.1,-0.1],[-1,1],label=None,color="k",linestyle="dotted",linewidth=0.5)
    
        # Plot the points
        
        for m in test_ms:
            r_color = (m+0.1)/0.2
            for c in test_cs:
                b_color = (c+0.1)/0.2
                
                # Skip the central point
                if m==0 and c==0:
                    continue
                
                color = (r_color,1-(r_color+b_color)/2,b_color)
                
                res = results[(c,m)]
                
                ax.errorbar([res["c"+m_tag+"_mean"]],[res["m"+m_tag+"_mean"]],
                            marker='.',markersize=markersize,label=str((c,m)),
                            linestyle="None",linewidth=2,color=color,
                            xerr=[res["c"+m_tag+"_sigma"]],
                            yerr=[res["m"+m_tag+"_sigma"]])
                
                eb1 = ax.errorbar([res["c"+m_tag+"p_mean"]],[res["m"+m_tag+"p_mean"]],
                            marker='o',markerfacecolor='None',markeredgecolor=color,
                            markersize=markersize,label=str((c,m))+"'",
                            linestyle="None",linewidth=2,color=color,
                            xerr=[res["c"+m_tag+"p_sigma"]],
                            yerr=[res["m"+m_tag+"p_sigma"]])
                
                eb1[-1][0].set_linestyle('dotted')
                eb1[-1][-1].set_linestyle('dotted')
    
        # Save the plot
        figname = (measured_label.lower() + "_" + calibration_set_size + "_"
                   + calibration_order + "_calibration." + file_format)
        if file_format == "eps":
            # Save in the paper location
            figname = join(paper_location,figname)
        pyplot.savefig(figname, format=file_format, bbox_inches="tight", pad_inches=0.05)
        
        fig.show()
        
    pyplot.show()
    
    return

if __name__ == "__main__":
    main(sys.argv)
