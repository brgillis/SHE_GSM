#!/usr/bin/env python

""" @file /disk2/brg/git/SHE_sim/shear_measurement_correction/plot_calibration_quality.py

    Created 11 July 2016

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
import numpy as np

import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from shear_calibration import magic_values as mv
    
markersize = 100
fontsize = 24
file_format = "eps"
paper_location = "/home/brg/Dropbox/gillis-comp-shared/Papers/Shear_Bias/"

m_target = 0.002
c_target = 0.00005

heatmap_m_min = -0.15
heatmap_m_max = 0.15

heatmap_log_N_min = 6
heatmap_log_N_max = 10

heatmap_points = 100

def main(argv):
    """ @TODO main docstring
    """
    
    for calibration_set_size in ("1e4","1e6"):
    
        if len(argv) > 1:
            calibration_set_size = argv[1]
            
        if calibration_set_size == "1e6":
            size_label = r"$n = 10^6$"
        else:
            size_label = r"$n = 10^4$"
        
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
                new_res["cpp_mean"] = new_res_2["cp_mean"]
                new_res["mpp_mean"] = new_res_2["mp_mean"]
                new_res["cpp_sigma"] = new_res_2["cp_sigma"]
                new_res["mpp_sigma"] = new_res_2["mp_sigma"]
                
                results[(c,m)] = new_res
            
        # Set up the figure
        
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel("Calibration order",fontsize=fontsize)
        ax.set_ylabel("d",fontsize=fontsize)
        
    #     if calibration_set_size=="1e4":
    #         lim_factor = 1
    #     else:
    #         lim_factor = 0.1
        ax.set_xlim(-0.5,2.5)
        if calibration_set_size == "1e6":
            ylolim = 6
            yuplim = 13
        else:
            ylolim = 60
            yuplim = 130
        ymin = ylolim + 0.05*(yuplim-ylolim)
        ymax = ylolim + 0.95*(yuplim-ylolim)
        ax.set_ylim(ylolim,yuplim)
        
        
        xticks = [0,1,2]
        ax.set_xticks(xticks)
        ax.set_xticklabels(["None","1st-order","2nd-order"])
        
        ax.text(0.95,0.95,size_label,horizontalalignment='right',
                verticalalignment='top',transform=ax.transAxes,
                fontsize=24)
    
        # Plot the points
        
        for m in test_ms:
            r_color = (m+0.1)/0.2
            for c in test_cs:
                b_color = (c+0.1)/0.2
                
                # Only label the central point, and give it the color black instead
                if m==0 and c==0:
                    color='k'
                else:
                    color = (r_color,1-(r_color+b_color)/2,b_color)
                
                res = results[(c,m)]
                
                Qm_square = (res["mm_mean"]**2 + res["mm_sigma"]**2)/ m_target**2
                Qc_square = (res["cm_mean"]**2 + res["cm_sigma"]**2)/ c_target**2
                Q_0 = np.sqrt(Qm_square + Qc_square)
                
                x_offset = c*0.5
                
                if(Q_0>=ymin and Q_0<=ymax):
                    y = Q_0
                    ax.scatter([0+x_offset],[y], marker='.',s=markersize,color=color)
                elif(Q_0>ymax):
                    y = ymax
                    dy = 0.8*(yuplim-ymax)
                    ax.errorbar([0+x_offset],[y],yerr=dy,marker='.',markersize=np.sqrt(markersize),color=color,
                                lolims=[True])
                elif(Q_0<ymin):
                    y = ymin
                    dy = 0.8*(ymin-ylolim)
                    ax.errorbar([0+x_offset],[y],yerr=dy,marker='.',markersize=np.sqrt(markersize),color=color,
                                uplims=[True])
                
                Qm_1 = np.sqrt(res["mp_mean"]**2 + res["mp_sigma"]**2)
                Qc_1 = np.sqrt(res["cp_mean"]**2 + res["cp_sigma"]**2)
                Q_1 = Qm_1 / m_target + Qc_1 / c_target
                ax.scatter([1+x_offset],[Q_1], marker='o',s=markersize,color=color,
                           facecolor='None',edgecolor=color)
                
                Qm_2 = np.sqrt(res["mpp_mean"]**2 + res["mpp_sigma"]**2)
                Qc_2 = np.sqrt(res["cpp_mean"]**2 + res["cpp_sigma"]**2)
                Q_2 = Qm_2 / m_target + Qc_2 / c_target
                ax.scatter([2+x_offset],[Q_2], marker='s',s=markersize,color=color,
                           facecolor='None',edgecolor=color)
                
                # Draw arrows connecting points
                ax.annotate("",xytext=(x_offset,y),xy=(1.15*(1+x_offset)/2.15,(y+1.15*Q_1)/2.15),
                            arrowprops=dict(arrowstyle="->",color=color))
                ax.annotate("",xytext=((1+x_offset)/2,(y+Q_1)/2),xy=(1+x_offset,Q_1),
                            arrowprops=dict(arrowstyle="->",color=color))
                
                ax.annotate("",xytext=(1+x_offset,Q_1),xy=((1+x_offset+1.15*(2+x_offset))/2.15,(Q_1+1.15*Q_2)/2.15),
                            arrowprops=dict(arrowstyle="->",color=color))
                ax.annotate("",xytext=((1+x_offset+2+x_offset)/2,(Q_1+Q_2)/2),xy=(2+x_offset,Q_2),
                            arrowprops=dict(arrowstyle="->",color=color))
    
        # Save the plot
        figname = ("both_" + calibration_set_size + "_qualities." + file_format)
        if file_format == "eps":
            # Save in the paper location
            figname = join(paper_location,figname)
        pyplot.savefig(figname, format=file_format, bbox_inches="tight", pad_inches=0.05)
        
        fig.show()
    
    # Now make a heatmap of predicted quality versus m and N
    m_vals = np.linspace(heatmap_m_min, heatmap_m_max, heatmap_points, endpoint=False)
    N_vals = np.array([np.logspace(heatmap_log_N_min, heatmap_log_N_max, heatmap_points, endpoint=False)]).transpose()
    
    S2 = (mv.default_shape_sigma/mv.default_shear_sigma)**2
    
    d_vals = np.sqrt( (3/(N_vals**2) * S2**2)*(1+2*m_vals+m_vals**2) +
                      1/N_vals * S2 * (1 - 2*m_vals - 3*m_vals**2 + 5*m_vals**3 + 5*m_vals**4) + m_vals**6)
    
    # Set up the figure now
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel("m",fontsize=fontsize)
    ax.set_ylabel("n",fontsize=fontsize)
    
    ax.axis([m_vals.min(),m_vals.max(),N_vals.min(),N_vals.max()])
    
    ax.set_yscale("log",nonposy="clip")
    
    pyplot.pcolormesh(m_vals,N_vals,d_vals,norm=matplotlib.colors.LogNorm(),figure=fig)
    pyplot.colorbar()
    
    req_N = S2 / (m_target**2 - m_vals**6) * ( 1 - 2*m_vals - 3*m_vals**2)
    
    # Handle bad values
    req_N = np.where((m_target**2 - m_vals**6)>0,req_N,1e99)
    
    ax.plot(m_vals,req_N,color='k')

    # Save the plot
    figname = ("projected_d." + file_format)
    if file_format == "eps":
        # Save in the paper location
        figname = join(paper_location,figname)
    pyplot.savefig(figname, format=file_format, bbox_inches="tight", pad_inches=0.05)
    
    fig.show()
    
    pyplot.show()
    
    return

if __name__ == "__main__":
    main(sys.argv)
