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
import numpy as np

import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from shear_calibration import magic_values as mv

sigma_gm = mv.default_shape_sigma
sigma_gs = mv.default_shear_sigma
    
markersize = 100
fontsize = 24
file_format = "eps"
paper_location = "/home/brg/Dropbox/gillis-comp-shared/Papers/Shear_Bias/"

m_min = -0.15
m_max = 0.15
m_points = 100

bayesian = False

def main(argv):
    """ @TODO main docstring
    """
    
    all_res = {}
    
    for calibration_set_size in ("1e4","1e6"):
    
        if len(argv) > 1:
            calibration_set_size = argv[1]
            
        if calibration_set_size == "1e6":
            size_label = r"$n = 10^6$"
        else:
            size_label = r"$n = 10^4$"
            
        results_dir_root = "calibration_results_" + calibration_set_size
        
        results_dir_1 = results_dir_root + "_1st"
        results_dir_2 = results_dir_root + "_2nd"
        results_dir_b = results_dir_root + "_bayesian"
        
        # Load in the calculated results
        test_ms = [-0.1, 0., 0.1]
        test_cs = [-0.1, 0., 0.1]
        
        results = {}
        
        for m in test_ms:
            for c in test_cs:
                
                if bayesian:
                
                    filename_1 = join(results_dir_b,"calibration_results_bayesian_m_" + str(m) + "_c_" + str(c) + ".bin")
                    
                    with open(filename_1,'rb') as fi:
                        new_res = cPickle.load(fi)
                    new_res["c_mean"] = c
                    new_res["m_mean"] = m
                    new_res["c_sigma"] = 0
                    new_res["m_sigma"] = 0
                    
                else:
                
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
                
        all_res[calibration_set_size] = results
                    
        # We'll plot both the change in measured values and actual values
            
        # Set up the figure
        
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
        ax.set_ylim(lim_factor*0.07,lim_factor*0.100)
        
        xticks = [lim_factor*0.002,lim_factor*0.0025,lim_factor*0.003,]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)
        
        yticks = [lim_factor*0.07,lim_factor*0.08,lim_factor*0.09,lim_factor*0.10,]
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
                    label_2 = "2nd-order correction"
                    if bayesian:
                        label_1 = "Bayesian correction"
                    else:
                        label_1 = "1st-order correction"
                    label_1p = "Predicted 1st-order correction"
                    color='k'
                else:
                    label_0 = None
                    label_1 = None
                    label_2 = None
                    label_1p = None
                    color = (r_color,1-(r_color+b_color)/2,b_color)
                
                
                res = results[(c,m)]
                
                ax.scatter([res["cm_sigma"]],[res["mm_sigma"]],
                            marker='.',s=markersize,label=label_0,color=color)
                
                ax.scatter([res["cp_sigma"]],[res["mp_sigma"]],
                            marker='o',facecolor='None',edgecolor=color,
                            s=markersize,label=label_1,color=color)
                
                if not bayesian:
                    ax.scatter([res["cpp_sigma"]],[res["mpp_sigma"]],
                                marker='s',facecolor='None',edgecolor=color,
                                s=markersize,label=label_2,color=color)
                
                    sigma_mm = 1/np.sqrt(float(calibration_set_size)) * ( sigma_gm/sigma_gs )
                    sigma_cm = sigma_mm*sigma_gs
                    
                    pred_cp_sigma = sigma_cm * (1 - m + m**2 + 1.5*sigma_mm**2)
                    pred_mp_sigma = sigma_mm * (1 - m - 2*m**2 + sigma_mm**2)
                    
                    ax.scatter([pred_cp_sigma],[pred_mp_sigma],
                                marker='x',s=markersize,label=label_1p,color=color)
                
                    # Draw lines connecting points for successive calibrations
                    ax.annotate("",xytext=(res["cp_sigma"],res["mp_sigma"]),xy=(res["cpp_sigma"],res["mpp_sigma"]),
                                arrowprops=dict(arrowstyle="->",color=color))
                    
                # Also draw half-arrows for the first correction
                ax.annotate("",xytext=(res["cm_sigma"],res["mm_sigma"]),
                            xy=((res["cm_sigma"]+res["cp_sigma"])/2,(res["mm_sigma"]+res["mp_sigma"])/2),
                            arrowprops=dict(arrowstyle="->",color=color))
                ax.annotate("",xytext=((1.15*res["cm_sigma"]+res["cp_sigma"])/2.15,(1.15*res["mm_sigma"]+res["mp_sigma"])/2.15),
                            xy=(res["cp_sigma"],res["mp_sigma"]),
                            arrowprops=dict(arrowstyle="->",color=color))
                
        ax.legend(loc='lower right',scatterpoints=1)
    
        # Save the plot
        figname = ("both_" + calibration_set_size + "_scatters." + file_format)
        if file_format == "eps":
            # Save in the paper location
            figname = join(paper_location,figname)
        pyplot.savefig(figname, format=file_format, bbox_inches="tight", pad_inches=0.05)
        
        fig.show()
        
        if not bayesian:
        
            # Save a table of the data as well
            tabname = join(paper_location,"table_data_" + calibration_set_size + ".tex")
            
            with open(tabname,'w') as fo:
                
                # Print the header portion
                line = ("\\begin{center}\n" +
                         "\\begin{tabular}{rrrrrrrrrrrrrr}\n" +
                         "$m$ & $c$ & $\\overline{\hat{m}}$ & $\\overline{m'}$ & $\\overline{m''}$ & $\\std{\hat{m}}$" + 
                         "& $\\std{m'}$ & $\\std{m''}$ & $\\overline{\hat{c}}$ & $\\overline{c'}$ & $\\overline{c''}$ & $\\std{\hat{c}}$ & $\\std{c'}$ & $\\std{c''}$ \\\\ \\hline\n")
                
                fo.write(line)
                
                for m in test_ms:
                    for c in test_cs:
                        res = results[(c,m)]
                        line = (" %1.1f & %1.1f & %1.4f & %1.4f & %1.4f & %1.4f & %1.4f & %1.4f & %1.4f & %1.4f & %1.4f & %1.4f & %1.4f & %1.4f" %
                                (m, c, res["mm_mean"], res["mp_mean"], res["mpp_mean"], res["mm_sigma"], res["mp_sigma"], res["mpp_sigma"],
                                res["cm_mean"], res["cp_mean"], res["cpp_mean"], res["cm_sigma"], res["cp_sigma"], res["cpp_sigma"]) )
            
                        # Unless it's the last one, add a linebreak
                        if (m,c) is not (test_ms[-1],test_cs[-1]):
                            line += " \\\\"
                            
                        line += "\n"
                        
                        fo.write(line)
                            
                # Print the footer portion
                fo.write("\\end{tabular}\n" +
                         "\\end{center}\n")
            
    # Now set up a figure of sigma(m) versus m
    
    for calibration_set_size in ("1e4", "1e6"):
        
        if calibration_set_size == "1e6":
            size_label = r"$n = 10^6$"
        else:
            size_label = r"$n = 10^4$"
        
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(r"$m$",fontsize=fontsize)
        ax.set_ylabel(r"$m$ or $\sigma[m]$",fontsize=fontsize)
        
        ax.set_xlim(-0.15,0.15)
        
        # ax.set_yscale("symlog",linthreshy=0.001)
        
        ax.text(0.05,0.85,size_label,horizontalalignment='left',
                verticalalignment='top',transform=ax.transAxes,
                fontsize=24)
            
        m_vals = np.linspace(m_min,m_max,m_points)
        
        S = mv.default_shape_sigma/mv.default_shear_sigma
        
        n = float(calibration_set_size)
        expec_mp_vals = S**2/n * (1 + m_vals) + m_vals**3
        sigma_mp_vals = S/np.sqrt(n) * (1 - m_vals - 2*m_vals**2 + S**2/n)
        
        ax.plot(m_vals,expec_mp_vals,color='r',linestyle="solid",label=r"$\left<m'\right>$")
        ax.plot(m_vals,sigma_mp_vals,color='b',linestyle="dashed",label=r"$\left<\sigma\left[m'\right]\right>$")
        ax.plot(m_vals,np.zeros_like(m_vals),color='k',linestyle="solid")
        
        test_mp_means = []
        test_sigma_m_means = []
        
        for m in test_ms:
            test_mp_means.append(all_res[calibration_set_size][(0.,m)]["mp_mean"])
            test_sigma_m_means.append(all_res[calibration_set_size][(0.,m)]["mp_sigma"])
            
        ax.scatter(test_ms,test_mp_means,color='r',label=r"$\overline{\hat{m}'}$",
                   marker='o',s=markersize,facecolor='None')
        ax.scatter(test_ms,test_sigma_m_means,color='b',label=r"$\sigma\left[\hat{m}'\right]$",
                   marker='o',s=markersize,facecolor='None')
        
        ax.legend(loc="center left",scatterpoints=1)
        
        # Save the plot
        figname = ("both_" + calibration_set_size + "_proj_m_scatters." + file_format)
        if file_format == "eps":
            # Save in the paper location
            figname = join(paper_location,figname)
        pyplot.savefig(figname, format=file_format, bbox_inches="tight", pad_inches=0.05)
        
        fig.show()
            
    pyplot.show()
    
    return

if __name__ == "__main__":
    main(sys.argv)
