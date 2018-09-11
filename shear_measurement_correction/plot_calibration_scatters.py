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
    
markersize = 400
s_markersize = np.sqrt(markersize)
fontsize = 24
file_format = "eps"
paper_location = "/home/brg/Dropbox/gillis-comp-shared/Papers/Shear_Bias/"

m_min = -0.25
m_max = 0.25
m_points = 100

m_target = 0.002

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
        test_ms = [-0.2, -0.1, 0., 0.1, 0.2]
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
        ax.set_xlim(lim_factor*0.002,lim_factor*0.0032)
        ax.set_ylim(lim_factor*0.055,lim_factor*0.111)
        
        xticks = [lim_factor*0.002,lim_factor*0.0025,lim_factor*0.003]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)
        
        yticks = [lim_factor*0.06,lim_factor*0.07,lim_factor*0.08,lim_factor*0.09,lim_factor*0.10,lim_factor*0.11,]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)
        
        ax.text(0.05,0.95,size_label,horizontalalignment='left',
                verticalalignment='top',transform=ax.transAxes,
                fontsize=24)
    
        # Plot the points
        
        for m in test_ms:
            r_color = (m+0.2)/0.4
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
                    label_2p = "Predicted 2nd-order correction"
                    color='k'
                else:
                    label_0 = None
                    label_1 = None
                    label_2 = None
                    label_1p = None
                    label_2p = None
                    color = (r_color,1-(r_color+b_color)/2,b_color)
                
                
                res = results[(c,m)]
                
                ax.scatter([res["cm_sigma"]],[res["mm_sigma"]],
                            marker='.',s=markersize,label=label_0,color=color)
                
                ax.scatter([res["cp_sigma"]],[res["mp_sigma"]],
                            marker=(3, 0, 0),facecolor='None',edgecolor=color,
                            s=markersize,label=label_1,color=color)
                
                if not bayesian:
                    ax.scatter([res["cpp_sigma"]],[res["mpp_sigma"]],
                                marker=(4, 0, 45),facecolor='None',edgecolor=color,
                                s=markersize,label=label_2,color=color)
                
                    sigma_mm = 1/np.sqrt(float(calibration_set_size)) * ( sigma_gm/sigma_gs )
                    sigma_cm = sigma_mm*sigma_gs
                    
                    pred_cp_sigma = sigma_cm * (1 - m + m**2 + 1.5*sigma_mm**2)
                    pred_mp_sigma = sigma_mm * (1 - m - 2*m**2 + sigma_mm**2)
                    
                    pred_cpp_sigma = pred_cp_sigma * (1 - sigma_mm**2 + 2*m*sigma_mm**2 - m**3)
                    pred_mpp_sigma = sigma_mm * (1 - m + m**2 + sigma_mm**2)
                    
                    ax.scatter([pred_cp_sigma],[pred_mp_sigma],
                                marker=(3, 2, 0),s=markersize,label=label_1p,color=color)
                    
                    ax.scatter([pred_cpp_sigma],[pred_mpp_sigma],
                                marker=(4, 2, 45),s=markersize,label=label_2p,color=color)
                
                    # Draw lines connecting points for successive calibrations
                    ax.annotate("",xytext=(res["cp_sigma"],res["mp_sigma"]),
                                xy=((res["cp_sigma"]+res["cpp_sigma"])/2,(res["mp_sigma"]+res["mpp_sigma"])/2),
                                arrowprops=dict(arrowstyle="->",color=color))
                    ax.annotate("",xytext=((1.7*res["cp_sigma"]+res["cpp_sigma"])/2.7,(1.7*res["mp_sigma"]+res["mpp_sigma"])/2.7),
                                xy=(res["cpp_sigma"],res["mpp_sigma"]),
                                arrowprops=dict(arrowstyle="->",color=color))
                    
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
            ylim = (-0.010,0.012)
        else:
            size_label = r"$n = 10^4$"
            ylim = (-0.010,0.12)
        
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(r"$m$",fontsize=fontsize)
        ax.set_ylabel(r"$m'$, $m''$, $\sigma[m']$, and $\sigma[m'']$",fontsize=fontsize)
        
        ax.set_xlim(m_min,m_max)
        ax.set_ylim(*ylim)
        
        # ax.set_yscale("symlog",linthreshy=0.001)
        
        ax.text(0.5,0.95,size_label,horizontalalignment='center',
                verticalalignment='top',transform=ax.transAxes,
                fontsize=24)
            
        m_vals = np.linspace(m_min,m_max,m_points)
        
        S = mv.default_shape_sigma/mv.default_shear_sigma
        
        n = float(calibration_set_size)
        expec_mp_vals = S**2/n * (1 + m_vals) + m_vals**3
        expec_mpp_vals = S**2/n * (S**2/n + 3*m_vals**2) + m_vals**6
        sigma_mp_vals = S/np.sqrt(n) * (1 - m_vals - 2*m_vals**2 + S**2/n)
        sigma_mpp_vals = S/np.sqrt(n) * (1 - m_vals + m_vals**2 + S**2/n)
        
        test_mp_means = []
        test_mpp_means = []
        test_sigma_mp_means = []
        test_sigma_mpp_means = []
        
        for m in test_ms:
            test_mp_means.append(all_res[calibration_set_size][(0.,m)]["mp_mean"])
            test_mpp_means.append(all_res[calibration_set_size][(0.,m)]["mpp_mean"])
            test_sigma_mp_means.append(all_res[calibration_set_size][(0.,m)]["mp_sigma"])
            test_sigma_mpp_means.append(all_res[calibration_set_size][(0.,m)]["mpp_sigma"])
        
        ax.plot(m_vals,expec_mp_vals,color='r',linestyle="solid",linewidth=1,label=r"$\left<m'\right>$")
        ax.plot(m_vals,expec_mpp_vals,color='r',linestyle="solid",linewidth=2,label=r"$\left<m''\right>$")
        ax.plot(test_ms,test_mp_means,markeredgecolor='r',linestyle='none',label=r"$\overline{\hat{m}'}$",
                   marker=(3, 0, 0),markersize=s_markersize,markerfacecolor='None')
        ax.plot(test_ms,test_mpp_means,markeredgecolor='r',linestyle='none',label=r"$\overline{\hat{m}''}$",
                   marker=(4, 0, 45),markersize=s_markersize,markerfacecolor='None')
        
        ax.plot(m_vals,sigma_mp_vals,color ='b',linestyle="dashed",linewidth=1,label=r"$\left<\sigma\left[m'\right]\right>$")
        ax.plot(m_vals,sigma_mpp_vals,color ='b',linestyle="dashed",linewidth=2,label=r"$\left<\sigma\left[m''\right]\right>$")
        ax.plot(test_ms,test_sigma_mp_means,markeredgecolor='b',linestyle='none',label=r"$\sigma\left[\hat{m}'\right]$",
                   marker=(3, 0, 0),markersize=s_markersize,markerfacecolor='None')
        ax.plot(test_ms,test_sigma_mpp_means,markeredgecolor='b',linestyle='none',label=r"$\sigma\left[\hat{m}''\right]$",
                   marker=(4, 0, 45),markersize=s_markersize,markerfacecolor='None')
        
        # Draw the zero line
        ax.plot(m_vals,np.zeros_like(m_vals),color='k',linestyle="solid")
        
        ax.legend(loc="center left",numpoints=1)
        
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
