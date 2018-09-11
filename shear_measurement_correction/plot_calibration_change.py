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

import cPickle
from os.path import join
import sys

import matplotlib

import matplotlib.pyplot as pyplot
import numpy as np
from shear_calibration import magic_values as mv


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


sigma_gm = mv.default_shape_sigma
sigma_gs = mv.default_shear_sigma

markersize = 12
pred_markersize = 64
fontsize = 24
ticksize = 16
file_format = "pdf"
paper_location = "/home/brg/Dropbox/gillis-comp-shared/Papers/Shear_Bias/"


def main(argv):
    """ @TODO main docstring
    """

    for calibration_set_size in ("1e4", "1e6"):
        for calibration_order in ("1st", "2nd"):
            # for calibration_order in ("bayesian",):

            if len(argv) > 1:
                calibration_set_size = argv[1]
            if len(argv) > 2:
                calibration_order = argv[2]

            if calibration_set_size == "1e6":
                size_label = r"$n = 10^6$"
                cal_scale = 0.1
            else:
                size_label = r"$n = 10^4$"
                cal_scale = 1.

            results_dir = "calibration_results_" + calibration_set_size + "_" + calibration_order

            # Load in the calculated results
            test_ms = [-0.2, -0.1, -0.01, -0.005, 0., 0.005, 0.01, 0.1, 0.2]
            c = 0.

            results = {}

            for m in test_ms:

                filename = join(results_dir, "calibration_results_m_" + str(m) + "_c_" + str(c) + ".bin")
                # filename = join(results_dir,"calibration_results_bayesian_m_" + str(m) + "_c_" + str(c) + ".bin")

                with open(filename, 'rb') as fi:
                    new_res = cPickle.load(fi)
                    new_res["m_mean"] = m
                    new_res["m_sigma"] = 0
                    results[(c, m)] = new_res

            # Set up the figure

            fig = pyplot.figure()
            fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

            ax = fig.add_subplot(1, 1, 1)
            ax.set_xlabel(r"Uncalibrated $m$", fontsize=fontsize)
            ax.set_ylabel(r"Post-calibration $m$", fontsize=fontsize)

            xlim = 0.31
            ylim = 0.12 * cal_scale

            ax.set_xlim(-xlim, xlim)
            ax.set_ylim(-ylim, ylim)

            xticks = [-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3]
            yticks = [-0.1 * cal_scale, -0.05 * cal_scale, 0, 0.05 * cal_scale, 0.1 * cal_scale]

            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks, fontsize=ticksize)

            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks, fontsize=ticksize)

            ax.text(0.06, 0.97, size_label, horizontalalignment='left',
                    verticalalignment='top', transform=ax.transAxes,
                    fontsize=24)

            # Draw x-axis
            ax.plot([-1, 1], [0, 0], label=None, color=(0, 0, 0), linestyle="solid", linewidth=0.5)

            # Plot the points

            for m in test_ms:
                r_color = (m + 0.2) / 0.4
                # b_color = (c + 0.1) / 0.2
                b_color = (m + 0.2) / 0.4

                if calibration_order == "2nd":
                    m2 = (4, 0, 45)
                else:
                    m2 = (3, 0, 0)

                # Create a dummy label for the center point, colored black
                if m == 0 and c == 0:
                    ax.errorbar([-10], [-10], marker=m2, markersize=markersize,
                                label=calibration_order + r"-order correction",
                                linestyle="None", linewidth=2, capthick=2,
                                markerfacecolor='None', markeredgecolor='k', color='k',
                                yerr=[1])

                color = (r_color, 1 - (r_color + b_color) / 2, b_color)

                res = results[(c, m)]

                ax.errorbar([res["mm_mean"]],
                            [res["mp_mean"]],
                            marker=m2, markerfacecolor='None', markeredgecolor=color,
                            markersize=markersize, label=None,
                            linestyle="None", linewidth=2,
                            capthick=2, capsize=4, color=color,
                            yerr=[res["mp_sigma"]])

            ax.legend(loc='lower right', fontsize=ticksize,
                      bbox_transform=ax.transAxes, numpoints=1, scatterpoints=1)

            # Save the plot
            figname = (calibration_set_size + "_"
                       + calibration_order + "_calibration." + file_format)
            if file_format == "eps" or file_format == "pdf":
                # Save in the paper location
                figname = join(paper_location, figname)
            pyplot.savefig(figname, format=file_format, bbox_inches="tight", pad_inches=0.05)

            fig.show()

    pyplot.show()
    return


if __name__ == "__main__":
    main(sys.argv)
