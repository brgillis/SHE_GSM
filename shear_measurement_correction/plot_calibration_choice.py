#!/usr/bin/env python

""" @file /disk2/brg/git/SHE_sim/shear_measurement_correction/plot_calibration_choice.py

    Created 12 Sep 2018

    @TODO: File docstring

    ---------------------------------------------------------------------

    Copyright (C) 2018 brg

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

from os.path import join
import sys

import matplotlib
from scipy.interpolate import InterpolatedUnivariateSpline as Spline

import matplotlib.pyplot as pyplot
import numpy as np
from shear_calibration import magic_values as mv


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

xlim = [-0.051, 0.051]
ylim = [-1, 50.]
m_true = 0.0075
m_sigma = 0.01

markersize = 12
pred_markersize = 64
fontsize = 24
ticksize = 16
file_format = "eps"
paper_location = "/home/brg/Dropbox/gillis-comp-shared/Papers/Shear_Bias/"


def gaus(x, m=0, s=1):
    return np.exp(-(x - m)**2 / (2 * s**2)) / np.sqrt(2 * np.pi * s**2)


def main(argv):
    """ @TODO main docstring
    """

    # Plot distribution of m_hat

    m_hat = np.linspace(xlim[0] * 1.1, xlim[1] * 1.1, 1000)
    dm_hat = m_hat[1:] - m_hat[:-1]
    p_m_hat = gaus(m_hat, m=m_true, s=m_sigma)

    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel(r"$\hat{m}$", fontsize=fontsize)
    ax.set_ylabel(r"$P(\hat{m})$", fontsize=fontsize)

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])

    #   Draw axes and limits
    ax.plot([-1, 1], [0, 0], label=None, color=(0, 0, 0), linestyle="solid", linewidth=0.5)
    ax.plot([0, 0], [-100, 100], label=None, color=(0, 0, 0), linestyle="solid", linewidth=0.5)
    ax.plot([-m_sigma, -m_sigma], [-100, 100], label=None, color=(0, 0, 0), linestyle="dashed", linewidth=1)
    ax.plot([m_sigma, m_sigma], [-100, 100], label=None, color=(0, 0, 0), linestyle="dashed", linewidth=1)

    #   Draw the histogram
    ax.plot(m_hat, p_m_hat, label=None, color=(0, 0, 0), linestyle="solid", linewidth=1)

    #   Shade in the regions
    ax.fill_between(m_hat[np.abs(m_hat) < m_sigma], p_m_hat[np.abs(m_hat) < m_sigma], color=(0.5, 0.5, 1.0))
    ax.fill_between(m_hat[m_hat <= -m_sigma], p_m_hat[m_hat <= -m_sigma], color=(1.0, 0.5, 0.5))
    ax.fill_between(m_hat[m_hat >= m_sigma], p_m_hat[m_hat >= m_sigma], color=(1.0, 0.5, 0.5))

    # Label the plot
    ax.text(0.06, 0.97, r"Measured $\hat{m}$", horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes,
            fontsize=24)

    # Save the plot
    figname = "calibration_choice_m_hat." + file_format
    if file_format == "eps" or file_format == "pdf":
        # Save in the paper location
        figname = join(paper_location, figname)
    pyplot.savefig(figname, format=file_format, bbox_inches="tight", pad_inches=0.05)

    fig.show()

    # Plot distribution of m'

    mp = (1 + m_true) * (1 - m_hat + 2 * m_hat**2) - 1
    dmp = mp[1:] - mp[:-1]
    p_mp = np.abs(p_m_hat[:-1] * dm_hat / dmp)

    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel(r"$m'$", fontsize=fontsize)
    ax.set_ylabel(r"$P(m')$", fontsize=fontsize)

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])

    #   Draw axes and limits
    ax.plot([-1, 1], [0, 0], label=None, color=(0, 0, 0), linestyle="solid", linewidth=0.5)
    ax.plot([0, 0], [-100, 100], label=None, color=(0, 0, 0), linestyle="solid", linewidth=0.5)
    ax.plot([-m_sigma, -m_sigma], [-100, 100], label=None, color=(0, 0, 0), linestyle="dashed", linewidth=1)
    ax.plot([m_sigma, m_sigma], [-100, 100], label=None, color=(0, 0, 0), linestyle="dashed", linewidth=1)

    #   Draw the histogram
    ax.plot(mp[:-1], p_mp, label=None, color=(0, 0, 0), linestyle="solid", linewidth=1)

    #   Shade in the regions
    ax.fill_between(mp[np.abs(m_hat) < m_sigma], p_mp[(np.abs(m_hat) < m_sigma)[:-1]], color=(0.5, 0.5, 1.0))
    ax.fill_between(mp[m_hat <= -m_sigma], p_mp[(m_hat <= -m_sigma)[:-1]], color=(1.0, 0.5, 0.5))
    ax.fill_between(mp[m_hat >= m_sigma][:-1], p_mp[(m_hat >= m_sigma)[:-1]], color=(1.0, 0.5, 0.5))

    # Label the plot
    ax.text(0.06, 0.97, r"Post-calibration $m'$", horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes,
            fontsize=24)

    # Save the plot
    figname = "calibration_choice_mp." + file_format
    if file_format == "eps" or file_format == "pdf":
        # Save in the paper location
        figname = join(paper_location, figname)
    pyplot.savefig(figname, format=file_format, bbox_inches="tight", pad_inches=0.05)

    fig.show()

    # Plot distribution of m^c

    mq = np.zeros_like(m_hat)
    mq[np.abs(m_hat) < m_sigma] = m_hat[np.abs(m_hat) < m_sigma]
    mq[np.abs(m_hat) >= m_sigma] = mp[np.abs(m_hat) >= m_sigma]

    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel(r"$m^c$", fontsize=fontsize)
    ax.set_ylabel(r"$P(m^c)$", fontsize=fontsize)

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])

    #   Draw axes and limits
    ax.plot([-1, 1], [0, 0], label=None, color=(0, 0, 0), linestyle="solid", linewidth=0.5)
    ax.plot([0, 0], [-100, 100], label=None, color=(0, 0, 0), linestyle="solid", linewidth=0.5)
    ax.plot([-m_sigma, -m_sigma], [-100, 100], label=None, color=(0, 0, 0), linestyle="dashed", linewidth=1)
    ax.plot([m_sigma, m_sigma], [-100, 100], label=None, color=(0, 0, 0), linestyle="dashed", linewidth=1)

    #   Calculations for the overlap region
    pmq_c_spline = Spline(mq[np.abs(m_hat) < m_sigma], p_mp[(np.abs(m_hat) < m_sigma)[:-1]])
    pmq_l_spline = Spline(np.flip(mq[m_hat <= -m_sigma], 0), np.abs(np.flip(p_mp[(m_hat <= -m_sigma)[:-1]], 0)))
    pmq_r_spline = Spline(np.flip(mq[m_hat >= m_sigma][:-1], 0), np.abs(np.flip(p_mp[(m_hat >= m_sigma)[:-1]], 0)))

    pmq_c_tot = np.where(np.logical_and(m_hat > mq[np.abs(m_hat) < m_sigma].min(),
                                        m_hat <= mq[np.abs(m_hat) < m_sigma].max()), pmq_c_spline(m_hat), 0).sum()
    pmq_c = np.where(np.abs(m_hat - m_true) < dm_hat[0], pmq_c_tot, 0)
    pmq_l = np.where(np.logical_and(m_hat > mq[m_hat <= -m_sigma].min(),
                                    m_hat <= mq[m_hat <= -m_sigma].max()), pmq_l_spline(m_hat), 0)
    pmq_r = np.where(np.logical_and(m_hat > mq[m_hat >= m_sigma].min(),
                                    m_hat <= mq[m_hat >= m_sigma].max()), pmq_r_spline(m_hat), 0)

    #   Draw the histogram
    ax.plot(m_hat, pmq_l + pmq_r + pmq_c, label=None, color=(0, 0, 0), linestyle="solid", linewidth=1)

    #   Shade in the regions
    ax.fill_between(m_hat, pmq_l + pmq_r + pmq_c, pmq_l + pmq_r, color=(0.5, 0.5, 1.0))
    ax.fill_between(m_hat, pmq_l, 0, color=(1.0, 0.5, 0.5))
    ax.fill_between(m_hat, pmq_r, 0, color=(1.0, 0.5, 0.5))

    # Label the plot
    ax.text(0.06, 0.97, r"Post-conditional-calibration $m^c$", horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes,
            fontsize=24)

    # Save the plot
    figname = "calibration_choice_mq." + file_format
    if file_format == "eps" or file_format == "pdf":
        # Save in the paper location
        figname = join(paper_location, figname)
    pyplot.savefig(figname, format=file_format, bbox_inches="tight", pad_inches=0.05)

    fig.show()

    pyplot.show()

    return


if __name__ == "__main__":
    main(sys.argv)
