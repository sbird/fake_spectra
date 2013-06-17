# -*- coding: utf-8 -*-
"""Make some plots of 2d histograms, velocity widths vs some other things"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import plot_spectra as ps
import os.path as path
import numpy as np
from save_figure import save_figure

base="/home/spb/scratch/Cosmo/"
outdir = base + "plots/mass/"
print "Plots at: ",outdir


def plot_mass_hists(sim, snap, ff=False,color="red"):
    """Plot mass histogram"""
    halo = "Cosmo"+str(sim)+"_V6"
    out = "cosmo"+str(sim)+"_mass_z"+str(snap)
    if ff:
        halo+="_512"
        out+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo)
    (mbins, pdf) = hspec.mass_hist()
    plt.semilogx(mbins,pdf,color=color)

def plot_mass_vs(sim, snap, ff=False):
    """Plot mass vs metallicity and vel width"""
    halo = "Cosmo"+str(sim)+"_V6"
    out = "cosmo"+str(sim)+"_met_mass_z"+str(snap)
    if ff:
        halo+="_512"
        out+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo)
    hspec.plot_Z_vs_mass()
    save_figure(path.join(outdir,out))
    plt.clf()
    out = "cosmo"+str(sim)+"_vel_mass_z"+str(snap)
    if ff:
        out+="_512"
    hspec.plot_vel_vs_mass("Si",2)
    save_figure(path.join(outdir,out))
    plt.clf()

if __name__ == "__main__":
    colors=["blue", "purple", "orange", "red"]
#     for ss in (0,2,3):
#         for zz in (54,60,68):
#             plot_mass_vs(ss, zz)

    for zz in (54,60,68):
        plot_mass_vs(0, zz,True)

    for ss in (0,2,3):
        plot_mass_hists(ss, 60, ff=False,color=colors[ss])
    save_figure(path.join(outdir,"cosmo_halos_feedback_z3"))
    plt.clf()

    for zz in (54,60,68):
        plot_mass_hists(0,zz,ff=True,color=colors.pop())
    save_figure(path.join(outdir,"cosmo_halos_0_512_zz"))
    plt.clf()
